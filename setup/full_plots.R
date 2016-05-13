library(ggplot2)
source('setup/trials.R')
source('setup/read_variants.R')
figs <- 'figs_Full_run'
db <- src_sqlite('Full_run.db')
db_res <- src_sqlite('Full_run_res.db')

dir.create(figs)

tmp_func <- function(x){
  if(x$area[1] %in% c('WI','EG')){
    tmp <- x
    tmp$area <- 'EG+WI'
    bind_rows(x,tmp)
  } else{
    x
  }
}


## ogives 
ogive <- function(a50,d){
  a <- seq(0,50,by=0.001)
  a50 <- d
  a95 <- d
  for(i in length(d)){
    tmp <- 1/(1 + exp(-(a-a50[i])/d[i]))
    a50[i] <- max(a[which(tmp<0.5)])
    a95[i] <- max(a[which(tmp<0.95)])
  }
  return(list(a50=a50,a95=a95))
}

tmp2 <- 
  tbl(db,'naf_all') %>%
  collect() %>% 
  filter(!grepl('NF-S',trial),grepl('disp',constant)|constant %in% c('rmsig','rfsig','rm50','rf50')) %>%
  spread(constant,value) %>% 
  group_by(trial) %>% 
  mutate(sel.male.50 = ogive(rmsig,rm50)$a50,
         sel.female.50 = ogive(rfsig,rf50)$a95,
         sel.male.95 = ogive(rmsig,rm50)$a50,
         sel.female.95 = ogive(rfsig,rf50)$a95) %>% 
  select(-c(rf50:rmsig)) %>% 
  gather(constant,value,-c(ref,trial)) %>%
  filter(!(constant %in% c('sel.female.50','sel.female.95'))| 
           (constant %in% c('sel.female.50','sel.female.95') & value < 10))
#  mutate(constant = plyr::revalue(constant,
#                                  c('disp1'='C1C2','disp2'='C2C3',
#                                  'sel.male.50'='Male a50','sel.female.50'='Female a50',
#                                  'sel.male.95'='Male a95','sel.female.95'='Female a95'))) 
  
  for(cons in unique(tmp2$constant)){
    pdf(file=sprintf('%s/%s.pdf',figs,cons),width=11,height = 11)
    print(ggplot(filter(tmp2,constant==cons),aes(value)) + geom_histogram() + 
      facet_wrap(~trial,scale='free') +
      geom_vline(data=filter(tmp2,ref=='0',constant==cons),aes(xintercept=value))+
      theme_bw() + ggtitle(cons) + ylab('Num trials') +
      xlab('Value'))
    dev.off()
  }

sight <- tbl(db,'naf_sight') %>%
  mutate(upper = exp(log(obs)+1.96*cv),
         lower = exp(log(obs)-1.96*cv)) %>%
  select(year,area,obs,lower,upper) %>%
  collect() %>%
  mutate(area=toupper(area),
         area=ifelse(area=='EI','EI/F',area))  %>%
  split(.$area) %>%
  map(~tmp_func(.)) %>% 
  bind_rows() %>%
  group_by(year,area) %>% 
  summarise(obs=sum(obs),lower=sum(lower),upper=sum(upper))


pop <- 
  tbl(db,'naf_pop') %>% 
  filter(pop_type == '1+') %>%
  collect(n=Inf) %>%
  group_by(year,ref,hypo,msyr,pop_id) %>%
  summarise(med = median(number),
            cond.975 = quantile(number,0.975),
            cond.0.25 = quantile(number,0.025),
            number = sum(ifelse(trial==0,number,0))) %>%
  ungroup() %>%
  mutate(msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref)) %>%
  rename(area=pop_id) %>%
  left_join(sight) %>%
  mutate(area = ordered(area,levels=c('EC','WG','EG','WI','EG+WI','EI/F','N','SP')),
         type = plyr::revalue(trialtype, 
                       c('A'='Pro rate abundance','B'='Baseline','C'='CPUE',
                         'F'='C2 to EG in 1985 - 2025 (opt. b)',
                         'G'='C2 to EG in 1985 (opt. a)',
                         'H'='High historical catch series',
                         'T'='Tag loss 20% yr 1, 10%/yr thereafter',
                         'J'='g(0) = 0.8',
                         'U'='Dome shaped selectivity',
                         'Y'='8 year future survey interval',
                         'D'='Upper bound on dispersal')))


pop.fem <- 
  tbl(db,'naf_pop') %>% 
  filter(pop_type == 'Mature females') %>%
  collect(n=Inf) %>%
  group_by(year,ref,hypo,msyr,pop_id) %>%
  summarise(med = median(number),
            cond.975 = quantile(number,0.975),
            cond.0.25 = quantile(number,0.025),
            number = sum(ifelse(trial==0,number,0))) %>%
  ungroup() %>%
  mutate(msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref)) %>%
  rename(pop=pop_id) %>%
  mutate(pop = ordered(pop,levels=c('W','C1','C2','C3','E','S'))) 

wi_tag <- 
  tbl(db,'naf_tag') %>%
  collect(n=Inf) %>%
  mutate(area = toupper(area)) %>% 
  filter(area %in% c('EC','WI','EGWI','EG+WI')) %>%
  mutate(msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         type = plyr::revalue(trialtype, 
                              c('A'='Pro rate abundance','B'='Baseline','C'='CPUE',
                                'F'='C2 to EG in 1985 - 2025 (opt. b)',
                                'G'='C2 to EG in 1985 (opt. a)',
                                'H'='High historical catch series',
                                'T'='Tag loss 20% yr 1, 10%/yr thereafter',
                                'J'='g(0) = 0.8',
                                'U'='Dome shaped selectivity'))) %>%
  group_by(ref,rela,area) %>%
  arrange(year) %>%
  mutate(prd1 = cumsum(prd),
         obs1 = cumsum(obs))



for(ref1 in unique(gsub('(..-..).+','\\1',pop$ref))){
  pdf(file=sprintf('%s/%s.pop.pdf',figs,ref1),width=7,height = 8)
  print(pop %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
          filter(ref %in% paste(ref1,c(1,4),sep='-')) %>%
          ggplot(aes(year,number,lty=msyr,fill=msyr)) + 
          geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.5) + 
          geom_line() + 
          geom_line(aes(year,med),col='red') + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,obs),col='black') + 
          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1) +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2015,by=20),
                             minor_breaks = seq(1860,2015,by=5))+
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle(ref1))
  dev.off()
}


for(ref1 in unique(gsub('(..-..).+','\\1',pop$ref))){
  pdf(file=sprintf('%s/%s.popfem.pdf',figs,ref1),width=7,height = 8)
  print(pop.fem %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
          filter(ref %in% paste(ref1,c(1,4),sep='-')) %>%
          ggplot(aes(year,number,lty=msyr,fill=msyr)) + 
          geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.5) + 
          geom_line() + 
          geom_line(aes(year,med),col='red') + 
          facet_wrap(~pop,scale='free_y',ncol=2) +
#          geom_point(aes(year,obs),col='black') + 
#          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('Mature Females') + xlab('Year') +
          geom_text(aes(label=pop),x=1960,y=Inf, vjust = 2,hjust = 1) +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2015,by=20),
                             minor_breaks = seq(1860,2015,by=5))+
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle(ref1))
  dev.off()
}





for(ms in c(1,4)){
pdf(file=sprintf('%s/baseline.%s.pop.pdf',figs,ms),width=11,height = 8)
print(pop %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
        filter(grepl(sprintf('NF-B.-%s',ms),ref)) %>%
        ggplot(aes(year,number,lty=as.factor(hypo))) + 
        #          geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.5) + 
        geom_line() + 
        #          geom_line(aes(year,med),col='red') + 
        facet_wrap(~area,scale='free_y',ncol=2) +
        geom_point(aes(year,obs),col='black') + 
        geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
        theme_bw() + ylab('1+ population') + xlab('Year') +
        geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1,col='black') +
        theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
              axis.text.x=element_text(size = 7),
              legend.position = c(0.6,0.1),
              panel.margin = unit(0.2,'cm'),
              plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
              strip.background = element_blank(),
              strip.text.x = element_blank())+
        scale_x_continuous(breaks = seq(1860,2015,by=20),
                           minor_breaks = seq(1860,2015,by=5))+
        expand_limits(y = 0)+
        #scale_fill_manual(values = c( "darkred", "darkgreen"))+
        #scale_fill_hue(l=40)+  
        scale_linetype_manual(name = NULL,
                             values = as.factor(c(1,2,3,5,6)),
                             labels = paste('Hypothesis',c(1,2,3,5,6))) +
          ggtitle(sprintf('Baseline %s%% (%s)',ms,ifelse(ms==1,'1+','Mature'))))
  dev.off()
}



age <- 
  tbl(db,'naf_age') %>%
  filter(trial ==0) %>%
  collect() %>% 
  distinct(year,area,age,ref,.keep_all = TRUE) %>% 
  gather(var,num, -c(year:age,trial,ref)) %>%
  separate(var,c('type','sex')) %>% 
  #filter(type!='obs') %>%
  spread(type,num) %>%
  mutate(sex = ifelse(sex=='m','Males','Females'),
         msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         type = plyr::revalue(trialtype, 
                              c('A'='Pro rate abundance','B'='Baseline','C'='CPUE',
                                'F'='C2 to EG in 1985 - 2025 (opt. b)',
                                'G'='C2 to EG in 1985 (opt. a)',
                                'H'='High historical catch series',
                                'T'='Tag loss 20% yr 1, 10%/yr thereafter',
                                'J'='g(0) = 0.8',
                                'U'='Dome shaped selectivity'))) %>%
#  left_join(age.dat) %>%
#  rename(prd= num) %>%
  group_by(ref,year,sex) %>%
  mutate(obs=obs/sum(obs,na.rm=TRUE),
         prd = prd/sum(prd,na.rm=TRUE),
         diff = obs-prd,
         dir = as.character(sign(diff+1e-10)))


for(ref1 in unique(gsub('(..-..).+','\\1',age$ref))){
  pdf(file=sprintf('%s/%s.age.pdf',figs,ref1),width=5,height = 6)
  print( 
    age %>% 
      filter(year<1990,ref %in% paste(ref1,c(1,4),sep='-')) %>%
      ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
      facet_grid(msyr~sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
      theme(legend.position='none') + ggtitle(ref1) +
      xlab('Year') + ylab('Age (in years)')
  )
  dev.off()
}


pdf(file=sprintf('%s/baseline1-4.age.pdf',figs),width=11,height = 8)
print( 
  
age %>% 
  filter(year<1990,grepl('NF-B[0-4]-.',ref)) %>%
  ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
  facet_grid(msyr~hypo+sex) + theme_minimal() + 
  scale_color_manual(values = c('lightblue','black')) + 
  theme(legend.position='none',
        axis.text.x = element_text(angle=90,size = 9,vjust = .3)) + 
  scale_x_continuous(breaks = seq(1970,1990,by=2))+
  ggtitle('Catch at age diagnostic for baseline hypotheses (1-4)') + 
  xlab('Year') + ylab('Age (in years)') + scale_size(range=c(0.5,3))
)
dev.off()

pdf(file=sprintf('%s/baseline5-8.age.pdf',figs),width=11,height = 8)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-B[5-8]-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~hypo+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for baseline hypotheses (5-8)') + 
    xlab('Year') + ylab('Age (in years)') 
)
dev.off()



pdf(file=sprintf('%s/other2A-G.age.pdf',figs),width=11,height = 8)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[(^B)|(^J-U)]2-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity  (A-G), hypotheses 2') + 
    xlab('Year') + ylab('Age (in years)') 
)
dev.off()

pdf(file=sprintf('%s/other2H-U.age.pdf',figs),width=11,height = 8)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[(^A-G)]2-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity trials, hypotheses 2') + 
    xlab('Year') + ylab('Age (in years)') 
)
dev.off()


pdf(file=sprintf('%s/other3A-G.age.pdf',figs),width=11,height = 8)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[(^B)|(^J-U)]3-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity  (A-G), hypotheses 2') + 
    xlab('Year') + ylab('Age (in years)') 
)
dev.off()

pdf(file=sprintf('%s/other3H-U.age.pdf',figs),width=11,height = 8)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[(^A-G)]3-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity trials, hypotheses 2') + 
    xlab('Year') + ylab('Age (in years)') 
)
dev.off()




# #for(ref1 in unique(gsub('(..-..).+','\\1',age$ref))){
# #  pdf(file=sprintf('%s/%s.tag.pdf',figs,ref1),width=5,height = 6)
# #  ref.tmp <- paste(ref1,c(1,4),sep='-')
# #  print( 
#  wi_tag %>%
#       filter(trialtype=='B',hypo == 8) %>%
#       ggplot(aes(year,ifelse(obs>0,obs1,NA),lty=msyr)) + geom_point() + 
#       geom_line(aes(year,prd1)) + 
#       facet_grid(rela~area,scale='free_y')+
#       theme_bw() + xlim(c(1965,2010)) + ylab('Num. tags recaptured in WI') +
#       xlab('Year') + ggtitle(ref1)
#   ) 
#   dev.off()
# }
# 

pdf(file=sprintf('%s/baseline.tag.pdf',figs),width=11,height = 8)
print(
  wi_tag %>% 
  filter(grepl('NF-B.-.',ref)) %>%
    ungroup() %>% 
  mutate(area = toupper(area),
         rela = toupper(rela)) %>% 
  filter(area == 'WI',rela =='WI') %>%
  ggplot(aes(year,ifelse(obs>0,obs1,NA),lty=msyr)) + geom_point() + 
  geom_line(aes(year,prd1)) + 
  facet_wrap(~hypo,scale='free_y')+
  theme_minimal() + xlim(c(1965,2010)) + ylab('Num. tags recaptured in WI') +
    guides(lty=guide_legend(title='MSYR'))+
  theme(legend.position = c(0.8,0.2), 
        axis.text.x = element_text(angle=45)) +
  xlab('Year') + ggtitle('Baseline trials tag recaptures within WI')
)
dev.off()



pdf(file=sprintf('%s/other.tag.pdf',figs),width=11,height = 8)
print(
  wi_tag %>% 
    filter(grepl('NF-[^B|^S].-.',ref)) %>%
    ungroup() %>% 
    mutate(area = toupper(area),
           rela = toupper(rela)) %>% 
    filter(area == 'WI',rela =='WI') %>%
    #filter(trialtype=='B',hypo == 8) %>%
    ggplot(aes(year,ifelse(obs>0,obs1,NA),lty=msyr)) + geom_point() + 
    geom_line(aes(year,prd1)) + 
    facet_wrap(~type+hypo,scale='free_y')+
    theme_minimal() + xlim(c(1965,2010)) + ylab('Num. tags recaptured in WI') +
    guides(lty=guide_legend(title='MSYR'))+
    theme(legend.position = c(0.6,0.1), 
          axis.text.x = element_text(angle=45)) +
    xlab('Year') + ggtitle('Other trials tag recaptures')
)
dev.off()




## first 10 trajectories
pdf(file=sprintf('%s/first10.pdf',figs),width=11,height = 8)
print(tbl(db,'naf_pop') %>%
        collect(n=Inf) %>% 
        filter(trial < 10) %>% 
        collect() %>% 
        ggplot(aes(year,number,group=interaction(as.factor(msyr),as.factor(trial),as.factor(hypo)),
                   col=as.factor(hypo),lty=as.factor(msyr))) + geom_line() + 
        facet_wrap(~pop_id) + theme_bw()
)
dev.off()


#system(sprintf('pdfjoin %s/NF*.pop.pdf',figs))
#system(sprintf('pdfjoin %s/NF*.age.pdf',figs))
#system(sprintf('pdfjoin %s/NF*.tag.pdf',figs))

#system('pdfjoin rep.pdf figs/baseline.*.pop.pdf figs_final/baseline.*.pop.pdf figs/baseline.tag.pdf figs_final/baseline1-4.age.pdf  figs_final/baseline5-8.age.pdf legends.pdf NF-W1.pop-joined.pdf NF-W1.age-joined.pdf NF-W1.tag-joined.pdf')

#system('pdfjoin figs/baseline.*.pop.pdf figs_final/baseline.*.pop.pdf figs/baseline.tag.pdf figs_final/baseline1-4.age.pdf  figs_final/baseline5-8.age.pdf')

tmp <- 
  tbl(db,'naf_all')  %>% 
  filter(constant == 'disp1',trial=='NF-B1-1',value>0.2)  %>% 
  select(trial = ref) %>% 
  collect(n=Inf)


pdf(file=sprintf('%s/cpue.pdf',figs),width=11,height = 8)
print(
  tbl(db,'naf_cpe') %>% 
    collect() %>% 
    ggplot(aes(year,obs,lty=ref)) + geom_point() + geom_line(aes(y=prd)) +
    facet_wrap(~series,scale='free') + 
    theme_minimal() + theme(legend.position ='none') +
    ylab('CPUE') + xlab('Year') +
    ggtitle('Fit to CPUE series')
)
dev.off()


pdf(file=sprintf('%s/cpue.xy.pdf',figs),width=11,height = 8)
print(
  tbl(db,'naf_cpe') %>% 
    collect() %>% 
    ggplot(aes(prd,obs,col=ref)) + geom_text(aes(label=year)) + geom_abline(slope=1,intercept=0) +
    facet_wrap(~series,scale='free') + 
    theme_minimal() + #theme(legend.position ='none') +
    ylab('CPUE') + xlab('Year') +
    ggtitle('Fit to CPUE series')
)
dev.off()



for(ms in c(1,4)){
  for(hyp in 2:3){
  pdf(file=sprintf('%s/other.%s.%s.pop.pdf',figs,hyp,ms),width=11,height = 8)
  print(pop %>% 
          filter(!(trialtype %in% c('B','S')), msyr == ms/100,hypo == hyp) %>%
          ggplot(aes(year,ifelse(trialtype=='J',0.8*number,number),col=type)) + 
          #          geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.5) + 
          geom_line() + 
          #          geom_line(aes(year,med),col='red') + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,obs),col='black') + 
          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                legend.position = c(0.75,0.1),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2015,by=20),
                             minor_breaks = seq(1860,2015,by=5))+
          expand_limits(y = 0)+
          guides(col=guide_legend(title=NULL,ncol=2))+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+  
#          scale_linetype_manual(name = NULL,
#                                values = as.factor(1:8),
#                                labels = paste('Hypothesis',1:8)) +
          ggtitle(sprintf('Other hypothesis %s trials %s%% (%s)',
                          hyp,ms,ifelse(ms==1,'1+','Mature'))))
  dev.off()
}
}

res.by.year <- 
  tbl(db_res,'naf_resbyyear') %>% 
  collect(n=Inf) %>% 
  filter(pop_type == 'area') %>%
  collect(n=Inf) %>%
  mutate(msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         type = plyr::revalue(trialtype, 
                              c('A'='Pro rate abundance','B'='Baseline','C'='CPUE',
                                'F'='C2 to EG in 1985 - 2025 (opt. b)',
                                'G'='C2 to EG in 1985 (opt. a)',
                                'H'='High historical catch series',
                                'T'='Tag loss 20% yr 1, 10%/yr thereafter',
                                'J'='g(0) = 0.8'))) %>%
  rename(area=pop_id) %>%
  left_join(sight) %>%
  mutate(area = ordered(area,levels=c('EC','WG','EG','WI','EG+WI','EI/F','N','SP')))


for(ref1 in unique(res.by.year$ref)){
  pdf(file=sprintf('%s/%s.res.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year %>% 
          filter(ref==ref1) %>% 
          ggplot(aes(year,abundance_med,lty=variant,fill=variant)) + 
          geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.5) + 
          geom_line() + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,obs),col='black') + 
          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1) +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=20),
                             minor_breaks = seq(1860,2115,by=5))+
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle(ref1))
  dev.off()
}



for(ref1 in unique(res.by.year$ref)){
  pdf(file=sprintf('%s/%s.res.pop.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year %>% 
          filter(ref==ref1) %>% 
          ggplot(aes(year,abundance_med,lty=variant,fill=variant)) + 
          geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          geom_line() + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,obs),col='black') + 
          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1) +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=20),
                             minor_breaks = seq(1860,2115,by=5))+
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle(ref1))
  dev.off()
}


for(ref1 in unique(res.by.year$ref)){
  pdf(file=sprintf('%s/%s.res.catbyar.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year %>% 
          filter(ref==ref1) %>% 
          ggplot(aes(year,catch_med,lty=variant,fill=variant)) + 
          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.2) + 
          geom_line() + 
          facet_wrap(~area,scale='free_y',ncol=2) +
#          geom_point(aes(year,obs),col='black') + 
#          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1) +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=20),
                             minor_breaks = seq(1860,2115,by=5))+
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle(ref1))
  dev.off()
}


thresh <- 
  tbl(db_res,'man_thresh') %>% 
  collect(n=Inf) %>% 
  rename(trialtype=type)

res.by.year.fem <- 
  tbl(db_res,'naf_resbyyear') %>% 
  collect(n=Inf) %>% 
  filter(pop_type != 'area') %>%
  mutate(msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         type = plyr::revalue(trialtype, 
                              c('A'='Pro rate abundance','B'='Baseline','C'='CPUE',
                                'F'='C2 to EG in 1985 - 2025 (opt. b)',
                                'G'='C2 to EG in 1985 (opt. a)',
                                'H'='High historical catch series',
                                'T'='Tag loss 20% yr 1, 10%/yr thereafter',
                                'J'='g(0) = 0.8'))) %>%
  rename(pop=pop_id) %>%
  mutate(pop = ordered(pop,levels=c('W','C1','C2','C3','E','S')))



for(ref1 in unique(res.by.year$ref)){
  pdf(file=sprintf('%s/%s.res.fem.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year.fem %>% 
          filter(ref==ref1) %>% 
          ggplot(aes(year,abundance_med,lty=variant,fill=variant)) + 
          geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.5) + 
          geom_line() + 
          facet_wrap(~pop,scale='free_y',ncol=2) +
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=pop),x=1960,y=Inf, vjust = 2,hjust = 1) +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=20),
                             minor_breaks = seq(1860,2115,by=5))+
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle(ref1))
  dev.off()
}




for(ref1 in unique(res.by.year$ref)){
  pdf(file=sprintf('%s/%s.res.femcat.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year.fem %>% 
          filter(ref==ref1) %>% 
          ggplot(aes(year,catch_med,lty=variant,fill=variant)) + 
          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.5) + 
          geom_line() + 
          facet_wrap(~pop,scale='free_y',ncol=2) +
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=pop),x=1960,y=Inf, vjust = 2,hjust = 1) +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=20),
                             minor_breaks = seq(1860,2115,by=5))+
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle(ref1))
  dev.off()
}


trial_stat <- NafPerformance('Full_run_res.db')

for(ref1 in unique(trial_stat$pop.res$ref)){
  pdf(file=sprintf('%s/greg.%s.fin.pdf',figs,ref1),width=7,height = 8)
  print( 
    trial_stat$pop.res %>% 
      ungroup() %>% 
      mutate(variant=as.numeric(gsub('V','',variant))) %>% 
      filter(ref==ref1) %>% 
      ggplot(aes(variant,dpl)) + geom_point() + facet_wrap(~pop_id) + 
      theme_minimal() + 
      geom_line(aes(variant,dplfin_60)) +
      geom_line(aes(variant,dplfin_72),lty=2) +
      expand_limits(y = 0) + 
      ggtitle(sprintf('%s final dpl',ref1))
  )
  dev.off()
}

for(ref1 in unique(trial_stat$pop.res$ref)){
  pdf(file=sprintf('%s/greg.%s.min.pdf',figs,ref1),width=7,height = 8)
  print( 
    
    trial_stat$pop.res %>% 
      ungroup() %>% 
      mutate(variant=as.numeric(gsub('V','',variant))) %>% 
      filter(ref==ref1) %>% 
      ggplot(aes(variant,pmin)) + geom_point() + facet_wrap(~pop_id) + 
      theme_minimal() + 
      geom_line(aes(variant,dplmin_60)) +
      geom_line(aes(variant,dplmin_72),lty=2) +
      expand_limits(y = 0) + 
      ggtitle(sprintf('%s minimum dpl',ref1))
  )
  dev.off()
}
