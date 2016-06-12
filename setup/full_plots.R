library(ggplot2)
source('setup/trials.R')
source('setup/read_variants.R')
figs <- 'figs_Full_run2'
db <- src_sqlite('Full_run2.db')
db_res <- src_sqlite('Full_run_res2.db')
type.key <- 
  c('A'='Pro rate abundance','B'='Baseline','C'='CPUE',
    'F'='C2 to EG in 1985 - 2025 (opt. b)',
    'G'='C2 to EG in 1985 (opt. a)',
    'H'='High historical catch series',
    'T'='Tag loss 20% yr 1, 10%/yr thereafter',
    'J'='g(0) = 0.8',
    'U'='Dome shaped selectivity',
    'Y'='8 year future survey interval',
    'D'='Upper bound on dispersal',
    'E'='Exclude 1987/9 abundance in WI, EG & EI/F',
    'Q'='Future WI & EI/F surveys exc. strata S 60N',
    'S'='Selectivity estimated for pre and post 2000 & use all age data')

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
  a95 <- rep(0,length(d))
  for(i in length(d)){ 
    tmp <- 1/(1 + exp(-(a-a50[i])/d[i]))
    a95[i] <- max(a[which(tmp<0.95)])
  }
  return(a95)
}

tmp2 <- 
  tbl(db,'naf_all') %>%
  collect(n=Inf) %>% 
  filter(!grepl('NF-S',trial),grepl('disp',constant)|constant %in% c('rmsig','rfsig','rm50','rf50')) %>%
  spread(constant,value) %>% 
 # head() %>% 
#  split(c(.$trial,.$ref)) %>% 
#  map(~
  mutate(sel.male.50 = rm50,
         sel.male.95 = ogive(rm50,rmsig),
         sel.female.50 = rf50,
         sel.female.95 = ogive(rf50,rfsig)) %>% 
  arrange(trial,ref)
  select(-c(rf50:rmsig)) %>% 
  gather(constant,value,-c(ref,trial)) %>%
  filter(!(constant %in% c('sel.female.50','sel.female.95'))| 
           (constant %in% c('sel.female.50','sel.female.95') & value < 10))
#  mutate(constant = plyr::revalue(constant,
#                                  c('disp1'='C1C2','disp2'='C2C3',
#                                  'sel.male.50'='Male a50','sel.female.50'='Female a50',
#                                  'sel.male.95'='Male a95','sel.female.95'='Female a95'))) 
pdf(file=sprintf('%s/parameter_distribution.pdf',figs),width=11,height = 11)

  for(cons in unique(tmp2$constant)){
    print(ggplot(filter(tmp2,constant==cons),aes(value)) + geom_histogram() + 
      facet_wrap(~trial,scale='free') +
      geom_vline(data=filter(tmp2,ref=='0',constant==cons),aes(xintercept=value))+
      theme_bw() + ggtitle(cons) + ylab('Num trials') +
      xlab('Value'))
  }
dev.off()

source('setup/load_data.R')

pdf(file=sprintf('%s/pop1plus.pdf',figs),width=11,height = 8)
source('setup/pop1plus.R')
dev.off()

pdf(file=sprintf('%s/age-bubble.pdf',figs),width=11,height = 8)
source('setup/age-bubble.R')
dev.off()

pdf(file=sprintf('%s/tagging.pdf',figs),width=11,height = 8)
source('setup/tagging.R')
dev.off()

pdf(file=sprintf('%s/performance-plots.pdf',figs),width=11,height = 8)
source('setup/performance-plots.R')
dev.off()

for(ref1 in unique(res.by.year.pop$ref)){
  pdf(file=sprintf('%s/res-%s.pdf',figs,ref1),width=11,height = 8)
  source('setup/performance-plots.R')
  dev.off()
}



for(ref1 in unique(gsub('(..-..).+','\\1',pop$ref))){
  pdf(file=sprintf('%s/%s.pop.pdf',figs,ref1),width=7,height = 8)
  print(pop %>% 
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
# 
# 
# pdf(file=sprintf('%s/cpue.pdf',figs),width=11,height = 8)
# print(
#   tbl(db,'naf_cpe') %>% 
#     collect() %>% 
#     ggplot(aes(year,obs,lty=ref)) + geom_point() + geom_line(aes(y=prd)) +
#     facet_wrap(~series,scale='free') + 
#     theme_minimal() + theme(legend.position ='none') +
#     ylab('CPUE') + xlab('Year') +
#     ggtitle('Fit to CPUE series')
# )
# dev.off()
# 
# 
# pdf(file=sprintf('%s/cpue.xy.pdf',figs),width=11,height = 8)
# print(
#   tbl(db,'naf_cpe') %>% 
#     collect() %>% 
#     ggplot(aes(prd,obs,col=ref)) + geom_text(aes(label=year)) + geom_abline(slope=1,intercept=0) +
#     facet_wrap(~series,scale='free') + 
#     theme_minimal() + #theme(legend.position ='none') +
#     ylab('CPUE') + xlab('Year') +
#     ggtitle('Fit to CPUE series')
# )
# dev.off()



res.by.year <- 
  tbl(db_res,'naf_resbyyear') %>% 
  filter(pop_type == 'area') %>%
  collect(n=Inf) %>%
  mutate(msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         type = plyr::revalue(trialtype, 
                             type.key)) %>%
  rename(area=pop_id) %>%
  left_join(sight) %>%
  mutate(area = ordered(area,levels=c('EC','WG','EG','WI','EG+WI','EI/F','N','SP')))


res.by.year.pop <- 
  tbl(db_res,'naf_resbyyear') %>% 
  collect(n=Inf) %>% 
  filter(pop_type != 'area') %>%
  collect(n=Inf) %>%
  mutate(msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         type = plyr::revalue(trialtype, 
                              type.key)) %>%
  rename(stock=pop_id) %>%
  mutate(area = ordered(stock,levels=c('W','C1','C2','C3','E','S'))) 
  


for(ref1 in unique(res.by.year$ref)){
  pdf(file=sprintf('%s/%s.res.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year %>% 
          filter(ref==ref1) %>% 
          ggplot(aes(year,abundance_med,lty=variant)) +
          geom_line(aes(y=abundance_upper)) +
          geom_line(aes(y=abundance_lower)) +
          geom_line() +
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0) + 
          geom_line() + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,obs),col='black') + 
          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=2080,y=-Inf, vjust = -1,hjust = 1) +
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



for(ref1 in unique(res.by.year.pop$ref)){
  pdf(file=sprintf('%s/%s.res.pop.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year.pop %>% 
          filter(ref==ref1,year>1995) %>% 
          ggplot(aes(year,dpl_med,lty=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          geom_line(aes(y=dpl_upper),col='gray') +
          geom_line(aes(y=dpl_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          theme_bw() + ylab('Depletion') + xlab('Year') +
          geom_text(aes(label=area),x=2000,y=Inf, vjust = 2,hjust = 1) +
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

          ggtitle(ref1))
  dev.off()
}


for(ref1 in unique(res.by.year$ref)){
  for(var1 in unique(res.by.year$variant)){
    pdf(file=sprintf('%s/%s.%s.res.catbyar.pdf',figs,ref1,var1),width=7,height = 8)
    print(res.by.year %>% 
            filter(ref==ref1,year>1930,variant == var1,
                   area %in% c('WI','EI/F')) %>% 
            ggplot(aes(year,catch_med)) + 
            #          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.2) + 
            geom_line(aes(y=catch_upper),col='gray') +
            geom_line(aes(y=catch_lower),col='gray') +
            geom_line() + 
            facet_wrap(~area,scale='free_y',ncol=1) +
            #          geom_point(aes(year,obs),col='black') + 
            #          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
            theme_bw() + ylab('1+ population') + xlab('Year') +
            geom_text(aes(label=area),x=2100,y=Inf, vjust = 2,hjust = 1) +
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
            ggtitle(sprintf('%s - variant %s',ref1,var1)))
    dev.off()
  }
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
                              type.key)) %>%
  rename(pop=pop_id) %>%
  mutate(pop = ordered(pop,levels=c('W','C1','C2','C3','E','S')))



for(ref1 in unique(res.by.year$ref)){
  pdf(file=sprintf('%s/%s.res.fem.pdf',figs,ref1),width=7,height = 8)
  print(res.by.year.fem %>% 
          filter(ref==ref1) %>% 
          ggplot(aes(year,abundance_med,lty=variant)) +
          geom_line(aes(y=abundance_upper),col='gray') +
          geom_line(aes(y=abundance_lower),col='gray') +
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.5) + 
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
          scale_x_continuous(breaks = seq(1860,2115,by=50),
                             minor_breaks = seq(1860,2115,by=25))+
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


