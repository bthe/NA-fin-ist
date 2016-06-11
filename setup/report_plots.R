## setup
library(ggplot2)
library(openxlsx)
source('setup/trials.R')
source('setup/read_variants.R')
figs <- 'figs_report'
db <- src_sqlite('Full_run.db')
db_res <- src_sqlite('Full_run_res.db')
type.key <- 
  c('A'='Pro rate abundance','B'='Baseline','C'='CPUE',
    'F'='C2 to EG in 1985 - 2025 (opt. b)',
    'G'='C2 to EG in 1985 (opt. a)',
    'H'='High historical catch series',
    'T'='Tag loss 20% yr 1, 10%/yr thereafter',
    'J'='g(0) = 0.8',
    'U'='Dome shaped selectivity',
    'Y'='8 year future survey interval',
    'D'='Upper bound on dispersal')

dir.create(figs)

tmp_func <- function(file){
  res <- readLines(file)[-c(2:28)]
  header <- scan(text=gsub(': ','.',res[2]),what='character',quiet = TRUE)
  ref <- gsub('\\s*(NF-..-.).+','\\1',res[1])
  hypo <- as.numeric(gsub('NF-.([0-9])-.','\\1',ref))
  variant <- gsub('\\s*NF-..-.\\s(V[0-9]).+','\\1',res[1])
  loc <- grep('Trial',res)
  tmp <- gsub('(^\\s*$|^\\s*[a-z].+|^-1|^\\s*0)','# \\1',tolower(res)) %>% 
    read.table(text=.)
  names(tmp) <- tolower(header)
  
  stock.names <- c('W','C1','C2','C3','E','S')
  area.names <- c('EC','WG','EG','WI','EI/F','N','SP')
  if(hypo %in% 7:8){
    stock.names <- c('W','C1','C2','E','S')
    area.names <- c('EC','WG','EG+WI','EI/F','N','SP')
  }
  if(hypo == 6){
    stock.names <- c('W','C1','C2','C3','S')
  }
  
  #   C    1  Total catch over management period by subarea (CTOT & C2TOT): 
  #   C    2  Initial stock size (P0) (for each stock)
  #   C    3  Final population size (PFIN) (for each stock)
  #   C    4  Minimum population sizes reached in each trial (PMIN)
  #   C    5  Average catch over the first 10yrs by subarea (CF10)
  #   C    6  Average catch over the last  10yrs by subarea (CL10)
  
  tmp %>%
    rename(year=yr) %>% 
    mutate(ref = ref,
           variant = variant,
           trial = cut(1:length(year),c(0,which(diff(year)<0),1e9),labels = FALSE)-1) %>% 
    gather(key,number,-c(year,ref,trial,variant)) %>% 
    separate(key,c('pop_type','pop_id')) %>% 
    mutate(pop_type = plyr::revalue(pop_type,c('fem'='Mature females','pk'='1+','ck'='Area catch','cj'='Stock catch')),
           pop_id = ifelse(pop_type %in% c('Mature females','Stock catch'),stock.names[as.numeric(pop_id)],
                           area.names[as.numeric(pop_id)]),
           pop_type = plyr::revalue(pop_type,c('Mature females'='pop abundance','1+'='area abundance',
                                               'Area catch'='area catch','Stock catch'='pop catch'))) %>%
    separate(pop_type,c('pop_type','tmp')) %>% 
    spread(tmp,number) %>%
    filter(trial!=0)  
}

res0 <- tmp_func('Full_run/NF-B1-1-V0.restest')
res1 <- tmp_func('Full_run/NF-B1-1-V1.restest')
res2 <- tmp_func('Full_run/NF-B1-1-V2.restest')
res3 <- tmp_func('Full_run/NF-B1-1-V3.restest')
res4 <- tmp_func('Full_run/NF-B1-1-V4.restest')
res5 <- tmp_func('Full_run/NF-B1-1-V5.restest')
res6 <- tmp_func('Full_run/NF-B1-1-V6.restest')
res7 <- tmp_func('Full_run/NF-B1-1-V7.restest')
res <- bind_rows(res0,res1,res3,res4,res5,res6,res7)

## performance plots for MSYR_1+=1%

trial_stat <- NafPerformance(db,db_res)

pdf(file=sprintf('%s/performance.pdf',figs),width=7,height = 8)

for(ref1 in unique(trial_stat$pop.res$ref)){
  if(substr(ref1,7,7)=='4')
    next
  print( 
    trial_stat$pop.res %>% 
      ungroup() %>% 
      mutate(variant=as.numeric(gsub('V','',variant))) %>% 
      filter(ref==ref1) %>% 
      ggplot(aes(variant,dpl)) + 
      geom_rect(aes(ymax=dplfin_60),xmin=-Inf,xmax=Inf,ymin=-Inf,
                fill = 'gray70') +
      geom_hline(aes(yintercept=dplfin_72)) +
      geom_point() + facet_wrap(~pop_id) + 
      theme_bw() + 
      expand_limits(y = 0) + 
      ggtitle(sprintf('%s final depletion',ref1)) + 
      ylab('Final depleption') + xlab('Variant') + 
      scale_x_continuous(breaks=0:7))
   print( 
    trial_stat$pop.res %>% 
      ungroup() %>% 
      mutate(variant=as.numeric(gsub('V','',variant))) %>% 
      filter(ref==ref1) %>% 
      ggplot(aes(variant,pmin)) + 
      geom_rect(aes(ymax=dplmin_60),xmin=-Inf,xmax=Inf,ymin=-Inf,
                fill = 'gray70') +
      geom_hline(aes(yintercept=dplmin_72)) +
      geom_point() + facet_wrap(~pop_id) + 
      theme_bw() + 
      expand_limits(y = 0) + 
      ggtitle(sprintf('%s minimum depletion ratio',ref1)) + 
      ylab('Minimum depleption ratio') + xlab('Variant') + 
      scale_x_continuous(breaks=0:7)
  )
}
dev.off()

## Example trial NF-B3-1..
res.by.year.pop <- 
  tbl(db_res,'naf_resbyyear') %>% 
  filter(ref=='NF-B3-1') %>% 
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

pdf(file=sprintf('%s/NF-B3-1.pdf',figs,ref1),width=7,height = 8)
print(res.by.year.pop %>% 
        ggplot(aes(year,abundance_med,lty=variant)) + 
        #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
        geom_line(aes(y=abundance_upper),col='gray') +
        geom_line(aes(y=abundance_lower),col='gray') +
        geom_line() + 
        facet_wrap(~stock,scale='free_y',ncol=2) +
        theme_bw() + ylab('Mature female population') + xlab('Year') +
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
        
        ggtitle('NF-B3-1 population trajectories'))


print(res.by.year.pop %>% 
        filter(year>1999) %>% 
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
        scale_x_continuous(breaks = seq(1860,2115,by=50),
                           minor_breaks = seq(1860,2115,by=25))+
        expand_limits(y = 0)+
        ggtitle('NF-B3-1 population depletion'))

print(res.by.year %>% 
        filter(year>1930,variant != 'V0',
               area %in% c('WI','EI/F')) %>% 
        ggplot(aes(year,catch_med,lty=variant)) + 
        #          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.2) + 
        #geom_line(aes(y=catch_upper),col='gray') +
        #geom_line(aes(y=catch_lower),col='gray') +
        geom_line() + 
        facet_wrap(~area,scale='free_y',ncol=1) +
        #          geom_point(aes(year,obs),col='black') + 
        #          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
        theme_bw() + ylab('Catches') + xlab('Year') +
        geom_text(aes(label=area),x=2100,y=Inf, vjust = 2,hjust = 1) +
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
        ggtitle('NF-B3-1 median catch trajectories'))


print(res.by.year %>% 
        filter(year>1999,variant != 'V0',
               area %in% c('WI','EI/F')) %>% 
        ggplot(aes(year,catch_med,lty=variant)) + 
        #          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.2) + 
        #geom_line(aes(y=catch_upper),col='gray') +
        #geom_line(aes(y=catch_lower),col='gray') +
        geom_line() + 
        facet_wrap(~area,scale='free_y',ncol=1) +
        #          geom_point(aes(year,obs),col='black') + 
        #          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
        theme_bw() + ylab('Catches') + xlab('Year') +
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
        ggtitle('NF-B3-1 median catch trajectories'))

for(var1 in unique(res$variant)){
  print(res %>% 
          filter(variant == var1, trial %in% 1:10,
                 pop_type!='area') %>%
          rename(stock=pop_id) %>% 
          ggplot(aes(year,abundance,group=trial)) + geom_line() +  
          geom_line() + 
          facet_wrap(~stock,scale='free_y') +
          #          geom_point(aes(year,obs),col='black') + 
          #          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('Mature female population') + xlab('Year') +
          geom_text(aes(label=stock),x=2000,y=Inf, vjust = 2,hjust = 1) +
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
          ggtitle(sprintf('NF-B3-1 10 populations trajectories - variant %s',var1) ))
}


dev.off()
    
## tables

write.xlsx(trial_stat,file=sprintf('%s/trialperformance.xlsx',figs))




