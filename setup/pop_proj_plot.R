  ## median population trajectories
  
  dat <- filter(res.by.year.pop,ref==ref1)
  print(dat %>% 
          ggplot(aes(year,abundance_med,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=abundance_upper),col='gray') +
          #geom_line(aes(y=abundance_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          geom_vline(xintercept=2015,lty=2,col='gray') +
          theme_bw() + ylab('Mature female population') + xlab('Year') +
          geom_text(aes(label=stock),x=2000,y=Inf, vjust = 2,hjust = 1,col='black') +
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
          scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
          guides(color = guide_legend(""), lty = guide_legend("")) +
          ggtitle(sprintf('%s population trajectories',ref1)))
  
  
  print(dat %>% 
          filter(variant %in% c('V0','V2','V3','V7')) %>% 
          ggplot(aes(year,abundance_med,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=abundance_upper),col='gray') +
          #geom_line(aes(y=abundance_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          geom_vline(xintercept=2015,lty=2,col='gray') +
          theme_bw() + ylab('Mature female population') + xlab('Year') +
          geom_text(aes(label=stock),x=2000,y=Inf, vjust = 2,hjust = 1,col='black') +
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
          scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#CC79A7")) +
          guides(color = guide_legend(""), lty = guide_legend("")) +
          ggtitle(sprintf('%s population trajectories, variants 0, 2, 3 and 7',ref1)))
  
  
  print(dat %>% 
          filter(year>1999) %>% 
          ggplot(aes(year,dpl_med,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=dpl_upper),col='gray') +
          #geom_line(aes(y=dpl_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          theme_bw() + ylab('Depletion') + xlab('Year') +
          geom_text(aes(label=stock),x=2020,y=Inf, vjust = 2,hjust = 1,col='black') +
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
          scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
          guides(color = guide_legend(""), lty = guide_legend(""))+
          ggtitle(sprintf('%s median population depletion',ref1)))


  print(dat %>% 
          filter(year>1999) %>% 
          ggplot(aes(year,dpl_lower,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=dpl_upper),col='gray') +
          #geom_line(aes(y=dpl_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          theme_bw() + ylab('Depletion') + xlab('Year') +
          geom_text(aes(label=stock),x=2020,y=Inf, vjust = 2,hjust = 1,col='black') +
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
          scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
          guides(color = guide_legend(""), lty = guide_legend(""))+
          ggtitle(sprintf('%s 5%%-ile population depletion',ref1)))
  
  print(dat %>% 
          filter(year>1999) %>% 
          ggplot(aes(year,dpl_upper,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=dpl_upper),col='gray') +
          #geom_line(aes(y=dpl_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          theme_bw() + ylab('Depletion') + xlab('Year') +
          geom_text(aes(label=stock),x=2020,y=Inf, vjust = 2,hjust = 1,col='black') +
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
          scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
          guides(color = guide_legend(""), lty = guide_legend(""))+
          ggtitle(sprintf('%s 95%%-ile population depletion',ref1)))
  
  
  print(dat %>% 
          filter(year>1999,variant %in% c('V0','V2','V3','V7')) %>% 
          ggplot(aes(year,dpl_med,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=dpl_upper),col='gray') +
          #geom_line(aes(y=dpl_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          theme_bw() + ylab('Depletion') + xlab('Year') +
          geom_text(aes(label=stock),x=2020,y=Inf, vjust = 2,hjust = 1,col='black') +
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
          scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#CC79A7")) +
          guides(color = guide_legend(""), lty = guide_legend(""))+
          ggtitle(sprintf('%s median population depletion, variants 0, 2, 3 and 7',ref1)))
  
  
  print(dat %>% 
          filter(year>1999,variant %in% c('V0','V2','V3','V7')) %>% 
          ggplot(aes(year,dpl_lower,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=dpl_upper),col='gray') +
          #geom_line(aes(y=dpl_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          theme_bw() + ylab('Depletion') + xlab('Year') +
          geom_text(aes(label=stock),x=2020,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=50),
                             minor_breaks = seq(1860,2115,by=25))+
          scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#CC79A7")) +
          expand_limits(y = 0)+
          guides(color = guide_legend(""), lty = guide_legend(""))+
          ggtitle(sprintf('%s 5%%-ile population depletion, variants 0, 2, 3 and 7',ref1)))
  
  print(dat %>% 
          filter(year>1999,variant %in% c('V0','V2','V3','V7')) %>% 
          ggplot(aes(year,dpl_upper,lty=variant,col=variant)) + 
          #geom_ribbon(aes(year,ymax=abundance_upper,ymin=abundance_lower),alpha=0.2) + 
          #geom_line(aes(y=dpl_upper),col='gray') +
          #geom_line(aes(y=dpl_lower),col='gray') +
          geom_line() + 
          facet_wrap(~stock,scale='free_y',ncol=2) +
          theme_bw() + ylab('Depletion') + xlab('Year') +
          geom_text(aes(label=stock),x=2020,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=50),
                             minor_breaks = seq(1860,2115,by=25))+
          scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#CC79A7")) +
          expand_limits(y = 0)+
          guides(color = guide_legend(""), lty = guide_legend(""))+
          ggtitle(sprintf('%s 95%%-ile population depletion, variants 0, 2, 3 and 7',ref1)))
  
  
    

  print(res.by.year %>% 
        filter(area %in% c('WI','EI/F'),ref==ref1) %>% 
        ggplot(aes(year,catch_med,lty=variant,col=variant)) + 
        #          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.2) + 
        #geom_line(aes(y=catch_upper),col='gray') +
        #geom_line(aes(y=catch_lower),col='gray') +
        geom_line() + 
        facet_wrap(~area,scale='free_y',ncol=1) +
        #          geom_point(aes(year,obs),col='black') + 
        #          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
        theme_bw() + ylab('Catches') + xlab('Year') +
        geom_text(aes(label=area),x=2100,y=Inf, vjust = 2,hjust = 1,col='black') +
        theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
              axis.text.x=element_text(size = 7),
              #legend.position = 'none',#c(0.7,0.2),
              panel.margin = unit(0.2,'cm'),
              plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
              strip.background = element_blank(),
              strip.text.x = element_blank())+
        scale_x_continuous(breaks = seq(1860,2115,by=50),
                           minor_breaks = seq(1860,2115,by=25))+
          scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
        expand_limits(y = 0)+
        #scale_fill_manual(values = c( "darkred", "darkgreen"))+
        #scale_fill_hue(l=40)+
        ggtitle('NF-B3-1 median catch trajectories'))

  
  print(res.by.year %>% 
          filter(area %in% c('WI','EI/F'),ref==ref1,year>1999) %>% 
          ggplot(aes(year,catch_med,lty=variant,col=variant)) + 
          #          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.2) + 
          #geom_line(aes(y=catch_upper),col='gray') +
          #geom_line(aes(y=catch_lower),col='gray') +
          geom_line() + 
          facet_wrap(~area,scale='free_y',ncol=1) +
          #          geom_point(aes(year,obs),col='black') + 
          #          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('Catches') + xlab('Year') +
          geom_text(aes(label=area),x=2100,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=50),
                             minor_breaks = seq(1860,2115,by=25))+
          scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle('NF-B3-1 median catch trajectories'))
  
  
  
  print(res.by.year %>% 
          filter(area %in% c('WI','EI/F'),ref==ref1,year>1999,
                 variant %in% c('V0','V2','V3','V7')) %>% 
          ggplot(aes(year,catch_med,lty=variant,col=variant)) + 
          #          geom_ribbon(aes(year,ymax=catch_upper,ymin=catch_lower),alpha=0.2) + 
          #geom_line(aes(y=catch_upper),col='gray') +
          #geom_line(aes(y=catch_lower),col='gray') +
          geom_line() + 
          facet_wrap(~area,scale='free_y',ncol=1) +
          #          geom_point(aes(year,obs),col='black') + 
          #          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('Catches') + xlab('Year') +
          geom_text(aes(label=area),x=2100,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=8),
                axis.text.x=element_text(size = 7),
                #legend.position = 'none',#c(0.7,0.2),
                panel.margin = unit(0.2,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2115,by=50),
                             minor_breaks = seq(1860,2115,by=25))+
          scale_color_manual(values = c("#000000", "#56B4E9", "#009E73", "#CC79A7")) +
          expand_limits(y = 0)+
          #scale_fill_manual(values = c( "darkred", "darkgreen"))+
          #scale_fill_hue(l=40)+
          ggtitle('NF-B3-1 median catch trajectories, variants 0, 2, 3 and 7'))
  
if(FALSE){
  
  
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
  
  res0 <- tmp_func('Full_run2/%s-V0.restest')
  res1 <- tmp_func('Full_run2/%s-V1.restest')
  res2 <- tmp_func('Full_run2/%s-V2.restest')
  res3 <- tmp_func('Full_run2/%s-V3.restest')
  res4 <- tmp_func('Full_run2/%s-V4.restest')
  res5 <- tmp_func('Full_run2/%s-V5.restest')
  res6 <- tmp_func('Full_run2/%s-V6.restest')
  res7 <- tmp_func('Full_run2/%s-V7.restest')
  res <- bind_rows(res0,res1,res3,res4,res5,res6,res7)
  
for(var1 in unique(res$variant)){
  print(res %>% 
          filter(variant == var1, trial %in% 1:10,
                 pop_type!='area') %>%
          rename(stock=pop_id) %>% 
          ggplot(aes(year,abundance,group=trial)) +   
          geom_line() + 
          geom_vline(xintercept = 2015,lty=2,col='gray') +
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
          ggtitle(sprintf('%s 10 populations trajectories - variant %s',ref1,var1) ))
}





}

  