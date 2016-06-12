

for(ms in c(1,4)){
  print(pop %>% 
          filter(grepl(sprintf('NF-B.-%s',ms),ref)) %>%
          mutate(hypo = paste('Hypothesis',hypo)) %>% 
          ggplot(aes(year,number)) + 
          geom_line(aes(col=hypo,lty=hypo)) + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,obs),col='black') + 
          geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=5),
                axis.text.x=element_text(size = 7),
                legend.position = c(0.7,0.1),
                panel.margin = unit(0.3,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2015,by=20),
                             minor_breaks = seq(1860,2015,by=5))+
          expand_limits(y = 0)+
          guides(color = guide_legend(""), lty = guide_legend(""))+
          ggtitle(sprintf('Baseline %s%% (%s)',ms,ifelse(ms==1,'1+','Mature'))))
}

for(ms in c(1,4)){
  for(hyp in 1:3){
    print(pop %>% 
            filter(!(trialtype %in% c('B','Y','Q','J','A','E')), msyr == ms/100,hypo == hyp) %>%
            ggplot(aes(year,number,col=type)) + 
            geom_line(aes(lty=type)) + 
            facet_wrap(~area,scale='free_y',ncol=2) +
            geom_point(aes(year,obs),col='black') + 
            geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
            theme_bw() + ylab('1+ population') + xlab('Year') +
            geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1,col='black') +
            theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=5),
                  axis.text.x=element_text(size = 7),
                  legend.position = c(0.75,0.1),
                  panel.margin = unit(0.3,'cm'),
                  plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                  strip.background = element_blank(),
                  strip.text.x = element_blank())+
            scale_x_continuous(breaks = seq(1860,2015,by=20),
                               minor_breaks = seq(1860,2015,by=5))+
            expand_limits(y = 0)+
            guides(col=guide_legend(title=NULL,ncol=2),lty=guide_legend(title=NULL))+
            ggtitle(sprintf('Other hypothesis %s trials %s%% (%s)',
                            hyp,ms,ifelse(ms==1,'1+','Mature'))))
  }
} 

# plot E trials

for(ms in c(1,4)){
  print(pop %>% 
          filter(trialtype =='E', msyr == ms/100) %>%
          mutate(hypo = paste('Hypothesis',hypo)) %>% 
          ggplot(aes(year,number,col=hypo)) + 
          geom_line(aes(lty=hypo)) + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,ifelse(year<1990 & area %in% c('EI/F','EG','WI'),NA,obs)),col='black') + 
          geom_errorbar(aes(year,ymax=ifelse(year<1990 & area %in% c('EI/F','EG','WI'),NA,upper),
                            ymin=ifelse(year<1990 & area %in% c('EI/F','EG','WI'),NA,lower)),
                        col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=5),
                axis.text.x=element_text(size = 7),
                legend.position = c(0.75,0.1),
                panel.margin = unit(0.3,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2015,by=20),
                             minor_breaks = seq(1860,2015,by=5))+
          expand_limits(y = 0)+
          guides(col=guide_legend(title=NULL,ncol=2),lty=guide_legend(title=NULL))+
          ggtitle(sprintf('Exclude 1987/9 abundance in EG, WI and EI/F %s%% (%s)',
                          ms,ifelse(ms==1,'1+','Mature'))))
  
} 


# plot J trials

for(ms in c(1,4)){
  print(pop %>% 
          filter(trialtype =='J', msyr == ms/100) %>%
          mutate(hypo = paste('Hypothesis',hypo)) %>% 
          ggplot(aes(year,number,col=hypo)) + 
          geom_line(aes(lty=hypo)) + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,obs/0.8),col='black') + 
          geom_errorbar(aes(year,ymax=upper/0.8,
                            ymin=lower/0.8),
                        col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=5),
                axis.text.x=element_text(size = 7),
                legend.position = c(0.75,0.1),
                panel.margin = unit(0.3,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2015,by=20),
                             minor_breaks = seq(1860,2015,by=5))+
          expand_limits(y = 0)+
          guides(col=guide_legend(title=NULL,ncol=2),lty=guide_legend(title=NULL))+
          ggtitle(sprintf('g(0) = 0.8 %s%% (%s)',
                          ms,ifelse(ms==1,'1+','Mature'))))
  
} 


# plot A trials

for(ms in c(1,4)){
  print(pop %>% 
          filter(trialtype =='A', msyr == ms/100) %>%
          mutate(hypo = paste('Hypothesis',hypo)) %>% 
          ggplot(aes(year,number,col=hypo)) + 
          geom_line(aes(lty=hypo)) + 
          facet_wrap(~area,scale='free_y',ncol=2) +
          geom_point(aes(year,pro.obs),col='black') + 
          geom_errorbar(aes(year,ymax=pro.upper,
                            ymin=pro.lower),
                        col='black') + 
          theme_bw() + ylab('1+ population') + xlab('Year') +
          geom_text(aes(label=area),x=1960,y=Inf, vjust = 2,hjust = 1,col='black') +
          theme(axis.text.y=element_text(angle = 90,hjust = 0.5,size=5),
                axis.text.x=element_text(size = 7),
                legend.position = c(0.75,0.1),
                panel.margin = unit(0.3,'cm'),
                plot.margin = unit(c(0.2,0.2,0.2,0.2),'cm'),
                strip.background = element_blank(),
                strip.text.x = element_blank())+
          scale_x_continuous(breaks = seq(1860,2015,by=20),
                             minor_breaks = seq(1860,2015,by=5))+
          expand_limits(y = 0)+
          guides(col=guide_legend(title=NULL,ncol=2),lty=guide_legend(title=NULL))+
          ggtitle(sprintf('Pro-rated abundance %s%% (%s)',
                          ms,ifelse(ms==1,'1+','Mature'))))
}

