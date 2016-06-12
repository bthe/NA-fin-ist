print( 
  
  age %>% 
    filter(year<1990,grepl('NF-B[0-4]-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~hypo+sex) + theme_minimal() + 
    scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3),
          strip.text = element_text(size=8)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    ggtitle('Catch at age diagnostic for baseline hypotheses (1-4)') + 
    xlab('Year') + ylab('Age (in years)') + scale_size(range=c(0.5,3))
)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-B[5-8]-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~hypo+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3),
          strip.text = element_text(size=8)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for baseline hypotheses (5-6)') + 
    xlab('Year') + ylab('Age (in years)') 
)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[A|D|E]2-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3),
          strip.text = element_text(size=8)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity  (A,E), hypotheses 2') + 
    xlab('Year') + ylab('Age (in years)') 
)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[^(A-E|Y|Q|F)]2-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3),
          strip.text = element_text(size=8)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity trials, hypotheses 2') + 
    xlab('Year') + ylab('Age (in years)') 
)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[A|D|E]3-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3),
          strip.text = element_text(size=8)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity  trials, hypotheses 3') + 
    xlab('Year') + ylab('Age (in years)') 
)
print( 
  
  age %>% 
    filter(year<1990,grepl('NF-[^(A-E|Y|Q|F|S|U)]3-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3),
          strip.text = element_text(size=8)) + 
    scale_x_continuous(breaks = seq(1970,1990,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for sensitivity trials, hypotheses 3') + 
    xlab('Year') + ylab('Age (in years)') 
)

print( 
  
  age %>% 
    filter(grepl('NF-[S|U]3-.',ref)) %>%
    ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
    facet_grid(msyr~type+sex) + theme_minimal() + scale_color_manual(values = c('lightblue','black','black')) + 
    theme(legend.position='none',
          axis.text.x = element_text(angle=90,size = 9,vjust = .3),
          strip.text = element_text(size=8)) + 
    scale_x_continuous(breaks = seq(1970,2010,by=2))+
    scale_size(range=c(0.5,3)) +
    ggtitle('Catch at age diagnostic for varying selectivity, hypotheses 3') + 
    xlab('Year') + ylab('Age (in years)') 
)
