print(
  wi_tag %>% 
    filter(grepl('NF-B.-.',ref)) %>%
    ungroup() %>% 
    mutate(area = toupper(area),
           rela = toupper(rela)) %>% 
    filter(area == 'WI',rela =='WI') %>%
    ggplot(aes(year,ifelse(obs>0,obs1,NA),lty=msyr)) + geom_point() + 
    geom_line(aes(year,prd1)) + 
    facet_wrap(~hypo)+
    theme_bw() + xlim(c(1965,2010)) + ylab('Num. tags recaptured in WI') +
   # geom_text(aes(label=hypo),x=1972,y=Inf, vjust = 1.5,hjust = 1,col='black') +
    guides(lty=guide_legend(title='MSYR'))+
    theme(legend.position = c(0.8,0.2), 
          axis.text.x = element_text(vjust=0.7,angle=45),
          strip.background = element_blank())+
          #strip.text.x = element_blank()) +
    xlab('Year') + ggtitle('Baseline trials tag recaptures within WI')
)


print(
  wi_tag %>% 
    filter(grepl('NF-[^(B|Y|Q|G)].-.',ref)) %>%
    ungroup() %>% 
    mutate(area = toupper(area),
           rela = toupper(rela),
           type = ifelse(trialtype=='S','Selectivity estimated pre & post 2000',type)) %>% 
    filter(area == 'WI',rela =='WI') %>%
    #filter(trialtype=='B',hypo == 8) %>%
    ggplot(aes(year,ifelse(obs>0,obs1,NA),lty=msyr)) + geom_point() + 
    geom_line(aes(year,prd1)) + 
    facet_wrap(~type+hypo)+
    theme_bw() + xlim(c(1965,2010)) + ylab('Num. tags recaptured in WI') +
    guides(lty=guide_legend(title='MSYR'))+
    theme(legend.position = c(0.6,0.1), 
          axis.text.x = element_text(vjust=0.7,angle=45),
          strip.background = element_blank(),
          strip.text = element_text(size=7)) +
    xlab('Year') + ggtitle('Other trials tag recaptures')
)


