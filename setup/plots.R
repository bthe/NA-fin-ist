library(ggplot2)
tmp_func <- function(x){
  if(x$area[1] %in% c('WI','EG')){
    tmp <- x
    tmp$area <- 'EG+WI'
    bind_rows(x,tmp)
  } else{
    x
  }
}

db <- src_sqlite('druna.db')
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
  


pop <- tbl(db,'naf_pop') 


pop %>% 
  filter(pop_type == '1+') %>%
  collect() %>%
  group_by(year,ref,hypo,msyr,pop_id) %>%
  summarise(med = median(number),
            cond.975 = quantile(number,0.975),
            cond.0.25 = quantile(number,0.025),
            number = sum(ifelse(trial==0,number,0))) %>%
  ungroup() %>%
  mutate(
    #agewt = gsub('T.-([0-9]).-[0-9]','\\1',ref),
    #tagwt = gsub('T([0-9])-..-[0-9]','\\1',ref),
    #         msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
    #         hypo = gsub('..-.([0-9]).+','\\1',ref),
    msyr = as.character(msyr/100),
    trialtype = gsub('..-([A-Z]).+','\\1',ref)) %>%
  rename(area=pop_id) %>%
  left_join(sight) %>%
  mutate(area = ordered(area,levels=c('EC','WG','EG','WI','EG+WI','EI/F','N','SP')))  -> tmp

tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(hypo == 3,trialtype%in%c('B','C')) %>% #,agewt==5,tagwt==7) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=msyr,fill=msyr,col=trialtype)) + 
  geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line(col='black') + 
  geom_line(aes(year,med),col='red') + 
  facet_wrap(~area) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')

full.ec <- 
tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(area == 'EC',trialtype %in% c('B')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  #ggplot(aes(year,number,lty=tagwt,col=agewt)) + 
  ggplot(aes(year,number,lty=msyr)) + 
  #geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~hypo) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')

red.ec <-
tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(area == 'EC',agewt%in%5:7,tagwt%in%7:8) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=tagwt,col=agewt)) + 
  #geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~hypo+msyr) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')

full.eif <- 
tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(area == 'EI/F') %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=tagwt,col=agewt)) + 
  #geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~hypo+msyr) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')

red.eif <- 
tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(area == 'EI/F',tagwt%in%5:7,agewt %in% 5:6) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=tagwt,col=agewt)) + 
  #geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~hypo+msyr) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')


full.eg <- 
  tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(area %in% c('EG','EG+WI')) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=tagwt,col=agewt)) + 
  #geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~hypo+msyr) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')

red.eg <- 
  tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(area %in% c('EG','EG+WI'),tagwt%in%5:7,agewt %in% 5:6) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=tagwt,col=agewt)) + 
  #geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~hypo+msyr) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')




pop %>% 
  filter(trial == 0,pop_type != '1+') %>%
  collect() %>%
  mutate(tagwt = gsub('T.-([0-9]).-[0-9]','\\1',ref),
         agewt = gsub('T([0-9])-..-[0-9]','\\1',ref)) %>%
  rename(stock=pop_id) %>%
  left_join(sight) -> tmp2
tmp2 %>% ## hypos 1,(-)2,3(-),(-)4,5(-),6(-),7,(-)8
  filter(hypo %in% c(3,7),msyr==4,tagwt==1,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=as.factor(hypo))) + geom_line() + 
  facet_wrap(~stock,scale='free_y')  + 
  theme_bw() + ylab('Mature females') + xlab('Year')


params <- tbl(db,'naf_all')

params %>% filter(trial == 'T9-16-4') %>% collect() %>% View()

age.dat <- read.table('data/catchbyage.dat',na.strings = '0') %>%
  gather(var,val,-c(V1:V2)) %>%
  na.omit() %>% 
  select(year=V1,age=V2,sex=var,obs=val) %>%
  mutate(sex = ifelse(sex=='V9','Males','Females')) %>% 
  group_by(year,sex) %>%
  mutate(obs=obs/sum(obs,na.rm=TRUE))

age <- 
  tbl(db,'naf_age') %>%
  filter(trial==0) %>% 
  collect() %>% 
#  tmp %>%
  gather(var,num, -c(year:age,trial,ref)) %>%
  separate(var,c('type','sex')) %>% 
#  filter(type!='obs') %>%
  spread(type,num) %>%
  mutate(sex = ifelse(sex=='m','Males','Females'),
         msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref)) %>%
#  left_join(age.dat) %>%
#  rename(prd= num) %>%
  group_by(ref,year,sex) %>%
  mutate(prd = prd/sum(prd,na.rm=TRUE),
         obs=obs/sum(obs,na.rm=TRUE),
         diff = obs-prd,
         dir = as.character(sign(diff+1e-10))) 

age %>% 
  filter(year<1990,trialtype=='B',hypo == 8) %>%
  ggplot(aes(year,age,size=abs(diff),col=dir)) + geom_point()+
  facet_grid(msyr~sex) + theme_bw() + scale_color_manual(values = c('lightblue','black')) + 
  theme(legend.position='none') 


tbl(db,'naf_tag') %>%
  filter(area %in% c('WI','EG+WI')) %>%
  collect() %>% 
  mutate(msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         rela = ordered(toupper(rela),levels=c('EC','WG','EG','WI','EG+WI','EI/F','N','SP'))) %>%
  group_by(ref,rela,area) %>%
  arrange(year) %>%
  mutate(prd1 = cumsum(prd),
         obs1 = cumsum(obs)) %>%
  filter(trialtype=='R',hypo == 4) %>%
    ggplot(aes(year,ifelse(obs>0,obs1,NA))) + geom_point() + 
  geom_line(aes(year,prd1,lty=msyr)) + 
  facet_wrap(~rela,scale='free_y',ncol=1)+
  theme_bw() + xlim(c(1965,2010)) + ylab('Num. tags recaptured in WI') +
  xlab('Year')
  
br4.plot <-   
  tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(hypo == 4,trialtype%in%c('B','R')) %>% #,agewt==5,tagwt==7) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=msyr,col=trialtype)) + 
  #  geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
#  geom_line(aes(year,med),col='red') + 
  facet_wrap(~area) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')


br3.plot <-   
  tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(hypo == 3,trialtype%in%c('B','R')) %>% #,agewt==5,tagwt==7) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=msyr,col=trialtype)) + 
  #  geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~area) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')


bg3.plot <-   
  tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(hypo == 3,trialtype%in%c('B','G')) %>% #,agewt==5,tagwt==7) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=msyr,col=trialtype)) + 
  #  geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~area) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')


bt3.plot <-   
  tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(hypo == 3,trialtype%in%c('B','T')) %>% #,agewt==5,tagwt==7) %>% #,trialtype %in% c('B','R')) %>% #,tagwt==9,agewt %in% c(9)) %>%
  ggplot(aes(year,number,lty=msyr,col=trialtype)) + 
  #  geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~area) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')

hypo3var <- 
tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(hypo == 3,trialtype%in%c('B','T','G','R','C'),area %in% c('EC','EG','WI','EI/F')) %>% 
  ggplot(aes(year,number,lty=msyr,col=trialtype)) + 
  #  geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~area) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')

hypo4var <- 
tmp %>% ## hypos 1,(-)2,3(-),(-)4,5(-),(-)6,7,(-)8
  filter(hypo == 4,trialtype%in%c('B','T','G','R','C'),area %in% c('EC','EG','WI','EI/F')) %>% 
  ggplot(aes(year,number,lty=msyr,col=trialtype)) + 
  #  geom_ribbon(aes(year,ymax=cond.975,ymin=cond.0.25),alpha=0.3) + 
  geom_line() + 
  #  geom_line(aes(year,med),col='red') + 
  facet_wrap(~area) +
  geom_point(aes(year,obs),col='black') + geom_errorbar(aes(year,ymax=upper,ymin=lower),col='black') + 
  theme_bw() + ylab('1+ population') + xlab('Year')



