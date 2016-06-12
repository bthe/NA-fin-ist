
sight <- tbl(db,'naf_sight') %>%
  mutate(upper = exp(log(obs)+1.96*cv),
         lower = exp(log(obs)-1.96*cv),
         pro.lower = exp(log(pro.obs)-1.96*cv),
         pro.upper = exp(log(pro.obs)+1.96*cv)) %>%
  select(year,area,obs,lower,upper,pro.obs,pro.lower,pro.upper) %>%
  collect() %>%
  mutate(area=toupper(area),
         area=ifelse(area=='EI','EI/F',area))  %>%
  split(.$area) %>%
  map(~tmp_func(.)) %>% 
  bind_rows() %>%
  group_by(year,area) %>% 
  summarise(obs=sum(obs),lower=sum(lower),upper=sum(upper),
            pro.obs=sum(pro.obs),pro.lower=sum(pro.lower),pro.upper=sum(pro.upper))


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
                              type.key))


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
                              type.key)) %>%
  group_by(ref,rela,area) %>%
  arrange(year) %>%
  mutate(prd1 = cumsum(prd),
         obs1 = cumsum(obs))



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
                              type.key)) %>%
  #  left_join(age.dat) %>%
  #  rename(prd= num) %>%
  group_by(ref,year,sex) %>%
  mutate(obs=obs/sum(obs,na.rm=TRUE),
         prd = prd/sum(prd,na.rm=TRUE),
         diff = obs-prd,
         dir = as.character(sign(diff+1e-10)))

trial_stat <- NafPerformance(db,db_res)

## Example trial NF-B3-1..
res.by.year.pop <- 
  tbl(db_res,'naf_resbyyear') %>% 
  filter(pop_type != 'area') %>%
  collect(n=Inf) %>%
  mutate(msyr = as.numeric(gsub('..-..-([0-9])','\\1',ref)),
         hypo = gsub('..-.([0-9]).+','\\1',ref),
         msyr = as.character(msyr/100),
         trialtype = gsub('..-([A-Z]).+','\\1',ref),
         type = plyr::revalue(trialtype, 
                              type.key)) %>%
  rename(stock=pop_id) %>%
  mutate(stock = ordered(stock,levels=c('W','C1','C2','C3','E','S'))) 

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

