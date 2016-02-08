db <- src_sqlite('trials.db')
sight <- tbl(db,'naf_sight') %>%
  mutate(upper = obs*(1+1.96*cv),
         lower = obs*(1-1.96*cv),
         ref='NF-B1-1') %>%
  select(year,area,number=obs,lower,upper,ref) %>%
  collect()

pop <- tbl(db,'naf_pop') 

area_map <- data.frame(pop_id=1:7,area = c('ec','wg','eg','wi','ei','n','sp'))

pop %>% 
  filter(trial == 0,pop_type == '1+',ref %in% c('NF-B1-1','NF-B1-4')) %>%
  left_join(area_map,copy=TRUE) %>%
  collect() %>%
  ggplot(aes(year,number,lty=ref)) + geom_line() + facet_wrap(~area,scale='free_y') +
  geom_point(data=sight) + geom_errorbar(data=sight,aes(year,ymax=upper,ymin=lower)) + 
  theme_bw()