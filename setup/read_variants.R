NafReadVariants <- function(dir='trials',search.string='NF-..-.-V[0-9].restest',db_name='trial.db'){
  
  db <- src_sqlite(db_name, create = T)
  
  print('read in management variant results')
  #db_drop_table(trials_db$con,'naf_pop')
  db_create_table(db$con,'naf_res',
                  types = c(ref = 'TEXT',variant = 'TEXT',
                            trial = 'INTEGER',pop_type = 'TEXT',
                            pop_id = 'TEXT',total_catch = 'REAL',
                            avg_catch = 'REAL',p0 = 'REAL',pfin = 'REAL',
                            pmin = 'REAL',cf10 = 'REAL',cl10 = 'REAL',final_dpl = 'REAL'))
  res <- 
    list.files(path = dir,pattern =search.string)%>%
    mclapply(., function(x) safely(NafReadVariants.restest)(sprintf('%s/%s',dir,x)),
             mc.cores = detectCores(logical=TRUE)) %>% 
    map('result') %>% 
    bind_rows() %>% 
    as.data.frame() %>% 
    db_insert_into(db$con,table='naf_res',values = .)
  
  
}

NafReadVariants.restest <- function(file='NAF.restest'){
  res <- readLines(file)[-c(2:28)]
  header <- scan(text=gsub(': ','.',res[2]),what='character')
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
    filter(trial!=0) %>% 
    group_by(ref,variant,trial,pop_type,pop_id) %>% 
    summarise(total_catch=sum(catch[year>2014]),avg_catch=mean(catch[year>2014]), 
              p0 = sum(abundance[year==1864]),pfin = sum(abundance[year==2115]),
              pmin = min(abundance[year>2014]), cf10=mean(catch[year %in% 2015:2024]),
              cl10 = mean(catch[year %in% 2105:2114]),
              final_dpl = pfin/p0) 
}


NafPerformance <- function(db_name='trials.db',dir='trials'){
  db <- src_sqlite(db_name)
  res <- tbl(db,'naf_res')
  thr <- tbl(db,'man_thresh')
    
  pop.res <- 
    res %>% 
    filter(pop_type =='pop') %>% 
    collect() %>%
    group_by(ref,variant,pop_id) %>%
    mutate(pmin = pmin/p0) %>% 
    summarise(dpl=quantile(final_dpl,0.05),
              pmin=quantile(pmin,0.05)) %>% 
    ungroup() %>% 
    mutate(hypo = as.numeric(gsub('NF-.([0-9])-.','\\1',ref)),
           type = gsub('NF-([A-Z]).-.','\\1',ref)) %>% 
    left_join(thr %>% select(-ref) %>% collect()) %>% 
    mutate(uab_fin = ifelse(dpl > dpl_72,'A',
                            ifelse(dpl > dpl_60,'B','U')),
           uab_min = ifelse(pmin > dpl_72,'A',
                            ifelse(pmin > dpl_60,'B','U')),
           uab_xcomb = ifelse(uab_fin == 'A' | uab_min == 'A','A',
                             ifelse(uab_fin == 'B' | uab_min == 'B','B','U')),
           pop_id = ordered(pop_id,levels= c('W','C1','C2','C3','E','S'))) %>% 
    gather(uab_stat,uab,uab_fin:uab_xcomb) %>%
    arrange(pop_id) %>% 
    group_by(ref,variant,uab_stat) %>% 
    summarise(uab = paste(uab,collapse=' ')) %>% 
    spread(uab_stat,uab) %>% 
    mutate(combined=ifelse(grepl('U',uab_xcomb),'U',
                           ifelse(grepl('B',uab_xcomb),'B','A')))
    
  area.res <- 
    res %>% 
    filter(pop_type =='area') %>% 
    collect() %>% 
    group_by(ref,variant,trial) %>%
    ## calculate catch statistics (- aboriginal catches)
    summarise(cf10 = sum(cf10)-19,
              cl10 = sum(cl10)-19,
              total_catch = sum(total_catch)-1900,
              avg_catch = sum(avg_catch)-19) %>% 
    ungroup() %>% 
    group_by(ref,variant) %>% 
    summarise(cf10 = median(cf10),
              cl10 = median(cl10),
              total_catch = median(total_catch),
              avg_catch = median(avg_catch)) 
  
  final_res <- 
    full_join(area.res,pop.res) %>% 
    arrange(ref,variant)
  return(final_res)
    
}