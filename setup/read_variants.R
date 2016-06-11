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

  db_create_table(db$con,'naf_resbyvariant',
                  types = c(ref = 'TEXT',variant = 'TEXT',
                            trial = 'INTEGER',
                            total_catch = 'REAL',
                            avg_catch = 'REAL',
                            cf10 = 'REAL',cl10 = 'REAL'))
  
  db_create_table(db$con,'naf_resbyyear',
                  types = c(ref = 'TEXT',variant = 'TEXT',
                            year = 'INTEGER',pop_type = 'TEXT',
                            pop_id = 'TEXT',
                            catch_lower = 'REAL',
                            catch_upper = 'REAL',
                            catch_med = 'REAL',
                            abundance_lower = 'REAL',
                            abundance_upper = 'REAL',
                            abundance_med = 'REAL',
                            dpl_lower = 'REAL',
                            dpl_upper = 'REAL',
                            dpl_med = 'REAL'))
  
  res <- 
    list.files(path = dir,pattern =search.string)%>%
    mclapply(., function(x) safely(NafReadVariants.restest)(sprintf('%s/%s',dir,x)),
             mc.cores = detectCores(logical=TRUE)) 
  
  res %>% 
    map('result') %>% 
    map('res.total') %>% 
    bind_rows() %>% 
    as.data.frame() %>% 
    db_insert_into(db$con,table='naf_res',values = .)
  
  res %>% 
    map('result') %>% 
    map('res.by.year') %>% 
    bind_rows() %>% 
    as.data.frame() %>% 
    db_insert_into(db$con,table='naf_resbyyear',values = .)
  
  res %>% 
    map('result') %>% 
    map('res.by.variant') %>% 
    bind_rows() %>% 
    as.data.frame() %>% 
    db_insert_into(db$con,table='naf_resbyvariant',values = .)
  
}

NafReadVariants.restest <- function(file='NAF.restest'){
  print(sprintf('Reading file %s',file))
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
  
#  res <- tmp_func(file)
#  res0 <- tmp_func(gsub('V[0-9]','V0',file))
  tmp2 <- 
    tmp_func(file) %>% 
    left_join(tmp_func(gsub('V[0-9]','V0',file)) %>% 
                select(ref,trial,year,pop_type,pop_id,zero_abundance=abundance))
  
  res.total <- 
    tmp2 %>% 
    group_by(ref,variant,trial,pop_type,pop_id) %>% 
    summarise(total_catch=sum(catch[year>2014]),avg_catch=mean(catch[year>2014]), 
              p0 = sum(abundance[year==1865]),pfin = sum(abundance[year==2115]),
              pmin = min(abundance/zero_abundance), cf10=mean(catch[year %in% 2015:2024]),
              cl10 = mean(catch[year %in% 2105:2114]),
              final_dpl = pfin/p0) 
  
  res.by.variant <- 
    tmp2 %>% 
    filter(pop_type=='area') %>% 
    group_by(ref,variant,trial) %>% 
    summarise(total_catch=sum(catch[year>2014]),
              avg_catch=mean(catch[year>2014]), 
              cf10=mean(catch[year %in% 2015:2024]),
              cl10 = mean(catch[year %in% 2105:2114])) 
  
  res.by.year <- 
    tmp2 %>% 
    group_by(ref,variant,pop_type,pop_id) %>%
    mutate(dpl = pmin(abundance/abundance[year==1865],1.0)) %>% 
    group_by(ref,variant,year,pop_type,pop_id) %>%
    summarise(catch_lower=quantile(catch,0.05),
              catch_upper=quantile(catch,0.95),
              catch_med = median(catch),
              abundance_lower=quantile(abundance,0.05),
              abundance_upper=quantile(abundance,0.95),
              abundance_med = median(abundance),
              dpl_lower=quantile(dpl,0.05),
              dpl_upper=quantile(dpl,0.95),
              dpl_med = median(dpl))
  

  return(list(res.total=res.total,res.by.year=res.by.year,
              res.by.variant = res.by.variant))
}


NafPerformance <- function(db,db_res){
  pop <- 
    tbl(db,'naf_pop') %>% 
    filter(pop_type =='Mature females', year %in% c(1865,2014)) %>% 
    collect(n=Inf) %>% 
    distinct(ref,year,trial,pop_type,pop_id,.keep_all=TRUE) %>% 
    mutate(year = sprintf('x%s',year)) %>% 
    spread(year,number) %>% 
    #mutate(dpl=x2014/x1865)%>% 
    mutate(dpl=pmin(x2014/x1865,0.99))%>% 
    select(ref,trial,pop_id,dpl,msyr) %>% 
    mutate(pop_id = ordered(pop_id, levels= c('W','C1','C2','C3','E','S'))) %>% 
    group_by(ref,pop_id) %>% 
    summarise(init_dpl_lower=quantile(dpl,0.05),
              init_dpl_med=quantile(dpl,0.5),
              inti_dpl_upper=quantile(dpl,0.95))
  
  res <- tbl(db_res,'naf_res')
  res.by.variant <- 
    tbl(db_res,'naf_resbyvariant') %>% 
    collect(n=Inf) %>% 
    mutate(avg_catch=7*avg_catch,
           cf10 = 7*cf10,
           cl10 = 7*cl10) %>% 
    group_by(ref,variant) %>% 
    summarise(catch_med = median(total_catch),
              cf10_lower = round(quantile(cf10,0.05)),
              cf10_med = round(quantile(cf10,0.5)),
              cf10_upper = round(quantile(cf10,0.95)),
              cl10_lower = round(quantile(cl10,0.05)),
              cl10_med = round(quantile(cl10,0.5)),
              cl10_upper = round(quantile(cl10,0.95)))
  
  thr <- tbl(db_res,'man_thresh') %>% 
    rename(clc=tuning,dplfin=pfin5,dplmin=pmin5,pop_id=stock) %>% 
    collect(n=Inf) %>% 
    gather(dpl_stat,dpl,dplfin:dplmin) %>% 
    unite(boundary,c(dpl_stat,clc)) %>%
    select(ref,pop_id,boundary,dpl) %>% 
    spread(boundary,dpl)
    
  pop.res <- 
    res %>% 
    filter(pop_type =='pop') %>% 
    collect(n=Inf) %>%
    group_by(ref,variant,pop_id) %>%
    summarise(final_dpl_lower=quantile(final_dpl,0.05),
              final_dpl_med=quantile(final_dpl,0.5),
              final_dpl_upper=quantile(final_dpl,0.95),
              pmin_lower=quantile(pmin,0.05),
              pmin_med=quantile(pmin,0.5),
              pmin_upper=quantile(pmin,0.95)) %>% 
    left_join(pop) %>% 
    left_join(thr) %>% 
    mutate(uab_fin = ifelse(final_dpl_lower > dplfin_72,'A',
                            ifelse(final_dpl_lower > dplfin_60,'B','U')),
           uab_min = ifelse(pmin_lower > dplmin_72,'A',
                            ifelse(pmin_lower > dplmin_60,'B','U')),
           uab_xcomb = ifelse(uab_fin == 'A' | uab_min == 'A','A',
                              ifelse(uab_fin == 'B' | uab_min == 'B','B','U')),
           pop_id = ordered(pop_id,levels= c('W','C1','C2','C3','E','S')))
  
  
  uab_stat <- 
    pop.res  %>% 
    gather(uab_stat,uab,uab_fin:uab_xcomb) %>%
    arrange(pop_id) %>% 
    group_by(ref,variant,uab_stat) %>% 
    summarise(uab = paste(uab,collapse=' ')) %>% 
    spread(uab_stat,uab) %>% 
    mutate(combined=ifelse(grepl('U',uab_xcomb),'U',
                           ifelse(grepl('B',uab_xcomb),'B','A')))
  
  uab_stat <- 
    res.by.variant %>% 
    left_join(uab_stat)
    
#   area.res <- 
#     res %>% 
#     filter(pop_type =='area') %>% 
#     collect(n=Inf) %>% 
#     group_by(ref,variant,trial) %>%
#     ## calculate catch statistics (- aboriginal catches)
#     summarise(cf10 = sum(cf10)-19,
#               cl10 = sum(cl10)-19,
#               total_catch = sum(total_catch)-1900,
#               avg_catch = sum(avg_catch)-19) %>% 
#     ungroup() %>% 
#     group_by(ref,variant) %>% 
#     summarise(cf10 = median(cf10),
#               cl10 = median(cl10),
#               total_catch = median(total_catch),
#               avg_catch = median(avg_catch)) 
  
  final_res <- 
    full_join(area.res,uab_stat) %>% 
    arrange(ref,variant)
  return(list(final_res=uab_stat,pop.res=pop.res))
    
}