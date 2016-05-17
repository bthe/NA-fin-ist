ManSetup <- function(db_name='druna.db',dir='outn',
                     template_file = 'setup/single_template.txt'){
  db <- src_sqlite(db_name)
  dir.create(dir)
  file.copy('data/random.num',sprintf('%s/random.num',dir))
  ## create the depletion files for each trial
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
    mutate(pop_id = ordered(pop_id, levels= c('W','C1','C2','C3','E','S')))

  tmp_func <- function(x){
    write(x,file=sprintf('%s/%s.dpl',dir,x))
    pop %>% 
      filter(ref==x) %>% 
      spread(pop_id,dpl) %>% 
      select(-c(ref,trial,msyr)) %>%
      ## pad the files with unecessary zeros to comply for mantst12 requirements
      mutate(ecol1=0,ecol2=0,ecol=0) %>% 
      write.table(.,file=sprintf('%s/%s.dpl',dir,x),
                  quote = FALSE, col.names =FALSE, row.names =FALSE,
                  append = TRUE)
  }
  unique(pop$ref) %>% 
    map(~tmp_func(.))

  ## create setup files for the single stock trials
  tmp_fuse <- function(x){
    for(st in 1:x$num.stocks){ 
      infuse(template_file,dplfile=sprintf('%s.dpl',x$ref),dplcol=st,
             msyr1=0.01*x$msyr,optf=ifelse(x$msyr==1,1,2)) %>% 
        write(file=sprintf('%s/%s-%s.ss',dir,x$ref,st))
    }
  }
  pop %>%
    group_by(ref,msyr) %>% 
    summarise(num.stocks = n()/101) %>% 
    split(.$ref) %>% 
    map(~tmp_fuse(.)) -> tmp
}


ManCall <- function(...){
  tmp <- list(...)
  if('run_dir' %in% names(tmp)){
    old.dir <- getwd()
    setwd(tmp$run_dir)
  }
  run.string <-
      c(paste('man-v12z -main {{copy|copy.dat}}',
              "-log {{copy|manlog}}.z0",
              '-res {{copy}}.{{clc_desc|72}}.res0'),
        paste('man-v12 -main {{copy|copy.dat}}',
              '-clc {{clc|CLC-N.PAR}}',
              '-log {{log|manlog}}',
              '-res {{copy}}.{{clc_desc|72}}.restest'),
        paste('manresv8 -res {{copy}}.{{clc_desc|72}}.restest',
              '-res0 {{copy}}.{{clc_desc|72}}.res0',
              '-traj {{copy}}.{{clc_desc|72}}.traj',
              '-thresh {{copy}}.{{clc_desc|72}}.thresh')) %>%
        map(~infuse(.,tmp))
  
  res <-
      for(run in run.string){
          tryCatch(system(run,
                          ignore.stdout = TRUE,
                          ignore.stderr = FALSE,
                          intern = TRUE),
                   error=function(e) sprintf('Trial %s was not completed',tmp$copy))
      }
  if('run_dir' %in% names(tmp)){
    setwd(old.dir)
  }
}


ManRun <- function(dir='outn',search.string='NF-..-.-[0-9].ss+$',
                   clc.60 = '../settings/CLC-N.60',
                   clc.72 = '../settings/CLC-N.72'){
  tmp <- mclapply(list.files(dir, search.string),
                  function(x) ManCall(run_dir=dir,copy=x,clc=clc.60,clc_desc=60),
                  mc.cores = detectCores(logical = TRUE))
  tmp <- mclapply(list.files(dir, search.string),
                  function(x) ManCall(run_dir=dir,copy=x,clc=clc.72,clc_desc=72),
                  mc.cores = detectCores(logical = TRUE))
}

ManResults <- function(db_name='trials.db',dir='outn'){
  db <- src_sqlite(db_name)
  list.files(dir,pattern = '*.restest') %>% 
    mclapply(., function(x) safely(ManResults.restest)(sprintf('%s/%s',dir,x)),
             mc.cores = detectCores(logical=TRUE)) %>% 
    map('result') %>% 
    bind_rows() %>% 
    as.data.frame() %>% 
    copy_to(db,df=.,name='man_thresh',temporary=FALSE)
}

ManResults.restest <- function(file='NAF.restest'){
  tmp_func <- function(file){
    res <- readLines(file)
    tuning <- gsub('^.+([0-9]{2}).+$','\\1',file) %>% 
      as.numeric()
    topinfo <- res[1:45]
    res <- res[-c(1:45)]
    ref <- gsub('\\s*(NF-..-.).+','\\1',topinfo[3])
    pop <- 
      gsub('\\s*NF-..-..+([0-9])','\\1',topinfo[3]) %>% 
      as.numeric()
    
    run.para <- 
      read.table(text = res[grepl('New Para',res)]) %>% 
      select(K=V6,init.dpl=V4) %>% 
      mutate(trial = 1:n())
    
    header <- c('year','fem','plus','catch') #scan(text=gsub(': ','.',res[2]),what='character')
    hypo <- as.numeric(gsub('NF-.([0-9])-.','\\1',ref))
    loc <- grep('Trial',res)
    tmp <- gsub('(^\\s*$|^\\s*[a-z].+|^-1)','# \\1',tolower(res)) %>% 
      read.table(text=.,fill=TRUE)
    names(tmp) <- tolower(header)
    if(class(tmp$fem)=='factor')
      tmp$fem <- as.numeric(as.character(tmp$fem))
    
    stock.names <- c('W','C1','C2','C3','E','S')
    if(hypo == 6){
      stock.names <- c('W','C1','C2','C3','S')
    }
    
    
    tmp %>%
      filter(!is.na(fem)) %>% 
      mutate(ref = ref,
             stock = stock.names[pop],
             tuning = tuning,
             trial = cut(1:length(year),c(0,which(diff(year)<0),1e9),labels = FALSE)) %>% 
      left_join(run.para) %>%
      select(-c(plus,catch)) %>% 
      rename(number=fem) %>% 
      group_by(trial) %>% 
      mutate(K = ifelse(abs(number[year==0]/init.dpl-K)>1,number[year==0]/init.dpl,K),
             dpl=number/K)
  }
  
  res <- tmp_func(file)
  res0 <- tmp_func(gsub('restest','res0',file))
  
  res %>% 
    left_join(res0 %>% select(year,stock,tuning,trial,ref,zero_dpl=dpl)) %>% 
    group_by(ref,stock,tuning,trial) %>% 
    summarise(pmin = min(dpl/zero_dpl),pfin=max(dpl[year==100])) %>% 
    group_by(ref,stock,tuning) %>% 
    summarise(pmin5 = quantile(pmin,0.05,na.rm=TRUE),
              pfin5 = quantile(pfin,0.05,na.rm=TRUE))
}