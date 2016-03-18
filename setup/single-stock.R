ManSetup <- function(db_name='druna.db',dir='outn',
                     template_file = 'setup/single_template.txt'){
  db <- src_sqlite(db_name)
  dir.create(dir)
  file.copy('data/random.num',sprintf('%s/random.num',dir))
  ## create the depletion files for each trial
  pop <- 
    tbl(db,'naf_pop') %>% 
    filter(pop_type =='Mature females', year %in% c(1864,2014)) %>% 
    collect() %>% 
    filter(grepl('NF-..-1',ref)) %>% 
    mutate(year = sprintf('x%s',year)) %>% 
    spread(year,number) %>% 
    mutate(dpl=x2014/x1864)%>% 
    select(ref,trial,pop_id,dpl) %>% 
    mutate(pop_id = ordered(pop_id, levels= c('W','C1','C2','C3','E','S')))

  tmp_func <- function(x){
    write(x,file=sprintf('%s/%s.dpl',dir,x))
    pop %>% 
      filter(ref==x) %>% 
      spread(pop_id,dpl) %>% 
      select(-c(ref,trial)) %>%
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
      infuse(template_file,dplfile=sprintf('%s.dpl',x$ref),dplcol=st) %>% 
        write(file=sprintf('%s/%s-%s.ss',dir,x$ref,st))
    }
  }
  pop %>%
    group_by(ref) %>% 
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
              '-res {{copy}}.res0'),
        paste('man-v12 -main {{copy|copy.dat}}',
              '-clc {{clc|CLC-N.PAR}}',
              '-log {{log|manlog}}',
              '-res {{copy}}.restest'),
        paste('manresv8 -res {{copy}}.restest',
              '-res0 {{copy}}.res0',
              '-traj {{copy}}.traj',
              '-thresh {{copy}}.thresh')) %>%
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


ManRun <- function(dir='outn',search.string='NF-..-1-[0-9].ss+$',
                   clc.60 = '../settings/CLC-N.60',
                   clc.72 = '../settings/CLC-N.72'){
  tmp <- mclapply(list.files(dir, search.string),
                  function(x) ManCall(run_dir=dir,copy=x,clc=clc.60),
                  mc.cores = detectCores(logical = TRUE))
  tmp <- mclapply(list.files(dir, search.string),
                  function(x) ManCall(run_dir=dir,copy=x,clc=clc.72),
                  mc.cores = detectCores(logical = TRUE))
}
    
ManResults <- function(db_name='trials.db',dir='outn'){
  db <- src_sqlite(db_name)
  stock.names <- c('W','C1','C2','C3','E','S')
  stock.names.78 <- c('W','C1','C2','E','S')
  stock.names.6 <- c('W','C1','C2','C3','S')
  
  list.files(dir,pattern='*.thresh') %>% 
    map(~safely(read.table)(file=sprintf('%s/%s',dir,.))) %>%
    map('result') %>% 
    bind_rows() %>% 
    select(ref=V1,pop_id=V2,dpl_60=V3,dpl_72=V4) %>% 
    mutate(hypo = as.numeric(gsub('NF-.([0-9])-.','\\1',ref)),
           type = gsub('NF-([A-Z]).-.','\\1',ref),
           pop_id = ifelse(hypo < 6,stock.names[pop_id],
                           ifelse(hypo==6,stock.names.6[pop_id],
                                  stock.names.78[pop_id]))) %>% 
    as.data.frame() %>% 
    copy_to(db,df=.,name='man_thresh',temporary=FALSE)
}

