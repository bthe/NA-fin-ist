NafResults <- function(dir='trials', db_name = 'trials.db',
                       search.string='NF-[A-Z][0-9]-[0-9]'){
  ## define an sqllite database to digest the wealth of output
  trials_db <- src_sqlite(db_name, create = T)
  print('Read in the fit to age')
  #  db_drop_table(trials_db$con,'naf_age')
  db_create_table(trials_db$con,'naf_age',
                  types = c(year = 'INTEGER',area = 'INTEGER',age = 'INTEGER',
                            obs.fem = 'REAL',prd.fem = 'REAL',obs.m = 'REAL',
                            prd.m = 'REAL',trial = 'INTEGER',ref = 'TEXT') )
  
  list.files(path = dir,pattern = sprintf('%s.age',search.string)) %>%
    map(~NafResults.readAge(sprintf('%s/%s',dir,.),
                            trials_db))
  db_create_index(trials_db$con, 'naf_age',
                  c('year','trial','age'))
  
  print('Read in the condition variables')
  #db_drop_table(trials_db$con,'naf_all')
  db_create_table(trials_db$con,'naf_all',
                  types=c(ref='TEXT',trial = 'INTEGER',constant='TEXT',
                          value = 'REAL'))
  list.files(path = dir,pattern = sprintf('%s.all',search.string)) %>%
    map(~NafResults.readAll(sprintf('%s/%s',dir,.),
                            trials_db))
  db_create_index(trials_db$con, 'naf_all','trial')
  
  print('Read in the catch ouput')
  #db_drop_table(trials_db$con,'naf_cat')
  db_create_table(trials_db$con,'naf_cat',
                  types = c(year = 'INTEGER',catt1 = 'INTEGER',catt2 = 'INTEGER',
                            catt3 = 'INTEGER',catt4 = 'INTEGER',catt5 = 'INTEGER',
                            catt6 = 'INTEGER',ck1 = 'INTEGER',ck2 = 'INTEGER',
                            ck3 = 'INTEGER',ck4 = 'INTEGER',ck5 = 'INTEGER',ck6 = 'INTEGER',
                            ck7 = 'INTEGER',cpk1 = 'INTEGER',cpk2 = 'INTEGER',
                            cpk3 = 'INTEGER',cpk4 = 'INTEGER',cpk5 = 'INTEGER',
                            cpk6 = 'INTEGER',cpk7 = 'INTEGER',ref = 'TEXT',
                            variant = 'TEXT',trial = 'INTEGER'))
  list.files(path = dir,pattern = sprintf('%s.cat',search.string)) %>%
    map(~NafResults.readCat(sprintf('%s/%s',dir,.),
                            trials_db))
  db_create_index(trials_db$con, 'naf_cat',c('trial','variant','year'))
  
  print('read in the population numbers')
  #db_drop_table(trials_db$con,'naf_pop')
  db_create_table(trials_db$con,'naf_pop',
                  types = c(year = 'INTEGER',ref = 'TEXT',variant = 'TEXT',
                            hypo = 'INTEGER', msyr = 'INTEGER',
                            trial = 'INTEGER',pop_type = 'TEXT',
                            pop_id = 'TEXT',number = 'INTEGER'))
  
  list.files(path = dir,pattern = sprintf('%s.pop',search.string)) %>%
    map(~NafResults.readPop(sprintf('%s/%s',dir,.),
                            trials_db))
  db_create_index(trials_db$con, 'naf_pop',c('ref','trial','variant','year'))
  
  print('Read tagging results')
  #db_drop_table(trials_db$con,'naf_tag')
  db_create_table(trials_db$con,'naf_tag',
                  types = c(rela = 'TEXT',year = 'INTEGER',ref = 'TEXT',
                            area = 'TEXT',obs = 'REAL',prd = 'REAL'))
  list.files(path = dir,pattern = 'NF-[A-Z][0-9]-[0-9].tag') %>%
    map(~NafResults.readTag(sprintf('%s/%s',dir,.),
                            trials_db))
  db_create_index(trials_db$con, 'naf_tag',c('ref','year'))
  
  print('Read survey data') 
  NafResults.readSurvey(file = sprintf('%s/surveyn.dat',dir),trials_db)
  
}


NafResults.readAge <- function(file='NAF.age', trials_db){
  
  age <- readLines(file)  
  if(length(age)>1){ 
    loc <- grep('YEAR',age)
    ## find the trial reference
    ref <- gsub('(..-.[0-9]-[0-9]).+','\\1',age[1],perl=TRUE)
    ## add conditioning trial number
    age <- paste(age[-c(1,loc[-1])],
                 c('trial',rep(0:(length(loc)-1),each = (diff(loc)[1]-1))))
    ## throw away the zero entries
    age <- age[!grepl('0.0     0.0     0.0     0.0',age)]
    age[1] <- tolower(age[1])
    tmp <- read.table(text = paste(age,collapse = '\n'), header = TRUE) %>%
      mutate(ref=ref) %>%
      db_insert_into(trials_db$con,table='naf_age',values = .)
  } else {
    warning(sprintf('File %s could not be read',file))    
  }
  
}

NafResults.readAll <- function(file='NAF.all', trials_db){
  all <- readLines(file)
  if(length(all)>2){
    ## find the trial reference
    ref <- gsub('(..-.[0-9]-[0-9]).+','\\1',all[1],perl=TRUE)
    all <- all[-1]
    ## remove uneccessary spaces in variable names
    all[1] <- gsub('+','',gsub('([A-Z][a-z]+) ([0-9])','\\1\\2',all[1]),fixed=TRUE)
    ## add selectivity parameter names to the header
    all[1] <- tolower(paste(all[1],'RM50 RMSIG RF50 RFSIG RM507 RMSIG7 RF507 RFSIG7'))
    ## if hypo 6-8 stuff needs to be added
    hypo <- as.numeric(gsub('..-.([0-9])-[0-9]','\\1',ref))
    if(hypo>5){
      all[1] <- paste(all[1],'kmat6 depn6')
    }
    tmp <- read.table(text=paste(all,collapse = '\n'),header=TRUE,fill = TRUE) %>%
      mutate(ref=ref) %>%
      rename(trial=n) %>%
      gather(constant,value,-c(ref,trial)) %>%
      db_insert_into(trials_db$con,table='naf_all',values = .)
    
  } else {
    warning(sprintf('File %s could not be read',file))    
  }
}

NafResults.readCat <- function(file='NAF.all', trials_db){
  catch <- tryCatch(read.table(file = file, skip = 2),
                    error = function(e) NULL)
  if(is.null(catch)){
    warning(sprintf('File %s could not be read',file))    
    return(catch)
  }
  header <- readLines(file,n=2)
  ref <- gsub('(..-.[0-9]-[0-9]).+','\\1',header[1],perl=TRUE)
  variant <- gsub('.+(V[0-9]).+','\\1',header[1])
  hypo <- as.numeric(gsub('..-.([0-9])-[0-9]','\\1',ref))
  nsuba <- 7
  nstk <- 6
  if(hypo>6){
    nsuba <- 6  
  } 
  if(hypo>5){
    nstk <- 5
  }
  
  names(catch) <- c('year',paste0('catt',1:nstk),
                    paste0('ck',1:nsuba),
                    paste0('cpk',1:nsuba))
  if(hypo>5){
    catch <- mutate(catch,catt6 = NA)
  }
  
  if(hypo>6){
    catch <- mutate(catch,ck7 = NA,cpk7=NA)
  }
  
  class(catch) <- 'data.frame'
  
  catch %>%
    mutate(ref = ref,
           variant = variant,
           trial = rep(0:(sum(year==min(year))-1),
                       each = diff(range(year))+1)) %>% 
    
    db_insert_into(trials_db$con,table='naf_cat',values = .)
}

NafResults.readPop <- function(file='NAF.pop', trials_db){
  pop <- tryCatch(read.table(file = file, skip = 2),
                  error = function(e) NULL)
  if(is.null(pop)){
    warning(sprintf('File %s could not be read',file))    
    return(pop)
  }
  
  header <- readLines(file,n=2)
  ref <- gsub('(..-.[0-9]-[0-9]).+','\\1',header[1],perl=TRUE)
  variant <- gsub('.+(V[0-9]).+','\\1',header[1])
  hypo <- gsub('..-.([0-9]).+','\\1',ref) %>% as.integer()
  msyr <- gsub('..-..-([0-9])','\\1',ref) %>% as.integer()
  
  names(pop) <- tolower(scan(text=gsub(': ','.',header[2]),
                             what='character', quiet=TRUE))
  
  class(pop) <- 'data.frame'
  
  stock.names <- c('W','C1','C2','C3','E','S')
  area.names <- c('EC','WG','EG','WI','EI/F','N','SP')
  if(hypo %in% 7:8){
    stock.names <- c('W','C1','C2','E','S')
    area.names <- c('EC','WG','EG+WI','EI/F','N','SP')
  }
  if(hypo == 6){
    stock.names <- c('W','C1','C2','C3','S')
  }
  
  pop %>%
    rename(year=yr) %>%
    mutate(ref = ref,
           variant = variant,
           hypo = hypo,
           msyr = msyr,
           trial = rep(0:(sum(year==min(year))-1),each = diff(range(year))+1)) %>%
    gather(population,number,-c(year,ref,variant,msyr,hypo,trial)) %>%
    separate(population,c('pop_type','pop_id')) %>%
    mutate(pop_type = ifelse(pop_type=='pf','Mature females','1+'),
           pop_id = ifelse(pop_type == '1+',area.names[as.numeric(pop_id)],
                           stock.names[as.numeric(pop_id)])) %>%
    db_insert_into(trials_db$con,table='naf_pop',values = .)
  
}

NafResults.readTag <- function(file='NAF.tag', trials_db){
  tag <- tryCatch(read.table(file = file, skip = 3),
                  error=function(e) NULL)
  header <- readLines(file,n=3)
  if(is.null(tag)){
    warning(sprintf('File %s could not be read',file))    
    return(tag)
  }
  
  ref <- gsub('(NF-[A-Z][0-9]-[0-9]).+','\\1',header[1],perl=TRUE)
  tmp <- tolower(scan(text=gsub('EGWIEI','EGWI EI',
                                gsub('Obs:|Pred:|\\+|\\*','',header[2])),
                      what='character', quiet=TRUE))
  areas <- unique(tmp[-c(1:2)])
  names(tag) <- c(tmp[1:2],paste(c(areas,areas),
                                 rep(c('obs','prd'),
                                     each=length(areas)),
                                 sep = '.'))
  class(tag) <- 'data.frame'
  tag %>%
    rename(year = yr) %>%
    mutate(ref = ref, 
           rela = areas[rela]) %>% 
    gather(type,value,-c(year,rela,ref)) %>% 
    separate(type,c('area','type')) %>% 
    spread(type,value) %>%
    db_insert_into(trials_db$con,table='naf_tag',values = .)
}

NafResults.readSurvey <- function(file='surveyn.dat', trials_db){
  sight <- read_lines(file = file, n_max = 22) %>%
    map(~stri_sub(.,to=70)) %>%
    unlist() %>%  stri_paste(collapse='\n') %>% 
    tolower() %>%
    read.table(text=.,skip=2,header=TRUE,stringsAsFactors = FALSE) %>% 
    mutate(area = gsub('/f|n/','',area)) %>%
    select(year=iyr,area, obs = sight,cv=cvcsla,
           pro.obs = prorated) %>%
    copy_to(trials_db,name='naf_sight',df = .,temporary = FALSE)
  
}