library(dplyr)
library(infuser)
library(purrr)
library(parallel)
library(RSQLite)
library(readr)
library(tidyr)
library(stringi)

#' Title
#'
#' @param hypo 
#' @param msyr 
#' @param trialtype 
#' @param template_file 
#'
#' @return
#' @export
#'
#' @examples
NafTrial <- function(hypo=1,msyr=0.01,trialtype='B',
                     template_file='setup/trial_template.txt',nf='NF',
                     ...){
  
  ## default values for variables
  nstk <- 6
  nsuba <- 7
  ngamma <- 6
  ndelta <- 2
  initdelta <- '0.15   0.02'
  initgamma <- 1.0
  ctype <- 'B'
  nummix <- 1
  mixmat2 <- ' '
  
  ## deviations  
  ## Fewer stocks in hypotheses 6 to 8
  if(hypo %in% c(6,7,8)){
    nstk <- 5
  }
  
  ## mixing parameter estimated between EI/F and SP in hypo 6
  if(hypo == 6){
    initgamma <- 0.4
  }
  
  ## WI and EG are merged
  if(hypo %in% c(7,8)){
    nsuba <- 6
  }
  
  if(hypo == 7){
    ndelta <- 1
    initgamma <- 0.4
  }
  ## mixing used instead of dispersion in hypotheses 4 and 8
  if(hypo %in% c(4,8)){
    ndelta <- 0
    initdelta <- ''
    initgamma <- 0.2
  }
  
  ## High catch series trials
  if(trialtype == 'H'){
    ctype <- 'H'
  }
  
  
  ## mixing matricies
  if(hypo==1){ 
    mixmat <- 
"1    1          1   1   0   0    0   0   0    1  4  0  0  0  0  0
2    1          0   0   1   0    0   0   0    0  0  0  0  0  0  0
3    1          0   0   0   1    0   0   0    0  0  0  0  0  0  0 {{mixmat2}}
4    1          0   0   0   0    1   0   0    0  0  0  0  0  0  0
5    1          0   0   0   0    0   1   0    0  0  0  0  0  0  0
6    1          0   0   0   0    0   0   1    0  0  0  0  0  0  0"
    
  } else if(hypo==2) {
    mixmat <-
"1    1         .88 .88 .1  .02   0   0   0    1  4  0  0  0  0  0 
2    1          0   0   1   0    0   0   0    0  0  0  0  0  0  0
3    1          0   0   0   1    0   0   0    0  0  0  0  0  0  0 {{mixmat2}}
4    1          0   0   0   0    1   0   0    0  0  0  0  0  0  0
5    1          0   0   0  .02  .1  .88  0    0  0  0  0  0  0  0
6    1          0   0   0   0    0   0   1    0  0  0  0  0  0  0"
  } else if(hypo==3) {
    mixmat <- 
"1    1          1   1   0   0   0   0   0     1  4  0  0  0  0  0
2    1         .1  .1  .9   0   0   0   0     1  4  0  0  0  0  0
3    1          0   0   0   1   0   0   0     0  0  0  0  0  0  0 {{mixmat2}}
4    1          0   0   0   0  .9  .1   0     0  0  0  0  0  0  0
5    1          0   0   0   0   0   1   0     0  0  0  0  0  0  0
6    1          0   0   0   0   0   0   1     0  0  0  0  0  0  0"
  } else if(hypo == 4){
    mixmat <- 
"1    1          1   1   0   0   0   0   0     1  4  0  0  0  0  0
2    1          1   1   1   1   0   0   0     2  6  5  3  0  0  0
3    1          0   0   1   1   1   0   0     0  0  3  5  3  0  0 {{mixmat2}}
4    1          0   0   0   1   1   1   0     0  0  0  3  5  3  0
5    1          0   0   0   0   0   1   0     0  0  0  0  0  0  0
6    1          0   0   0   0   0   0   1     0  0  0  0  0  0  0"
  } else if(hypo == 5){
    mixmat <- 
"1    1          1   1   0   0   0   0   0     1  4  0  0  0  0  0
2    1          0   0   1   0   0   0   0     0  0  0  0  0  0  0
3    1          0   0   0   1   0   0   0     0  0  0  0  0  0  0 {{mixmat2}}
4    1          0   0   0   0   1   0   0     0  0  0  0  0  0  0
5    1          0   0   0   0   0   1   0     0  0  0  0  0  0  0
6    1          0   0   0   0  .02 .1  .88    0  0  0  0  0  0  0"
  } else if(hypo == 6){
    mixmat <- 
"1    1         .88 .88 .1  .02   0   0   0    1  4  0  0  0  0  0
2    1          0   0   1   0    0   0   0    0  0  0  0  0  0  0
3    1          0   0   0   1    0   0   0    0  0  0  0  0  0  0 {{mixmat2}}
4    1          0   0   0   0    1   1   0    0  0  0  0  3  2  0
5    1          0   0   0   0    0   0   1    0  0  0  0  0  0  0"
  } else if(hypo == 7){
    mixmat <- 
"1    1          1   1   0   0    0   0        1  4  0  0  0  0
2    1          1   1   1   0    0   0        5  6  3  0  0  0 
3    1          0   0   0  .9   .1   0        0  0  0  0  0  0 {{mixmat2}}
4    1          0   0   0   0    1   0        0  0  0  0  0  0
5    1          0   0   0   0    0   1        0  0  0  0  0  0"
  } else if(hypo == 8){
    mixmat <- 
"1    1          1   1   0   0    0   0        1  4  0  0  0  0 
2    1          1   1   1   1    0   0        2  6  5  3  0  0 
3    1          0   0   1   1    1   0        0  0  3  5  3  0 {{mixmat2}}
4    1          0   0   0   0    1   0        0  0  0  0  0  0 
5    1          0   0   0   0    0   1        0  0  0  0  0  0"
  }

  ## Change in mixing matricies
  if(trialtype %in% c('G','F')){
    nummix <- 2
    mixmat2 <- '\n3   2       0  0   .3  .7  0  0  0    0  0  0  0  0  0  0'
  } else {
    mixmat2 <- ''
  }
  mixmat <- infuse(mixmat,mixmat2=mixmat2)
  
    
  ## trial code generated
  ref <- '{{nf}}-{{trialtype}}{{hypo}}-{{msyr}} H{{hypo}}-{{msyr}}' %>%
    infuse(hypo=hypo,msyr=100*msyr,trialtype=trialtype,nf=nf)
  
  ## assign msyr rates to nstk stocks
  msyr.txt <- 
    sprintf('MSYR OF STOCK 1 (MSYR%1$s(%1$s))                   %2$s',
            1:nstk,msyr) %>%
    paste(collapse='\n')
  
  ## change the density dependance from 1+ for 1% to mature for 4%
  if(msyr==0.01){
    optf <- 0
  } else {
    optf <- 1
  }
  
  ## produce the trial setup
  txt <- infuse(template_file,
                c(list(ref=ref,
                       msyrtxt=msyr.txt,
                       optf=optf,
                       nstk=nstk,
                       mixmat=mixmat,
                       ndelta=ndelta,
                       nsuba = nsuba,
                       initdelta = initdelta,
                       initgamma = initgamma,
                       disp=ifelse(ndelta==0,0,1),
                       disp2=ifelse(ndelta==0,'0',ifelse(hypo==7,'','1')),
                       ctype = ctype,
                       nummix = nummix),
                  list(...)))
  attr(txt,'ref') <- '{{nf}}-{{trialtype}}{{hypo}}-{{msyr}}' %>%
    infuse(hypo=hypo,msyr=100*msyr,trialtype=trialtype,nf=nf)
  class(txt) <- c('NafTrial',class(txt))
  return(txt)
}

#' Title
#'
#' @param obj 
#' @param dir 
#'
#' @return
#' @export
#'
#' @examples
NafWrite <- function(obj,dir='.'){
  write(as.character(obj),file = sprintf('%s/%s.dat',dir,attr(obj,'ref')))
}

#' Run naf-ist
#'
#' Wrapper function that allows the naf-ist program to be called directly from R.
#' 
#' @param ... list of variables used in the run
#'
#' @return NULL
#' @export
#'
#' @examples
#' callNaf()
NafCall <- function(...){
  tmp <- list(...)
  if('run_dir' %in% names(tmp)){
    old.dir <- getwd()
    setwd(tmp$run_dir)
  }
  run.string <- 'naf-ist -main {{copyna|copyna.dat}} -par {{nafpar|nafpar.par}}' %>%
    infuse(tmp)
  res <- tryCatch(system(run.string,
                         ignore.stdout = TRUE,
                         ignore.stderr = FALSE,
                         intern = TRUE),
                  error=function(e) sprintf('Trial %s was not completed',tmp$copyna))
  if('run_dir' %in% names(tmp)){
    setwd(old.dir)
  }
}

NafWeights <- function(dir = 'weights_trials'){
  dir.create(dir,showWarnings = FALSE)
  tmp <- 
    list.files('data') %>% 
    map(~file.copy(from=sprintf('data/%s',.),
                   to=sprintf('%s/%s',dir,.),
                   overwrite = TRUE))
  tmp <- file.copy('settings/manage-single.dat',sprintf('%s/manage.dat',dir))
  for(hypo in 1:8){
    for(msyr in c(0.01,0.04)){
      for(tagwt in 0.1*(1:9)){
        for(wtagei in 0.01*(1:9)){
            NafWrite(NafTrial(hypo = hypo,msyr = msyr,
                              tagwt=tagwt,wtagei = wtagei,
                              nf= sprintf('T%s',10*tagwt),
                              trialtype = 100*wtagei),dir)
        }
      }
    }
  }
}
  
NafSetup <- function(dir = 'trials'){
  dir.create(dir,showWarnings = FALSE)
  tmp <- 
    list.files('data') %>% 
    map(~file.copy(from=sprintf('data/%s',.),
                   to=sprintf('%s/%s',dir,.),
                   overwrite = TRUE))
  tmp <- file.copy('settings/manage.dat',sprintf('%s/manage.dat',dir))
  ## NF-B trials (baseline)
  for(hypo in 1:8){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo,msyr = msyr),dir)
    }
  }
  ## NF-H trials (high catch series)
  for(hypo in c(1,3,4)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='H'),dir)
    }
  }
  
  ## NF-T trials (tag loss 20% in y1, 10 % there afther)
  for(hypo in c(1,3,4)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='T', tloss1 = 0.2, 
                        tloss2 = 0.1),dir)
    }
  }
  
  ## NF-W1 (tags upweighted x10)
  for(msyr in c(0.01,0.04)){
    NafWrite(NafTrial(hypo = 1, msyr = msyr,
                      trialtype='W', tagwgt=10),dir)
  }
  
  ## NF-S trials (selectivity before and after 2007 estimated)
  for(hypo in c(3,4)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='S', selyr = 2007),dir)
    }
  }
  
  ## NF-U1 (selectivy decrease)
  for(msyr in c(0.01,0.04)){
    NafWrite(NafTrial(hypo = 1, msyr = msyr,
                      trialtype='U', seldec = 0.04),dir)
  }
  
  
  ## NF-R trials (tags excluded after year 1)
  for(hypo in c(3,4)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='R', 
                        tagexc1 = 1),dir)
    }
  }

  ## NF-A trials (pro-rate abundance)
  for(hypo in c(3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='A', 
                        optsgt = 1),dir)
    }
  }
  
  ## NF-G trials changes in mixing matricies
  for(hypo in c(1,3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='G', optmix = 2),dir)
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='F', optmix = 3),dir)
    }
  }
  ## NF-C trials adds a cpue likelihood 
  for(msyr in c(0.01,0.04)){
    NafWrite(NafTrial(hypo = 3, msyr = msyr,
                      trialtype='C', optcpe = 1),dir)
  }

  
  
    
  ## TODO: 
  ## - NF-X, NF-C (should we skip them?)
  ## - Missing data: NF-H, NF-X NF-P
  ## - Future trials: NF-Q, A, Y, F 
  
}

NafCond <- function(dir='trials',
                    nafpar='../variants/Naf-v0.par',
                    search.string = 'NF-[A-Z][0-9]-[0-9].dat'){
  mclapply(list.files(dir, search.string),
           function(x) NafCall(run_dir=dir,copyna=x,nafpar=nafpar),
           mc.cores = detectCores(logical = TRUE))
} 


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

NafResults.readTag <- function(file='NAF.pop', trials_db){
  tag <- tryCatch(read.table(file = file, skip = 3),
                  error=function(e) NULL)
  header <- readLines(file,n=3)
  if(is.null(tag)){
    warning(sprintf('File %s could not be read',file))    
    return(tag)
  }
  
  ref <- gsub('(NF-[A-Z][0-9]-[0-9]).+','\\1',header[1],perl=TRUE)
  tmp <- tolower(scan(text=gsub('Obs:|Pred:|\\+|\\*','',header[2]),
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

