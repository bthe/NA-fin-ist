library(dplyr)
library(infuser)
library(purrr)
library(parallel)
library(RSQLite)
library(readr)
library(tidyr)
library(stringi)
library(ggplot2)

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
                     disp1 = 1,
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
  initgamma2 <- ' '
  
  ## deviations  
  ## Fewer stocks in hypotheses 6 to 8
  if(hypo %in% c(6,7,8)){
    nstk <- 5
  }
  
  ## mixing parameter estimated between EI/F and SP in hypo 6
  if(hypo %in% c(5,6)){
    initgamma <- 0.4
  }
  
  ## WI and EG are merged
  if(hypo %in% c(7,8)){
    nsuba <- 6
  }
  
  if(hypo == 7){
    ndelta <- 1
    initgamma <- 0.4
    ngamma <- 8
    initgamma2 <- '1.  .35'
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
4    1          0   0   0   0   1   1   0     0  0  0  0  3  2  0
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
6    1          0   0   0   0  .23 .23  .77   0  0  0  0  3  2  0"
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
3    1          0   0   0   1    1   0        0  0  0  7  8  0 {{mixmat2}}
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
  ref <- 
    '{{nf}}-{{trialtype}}{{hypo}}-{{msyr}} H{{hypo}}-{{msyr}}' %>%
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
                       initgamma2 = initgamma2,
                       disp=ifelse(ndelta==0,0,disp1),
                       disp2=ifelse(ndelta==0,'0',ifelse(hypo==7,'','1')),
                       ctype = ctype,
                       nummix = nummix,
                       ngamma = ngamma,
                       wtagei = ifelse(msyr==0.01,0.1,1)),
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
  run.string <- 
    paste('{{nafist|naf-ist}} -main {{copyna|copyna.dat}}',
          '-par {{nafpar|nafpar.par}}',
          '-man {{man|manage.dat}}',
          '-con {{con|nafcon.dat}}',
          '-survey {{survey|surveyn.dat}}') %>%
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
  
NafSetup <- function(dir = 'trials',hypos=c(1,2,3,5,6)){
  dir.create(dir,showWarnings = FALSE)
  tmp <- 
    list.files('data') %>% 
    map(~file.copy(from=sprintf('data/%s',.),
                   to=sprintf('%s/%s',dir,.),
                   overwrite = TRUE))
  tmp <- file.copy('settings/manage.dat',sprintf('%s/manage.dat',dir))
  tmp <- file.copy('settings/CLC-N.PAR',sprintf('%s/clc-n.par',dir))
  
  ## NF-B trials (baseline)
  for(hypo in hypos){
    for(msyr in c(0.01,0.04)){
      if(hypo == 6 & msyr ==0.01){
        next
      } else {
        NafWrite(NafTrial(hypo = hypo,msyr = msyr),dir)
        tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-B%s-%s.survey',dir,hypo,100*msyr))
      }
    }
  } 
  
  ## NF-H trials (high catch series)
  for(hypo in c(2,3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='H'),dir)
      tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-H%s-%s.survey',dir,hypo,100*msyr))
    }
  }
  
  ## NF-Q trials (a copy of NF-B3 trials with future surveys excluding 60 deg north)
  for(msyr in c(0.01,0.04)){
    NafWrite(NafTrial(hypo = 3, msyr = msyr,
                      trialtype='Q',sratio1 = 0.78, sratio2 =0.93,
                      vratio1 = 0.162, vratio2 =0.085),dir)
    tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-Q%s-%s.survey',dir,3,100*msyr))
  }
  
  ## NF-A trials (pro-rate abundance)
  for(hypo in c(2,3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='A', 
                        optsgt = 1),dir)
      tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-A%s-%s.survey',dir,hypo,100*msyr))
    }
  }
  
  ## NF-U3 (selectivy decrease)
  for(msyr in c(0.01,0.04)){
    NafWrite(NafTrial(hypo = 3, msyr = msyr,
                      trialtype='U', seldec = 0.04,
                      mort3=0.04),dir)
    tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-U%s-%s.survey',dir,3,100*msyr))
  }
  
  ## NF-G trials changes in mixing matricies
  for(hypo in c(2,3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='G', optmix = 2),dir)
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='F', optmix = 3),dir)
      tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-G%s-%s.survey',dir,hypo,100*msyr))
      tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-F%s-%s.survey',dir,hypo,100*msyr))
    }
  }
  
  ## NF-S trials (selectivity before and after 2007 estimated)
  for(hypo in c(3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='S', selyr = 2007),dir)
      tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-S%s-%s.survey',dir,hypo,100*msyr))
      
    }
  }
  
  ## NF-Y trials (baseline with 8 year survey intervals)
  for(hypo in hypos){
    for(msyr in c(0.01,0.04)){
      if(hypo == 6 & msyr ==0.01){
        next
      } else {
        NafWrite(NafTrial(hypo = hypo,msyr = msyr,
                          trialtype='Y'),dir)
      }
      tmp <- file.copy('data/surveyn8.dat',sprintf('%s/NF-Y%s-%s.survey',dir,hypo,100*msyr))
    }
  }  
  
  ## NF-E trials exclude 1987/1989 abundance estimates in WI, EG & EI/F
  for(hypo in c(2,3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='E'),dir)
      tmp <- file.copy('data/survey-red.dat',sprintf('%s/NF-E%s-%s.survey',dir,hypo,100*msyr))
    }
  }
  
  ## NF-D trials put a lower upper bound on dispersion
  for(hypo in c(1,3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='D',dispb=0.2),dir)
      tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-D%s-%s.survey',dir,hypo,100*msyr))
    }
  }
  
  ## NF-J assumes a gzero of 0.8
  for(hypo in c(2,3)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = hypo, msyr = msyr,
                        trialtype='J', gzero = 0.8),dir)
      tmp <- file.copy('data/surveyn.dat',sprintf('%s/NF-J%s-%s.survey',dir,hypo,100*msyr))
      }
  }
  
  
    
  ## TODO: 
  ## - NF-X, NF-C (should we skip them?)
  ## - Missing data: NF-X NF-P
  ## - Future trials: NF-Q, Y
  
}

NafCond <- function(dir='trials',
                    nafpar='../variants/Naf-v0.par',
                    search.string = 'NF-[A-Z][0-9]-[0-9].dat',
                    ...){
  mclapply(list.files(dir, search.string),
           function(x){ 
             survey <- gsub('dat','survey',x)
             NafCall(run_dir=dir,copyna=x,nafpar=nafpar,survey=survey,...)
           },
           mc.cores = detectCores(logical = TRUE))
} 


NafVariants <- function(dir='trials',
                        var.dir ='variants',
                        var.search = 'Naf-v[0-9].par',
                        search.string = 'NF-[A-Z][0-9]-[0-9].dat',
                        man='../settings/manage.run', 
                        ...){
  file.copy('data/sur-v4.dat',sprintf('%s/SUR-V4.DAT',dir))  
  file.copy('data/CLC-N.PAR',sprintf('%s/clc-n.par',dir))
  
  expand.grid(ref = list.files(dir, pattern = search.string),
              variant = list.files(var.dir,pattern = var.search)) %>%
    split(rownames(.)) %>% 
  mclapply(function(x){
    file.copy(sprintf('%s/%s',dir,gsub('.dat','.all',x$ref)),
              sprintf('%s/%s',dir,gsub('.dat','.con',x$ref)))
    survey <- gsub('dat','survey',x$ref)
    NafCall(run_dir=dir,copyna=x$ref,
            nafpar=sprintf('../%s/%s',var.dir,x$variant),
            man=man,con=gsub('.dat','.con',x$ref),
            survey=survey,
            ...)
    },
    mc.cores = detectCores(logical = TRUE))
} 



NafDisp <- function(dir='disp-test'){
  dir.create(dir,showWarnings = FALSE)
  tmp <- 
    list.files('data') %>% 
    map(~file.copy(from=sprintf('data/%s',.),
                   to=sprintf('%s/%s',dir,.),
                   overwrite = TRUE))
  tmp <- file.copy('settings/manage-single.dat',sprintf('%s/manage.dat',dir))
  for(i in seq(0.05,0.5,by=0.02)){
    for(msyr in c(0.01,0.04)){
      NafWrite(NafTrial(hypo = 3, msyr = msyr,
                        trialtype=, disp1=0,
                        nf = round(100*i),
                        egwi_disp=i),dir)
    }
  }
}

