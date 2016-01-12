library(dplyr)
library(infuser)
library(purrr)

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
trial <- function(hypo=1,msyr=0.01,trialtype='B',
                  template_file='setup/trial_template.txt'){
  
  ## default values for variables
  nstk <- 6
  nsuba <- 7
  ngamma <- 6
  ndelta <- 2
  initdelta <- '0.15   0.02'
  initgamma <- 1.0
  
  ## deviations  
  ## Fewer stocks in hypotheses 6 to 8
  if(hypo %in% c(6,7,8)){
    nstk <- 5
  }
  ## WI and EG are merged
  if(hypo %in% c(7,8)){
    nsuba <- 6
  }
  
  ## mixing used instead of dispersion in hypotheses 4 and 8
  if(hypo %in% c(4,8)){
    ndelta <- 0
    initdelta <- ''
    initgamma <- 0.2
  }
  
  
  ## mixing matricies
  if(hypo==1){ 
    mixmat <- 
"1    1          1   1   0   0    0   0   0    1  4  0  0  0  0  0
2    1          0   0   1   0    0   0   0    0  0  0  0  0  0  0
3    1          0   0   0   1    0   0   0    0  0  0  0  0  0  0
4    1          0   0   0   0    1   0   0    0  0  0  0  0  0  0
5    1          0   0   0   0    0   1   0    0  0  0  0  0  0  0
6    1          0   0   0   0    0   0   1    0  0  0  0  0  0  0"
    
  } else if(hypo==2) {
    mixmat <-
"1    1         .88 .88 .1  .02   0   0   0    1  4  0  0  0  0  0 
2    1          0   0   1   0    0   0   0    0  0  0  0  0  0  0
3    1          0   0   0   1    0   0   0    0  0  0  0  0  0  0
4    1          0   0   0   0    1   0   0    0  0  0  0  0  0  0
5    1          0   0   0  .02  .1  .88  0    0  0  0  0  0  0  0
6    1          0   0   0   0    0   0   1    0  0  0  0  0  0  0"
  } else if(hypo==3) {
    mixmat <- 
"1    1          1   1   0   0   0   0   0     1  4  0  0  0  0  0
2    1         .1  .1  .9   0   0   0   0     1  4  0  0  0  0  0
3    1          0   0   0   1   0   0   0     0  0  0  0  0  0  0
4    1          0   0   0   0  .9  .1   0     0  0  0  0  0  0  0
5    1          0   0   0   0   0   1   0     0  0  0  0  0  0  0
6    1          0   0   0   0   0   0   1     0  0  0  0  0  0  0"
  } else if(hypo == 4){
    mixmat <- 
"1    1          1   1   0   0   0   0   0     1  4  0  0  0  0  0
2    1          1   1   1   1   0   0   0     2  6  5  3  0  0  0
3    1          0   0   1   1   1   0   0     0  0  3  5  3  0  0
4    1          0   0   0   1   1   1   0     0  0  0  3  5  3  0
5    1          0   0   0   0   0   1   0     0  0  0  0  0  0  0
6    1          0   0   0   0   0   0   1     0  0  0  0  0  0  0"
  } else if(hypo == 5){
    mixmat <- 
"1    1          1   1   0   0   0   0   0     1  4  0  0  0  0  0
2    1          0   0   1   0   0   0   0     0  0  0  0  0  0  0
3    1          0   0   0   1   0   0   0     0  0  0  0  0  0  0
4    1          0   0   0   0   1   0   0     0  0  0  0  0  0  0
5    1          0   0   0   0   0   1   0     0  0  0  0  0  0  0
6    1          0   0   0   0  .02 .1  .88    0  0  0  0  0  0  0"
  } else if(hypo == 6){
    mixmat <- 
"1    1         .88 .88 .1  .02   0   0   0    1  4  0  0  0  0  0
2    1          0   0   1   0    0   0   0    0  0  0  0  0  0  0
3    1          0   0   0   1    0   0   0    0  0  0  0  0  0  0
4    1          0   0   0   0    1   1   0    0  0  0  0  3  2  0
5    1          0   0   0   0    0   0   1    0  0  0  0  0  0  0"
  } else if(hypo == 7){
    mixmat <- 
"1    1          1   1   0   0    0   0        1  4  0  0  0  0
2    1          1   1   1   0    0   0        5  6  3  0  0  0
3    1          0   0   0  .9   .1   0        0  0  0  0  0  0
4    1          0   0   0   0    1   0        0  0  0  0  0  0
5    1          0   0   0   0    0   1        0  0  0  0  0  0"
  } else if(hypo == 8){
    mixmat <- 
"1    1          1   1   0   0    0   0        1  4  0  0  0  0 
2    1          1   1   1   1    0   0        2  6  5  3  0  0 
3    1          0   0   1   1    1   0        0  0  3  5  3  0 
4    1          0   0   0   0    1   0        0  0  0  0  0  0 
5    1          0   0   0   0    0   1        0  0  0  0  0  0"
  }
  
  ## trial code generated
  ref <- 'NF-{{trialtype}}{{hypo}}-{{msyr}} H{{hypo}}-{{msyr}}' %>%
    infuse(hypo=hypo,msyr=100*msyr,trialtype=trialtype)
  
  ## assign msyr rates to nstk stocks
  msyr.txt <- 
    sprintf('MSYR OF STOCK 1 (MSYR%1$s(%1$s))                   %2$s',
            1:nstk,msyr) %>%
    paste(collapse='\n')
  
  ## change the density dependance from 1+ for 1% to mature for 4%
  if(msyr==0.01){
    optf <- 1
  } else {
    optf <- 0
  }
  
  ## produce the trial setup
  txt <- infuse(template_file,
         ref=ref,
         msyrtxt=msyr.txt,
         optf=optf,nstk=nstk,
         mixmat=mixmat,
         ndelta=ndelta,
         initdelta = initdelta,
         initgamma = initgamma,
         disp=ifelse(ndelta==0,0,1))
  attr(txt,'ref') <- 'NF-{{trialtype}}{{hypo}}-{{msyr}}' %>%
    infuse(hypo=hypo,msyr=100*msyr,trialtype=trialtype)
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
Nafwrite <- function(obj,dir='.'){
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
callNaf <- function(...){
  tmp <- list(...)
  if('run_dir' %in% names(tmp)){
    old.dir <- getwd()
    setwd(tmp$run_dir)
  }
  run.string <- 'naf-ist -main {{copyna|copyna.dat}} -par {{nafpar|nafpar.par}}' %>%
    infuse(tmp)
  res <- system(run.string,ignore.stdout = TRUE,ignore.stderr = FALSE,
                intern = FALSE)
  if('run_dir' %in% names(tmp)){
    setwd(old.dir)
  }
}


NafSetup <- function(dir = 'trials'){
  dir.create(dir,showWarnings = FALSE)
  tmp <- 
    list.files('data') %>% map(~file.copy(from=sprintf('data/%s',.),
                                        to=sprintf('%s/%s',dir,.),
                                        overwrite = TRUE))
  ## NF-B1-1%
  Nafwrite(trial(hypo = 1,msyr = 0.01),dir)
  ## NF-B1-4%
  Nafwrite(trial(hypo = 1,msyr = 0.04),dir)
  ## NF-B2-1%
  Nafwrite(trial(hypo = 2,msyr = 0.01),dir)
  ## NF-B2-4%
  Nafwrite(trial(hypo = 2,msyr = 0.04),dir)
  ## NF-B3-1%
  Nafwrite(trial(hypo = 3,msyr = 0.01),dir)
  ## NF-B3-4%
  Nafwrite(trial(hypo = 3,msyr = 0.04),dir)
  ## NF-B4-1%
  Nafwrite(trial(hypo = 4,msyr = 0.01),dir)
  ## NF-B4-4%
  Nafwrite(trial(hypo = 4,msyr = 0.04),dir)
  ## NF-B5-1%
  Nafwrite(trial(hypo = 5,msyr = 0.01),dir)
  ## NF-B5-4%
  Nafwrite(trial(hypo = 5,msyr = 0.04),dir)
  ## NF-B6-1%
  Nafwrite(trial(hypo = 6,msyr = 0.01),dir)
  ## NF-B6-4%
  Nafwrite(trial(hypo = 6,msyr = 0.04),dir)
  ## NF-B7-1%
  Nafwrite(trial(hypo = 7,msyr = 0.01),dir)
  ## NF-B7-4%
  Nafwrite(trial(hypo = 7,msyr = 0.04),dir)
  ## NF-B8-1%
  Nafwrite(trial(hypo = 8,msyr = 0.01),dir)
  ## NF-B8-4%
  Nafwrite(trial(hypo = 8,msyr = 0.04),dir)
}

