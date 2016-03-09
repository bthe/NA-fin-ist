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
    mutate(dpl=1-x2014/x1864)%>% 
    select(ref,trial,pop_id,dpl) 

  tmp_func <- function(x){
    write(x,file=sprintf('%s/%s.dpl',dir,x))
    pop %>% 
      filter(ref==x) %>% 
      spread(pop_id,dpl) %>% 
      select(-c(ref,trial)) %>%
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
    map(~tmp_fuse(.))
}


  
  