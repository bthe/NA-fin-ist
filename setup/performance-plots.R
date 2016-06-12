dat <-
  trial_stat$pop.res %>% 
  ungroup() %>% 
  mutate(variant=as.numeric(gsub('V','',variant))) %>% 
  filter(pop_id%in%c('C1','C2','C3')) %>%
  select(ref,variant,pop_id,res_fin=final_dpl_lower,
         res_min=pmin_lower) %>% 
  gather(type,dpl,c(res_min,res_fin)) %>% 
  mutate(type=ifelse(type=='res_min','P-min','P-fin')) %>% 
  left_join(trial_stat$pop.res %>%
              ungroup() %>% 
              filter(variant == 'V0') %>% 
              select(ref,pop_id,dplfin_60,dplfin_72,dplmin_60,dplmin_72) %>% 
              gather(type,boundary,-c(ref,pop_id)) %>% 
              separate(type,c('type','stat')) %>% 
              mutate(type = ifelse(type=='dplfin','P-fin','P-min'),
                     stat = ifelse(stat == '60','B','A')) %>% 
              spread(stat,boundary)) 

## baseline trials
print(
dat %>% 
  filter(grepl('NF-B.-1',ref)) %>% 
  ggplot(aes(variant,dpl)) + 
  geom_rect(aes(ymax=B),xmin=-Inf,xmax=Inf,ymin=-Inf,
            fill = 'gray70') +
  geom_hline(aes(yintercept=A),lty = 2) +
  geom_hline(aes(yintercept=B),lty = 2) +
  geom_point() + facet_grid(ref~type+pop_id) + 
  theme_bw() + 
  expand_limits(y = 0) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=7))+
  ggtitle('Depletion indicators for the Baseline trials') + 
  ylab('Depleption indicator') + xlab('Management variant') 
  )

print(
  dat %>% 
    filter(grepl('NF-Y.-1',ref)) %>% 
    ggplot(aes(variant,dpl)) + 
    geom_rect(aes(ymax=B),xmin=-Inf,xmax=Inf,ymin=-Inf,
              fill = 'gray70') +
    geom_hline(aes(yintercept=A),lty = 2) +
    geom_hline(aes(yintercept=B),lty = 2) +
    geom_point() + facet_grid(ref~type+pop_id) + 
    theme_bw() + 
    expand_limits(y = 0) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size=7))+
    ggtitle('Depletion indicators for the Baseline trials w. 8yr survey interval') + 
    ylab('Depleption indicator') + xlab('Management variant') 
)

print(
  dat %>% 
    filter(grepl('NF-[^(B|E-Y)].-1',ref)) %>% 
    ggplot(aes(variant,dpl)) + 
    geom_rect(aes(ymax=B),xmin=-Inf,xmax=Inf,ymin=-Inf,
              fill = 'gray70') +
    geom_hline(aes(yintercept=A),lty = 2) +
    geom_hline(aes(yintercept=B),lty = 2) +
    geom_point() + facet_grid(ref~type+pop_id) + 
    theme_bw() + 
    expand_limits(y = 0) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size=7))+
    ggtitle('Depletion indicators for trials A and D ') + 
    ylab('Depleption indicator') + xlab('Management variant') 
)

print(
  dat %>% 
    filter(grepl('NF-[^(A-D|G-Y)].-1',ref)) %>% 
    ggplot(aes(variant,dpl)) + 
    geom_rect(aes(ymax=B),xmin=-Inf,xmax=Inf,ymin=-Inf,
              fill = 'gray70') +
    geom_hline(aes(yintercept=A),lty = 2) +
    geom_hline(aes(yintercept=B),lty = 2) +
    geom_point() + facet_grid(ref~type+pop_id) + 
    theme_bw() + 
    expand_limits(y = 0) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size=7))+
    ggtitle('Depletion indicators for trials E and F') + 
    ylab('Depleption indicator') + xlab('Management variant') 
)


print(
  dat %>% 
    filter(grepl('NF-[^(A-F|I-Y)].-1',ref)) %>% 
    ggplot(aes(variant,dpl)) + 
    geom_rect(aes(ymax=B),xmin=-Inf,xmax=Inf,ymin=-Inf,
              fill = 'gray70') +
    geom_hline(aes(yintercept=A),lty = 2) +
    geom_hline(aes(yintercept=B),lty = 2) +
    geom_point() + facet_grid(ref~type+pop_id) + 
    theme_bw() + 
    expand_limits(y = 0) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size=7))+
    ggtitle('Depletion indicators for trials G and H') + 
    ylab('Depleption indicator') + xlab('Management variant') 
)


print(
  dat %>% 
    filter(grepl('NF-[^(A-H|Y)].-1',ref)) %>% 
    ggplot(aes(variant,dpl)) + 
    geom_rect(aes(ymax=B),xmin=-Inf,xmax=Inf,ymin=-Inf,
              fill = 'gray70') +
    geom_hline(aes(yintercept=A),lty = 2) +
    geom_hline(aes(yintercept=B),lty = 2) +
    geom_point() + facet_grid(ref~type+pop_id) + 
    theme_bw() + 
    expand_limits(y = 0) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size=7))+
    ggtitle('Depletion indicators for trials J, Q, S and U') + 
    ylab('Depleption indicator') + xlab('Management variant') 
)
