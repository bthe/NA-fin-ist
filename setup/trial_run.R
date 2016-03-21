source('setup/trials.R')
source('setup/single-stock.R')
source('setup/read_variants.R')
source('setup/read_results.R')


NafSetup('FullTrials')
NafCond(dir='FullTrials',nafist='../naf-ist')
NafResults('FullTrials',db_name='FullTrials.db')
NafVariants(dir='FullTrials',var.search = 'Naf-v[0-2].par')
file.copy('settings/clc-n-rev.par',to = 'FullTrials/clc-n.par',overwrite = TRUE)
NafVariants(dir='FullTrials',var.search = 'Naf-v[3-7].par')
ManSetup(dir='outn-full',db_name ='FullTrials.db')
ManRun(dir='outn-full')
NafReadVariants(dir='FullTrials',search.string='NF-..-.-V[0-9].restest',db_name='FullTrials_res.db')
ManResults(dir='outn-Trials',db_name = 'FullTrials_res.db')

