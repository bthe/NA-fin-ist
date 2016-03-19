source('setup/trials.R')
source('setup/single-stock.R')
source('setup/read_variants.R')
source('setup/read_results.R')


NafSetup('FullAbund')
NafCond(dir='FullAbund')
NafResults('FullAbund',db_name='FullAbund.db')
NafVariants(dir='FullAbund',var.search = 'Naf-v[0-2].par')
file.copy('settings/clc-n-rev.par',to = 'FullAbund/clc-n.par',overwrite = TRUE)
NafVariants(dir='FullAbund',var.search = 'Naf-v[3-7].par')
ManSetup(dir='outn-full',db_name ='FullAbund.db')
ManRun(dir='outn-full')
NafReadVariants(dir='FullAbund',search.string='NF-..-.-V[0-9].restest',db_name='FullAbund_res.db')
ManResults(dir='outn-full',db_name = 'FullAbund_res.db')

