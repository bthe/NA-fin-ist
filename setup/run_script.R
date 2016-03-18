source('setup/trials.R')
source('setup/single-stock.R')
source('setup/read_variants.R')

ManSetup(dir='outn',db_name ='druna.db')
ManRun(dir='outn')
NafReadVariants(dir='druna_cond',search.string='NF-..-.-V[0-9].restest',db_name='druna_res.db')
ManResults(dir='outn',db_name = 'druna_res.db')


ManSetup(dir='outn-final',db_name ='druna_final.db')
ManRun(dir='outn-final')
NafReadVariants(dir='druna_final',search.string='NF-..-.-V[0-9].restest',db_name='druna_final_res.db')
ManResults(dir='outn-final',db_name = 'druna_final_res.db')

