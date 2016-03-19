source('setup/trials.R')
source('setup/single-stock.R')
source('setup/read_variants.R')
source('setup/read_results.R')


NafSetup('Abund')
NafCond(dir='RedAbund',survey='../data/survey-red.dat',
        search.string = 'NF-.[0-6]-..dat')
NafCond(dir='RedAbund',survey='../data/survey-red7.dat',
        search.string = 'NF-.[7|8]-..dat')
NafResults('RedAbund',db_name='RedAbund.db')
NafVariants(dir='RedAbund',var.search = 'Naf-v[0-2].par',
            search.string = 'NF-.[0-6]-..dat',
            survey='../data/survey-red.dat')
file.copy('settings/clc-n-rev.par',to = 'RedAbund/clc-n.par',overwrite = TRUE)
NafVariants(dir='RedAbund',var.search = 'Naf-v[3-7].par',
            search.string = 'NF-.[7-8]-..dat',
            survey='../data/survey-red7.dat')
ManSetup(dir='outn-Red',db_name ='RedAbund.db')
ManRun(dir='outn-Red')
NafReadVariants(dir='RedAbund',search.string='NF-..-.-V[0-9].restest',db_name='RedAbund_res.db')
ManResults(dir='outn-Red',db_name = 'RedAbund_res.db')

