templates = {}
## templates['CarboxyGroup']={'O10':{'sybylType':'O.co2','neighbours':['C3'],'alreadyvisited':False},
##                            'O20':{'sybylType':'O.co2','neighbours':['C3'],'alreadyvisited':False},
##                            'C3': {'sybylType':'C.2',  'neighbours':['O10','O20'],'alreadyvisited':False}}
templates['AceticAcid']={'O1' :{'sybylType':'O.co2','neighbours':['C'],'alreadyvisited':False},
                         'O2' :{'sybylType':'O.co2','neighbours':['C'],'alreadyvisited':False},
                         'C'  :{'sybylType':'C.2'  ,'neighbours':['O1','O2','C40'],'alreadyvisited':False},
                         'C40':{'sybylType':'C.3'  ,'neighbours':['C'],'alreadyvisited':False}}
templates['IDD594']={'O1':  {'sybylType':'O.co2','neighbours':['C'],'alreadyvisited':False},
                     'O2':  {'sybylType':'O.co2','neighbours':['C'],'alreadyvisited':False},
                     'C' :  {'sybylType':'C.2',  'neighbours':['O1','O2','C40'],'alreadyvisited':False},
                     'C40': {'sybylType':'C.3',  'neighbours':['C','O50'],'alreadyvisited':False},
                     'O50': {'sybylType':'O.3',  'neighbours':['C40','C60'],'alreadyvisited':False},
                     'C60': {'sybylType':'C.ar', 'neighbours':['O50'],'alreadyvisited':False}}
templates['acetylsalicylicacid']={'O1' :{'sybylType':'O.co2','neighbours':['C'],'alreadyvisited':False},
                                   'O2' :{'sybylType':'O.co2','neighbours':['C'],'alreadyvisited':False},
                                   'C'  :{'sybylType':'C.2'  ,'neighbours':['O1','O2','C40'],'alreadyvisited':False},
                                   'C20'  :{'sybylType':'C.3'  ,'neighbours':['C','C40','C30'],'alreadyvisited':False},
                                   'C30'  :{'sybylType':'C.3'  ,'neighbours':['C20','C50'],'alreadyvisited':False},
                                   'C50'  :{'sybylType':'C.3'  ,'neighbours':['C30','C80'],'alreadyvisited':False},
                                   'C80'  :{'sybylType':'C.3'  ,'neighbours':['C50','C60'],'alreadyvisited':False},
                                   'C60'  :{'sybylType':'C.3'  ,'neighbours':['C80','C40'],'alreadyvisited':False},
                                   'C40'  :{'sybylType':'C.3'  ,'neighbours':['C','C60'],'alreadyvisited':False}}
#templates['CRAP']={'O10':{'sybylType':'O.co2','neighbours':['C30'],'alreadyvisited':False},
#                   'O20':{'sybylType':'O.co2','neighbours':['C30'],'alreadyvisited':False},
#                   'C30':{'sybylType':'C.2'  ,'neighbours':['O10','O20','C40'],'alreadyvisited':False},
#                   'C40':{'sybylType':'C.3'  ,'neighbours':['C30','O50'],'alreadyvisited':False},
#                   'O50':{'sybylType':'O.3',  'neighbours':['C40'],'alreadyvisited':False}}
## templates['PropanoicAcid']={'O10':{'sybylType':'O.co2','neighbours':['C30'],'alreadyvisited':False},
##                             'O20':{'sybylType':'O.co2','neighbours':['C30'],'alreadyvisited':False},
##                             'C30':{'sybylType':'C.2'  ,'neighbours':['O10','O20','C40'],'alreadyvisited':False},
##                             'C40':{'sybylType':'C.3'  ,'neighbours':['C30','C50'],'alreadyvisited':False},
##                             'C50':{'sybylType':'C.3'  ,'neighbours':['C40'],'alreadyvisited':False}}
## templates['piperidine']={'N1':{ 'sybylType':'N.4','neighbours':['C1','C5'],'alreadyvisited':False},
##                          'C1':{ 'sybylType':'C.3','neighbours':['N1','C2'],'alreadyvisited':False},
##                          'C2':{ 'sybylType':'C.3','neighbours':['C1','C3'],'alreadyvisited':False},
##                          'C3':{ 'sybylType':'C.3','neighbours':['C2','C4'],'alreadyvisited':False},
##                          'C4':{ 'sybylType':'C.3','neighbours':['C3','C5'],'alreadyvisited':False},
##                          'C5':{ 'sybylType':'C.3','neighbours':['N1','C4'],'alreadyvisited':False}}
## templates['imidazole']={'N10':{ 'sybylType':'N.pl3','neighbours':['C30','C10'],'alreadyvisited':False},
##                         'C10':{ 'sybylType':'C.2',  'neighbours':['N10','N20'],'alreadyvisited':False},
##                         'N20':{ 'sybylType':'N.2',  'neighbours':['C10','C20'],'alreadyvisited':False},
##                         'C20':{ 'sybylType':'C.2',  'neighbours':['N20','C30'],'alreadyvisited':False},
##                         'C30':{ 'sybylType':'C.2',  'neighbours':['C20','N10'],'alreadyvisited':False}}
## templates['Ethine']={'C10':{ 'sybylType':'C.1','neighbours':['C20'],'alreadyvisited':False},
##                      'C20':{ 'sybylType':'C.1','neighbours':['C10'],'alreadyvisited':False}}
