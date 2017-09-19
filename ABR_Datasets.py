#current 10 June 2017
# This table should eventually be converted to an excel spreadsheet organized
# with a worksheet for each data set, and within each worksheet, each subject's information
# placed in a single row. Make the table pandas-readable.

# on the other hand, for now, it is useful as is.

ABR_Datasets = { 'NrCAMKO': {'dir': "Tessa/Coate-NrCAM-ABRs/KO", 
                    'invert': True, 
                    'clickselect: ': [['0849'], None, None, None, None, None, None, None],
                    'toneselect': [['0901'], ['1446', '1505'], None, None, None, None, None, None],
                    'term': '\r', 'minlat': 2.2},
                 'NrCAMWT': {'dir': "Tessa/Coate-NrCAM-ABRs/WT", 
                                        'invert': True, 
                                        'term': '\r', 'minlat': 2.2},
                 'TessaCBA': {'dir': 'Tessa/CBA', 'term': '\r', 'minlat': 2.2, 'invert': True},
                 'TessaCBANE': {'dir': 'Tessa/CBA_NoiseExposed', 'term': '\r', 'minlat': 2.2, 'invert': True},
                 'CNTNAP2X': {'dir': 'Tessa/CNTNAP2', 'term': '\r', 'minlat': 2.2, 'invert': True},
                 'CNTNAP2Het': {'dir': 'Tessa/CNTNAP2_Het', 'term': '\r', 'minlat': 2.2, 'invert': True},
                 'CNTNAP2HG': {'dir': 'Tessa/CNTNAP2_Het_GP4.3', 'term': '\r', 'minlat': 2.2, 'invert': True},
                 'CNTNAP2KO': {'dir': 'Tessa/CNTNAP2_KO',
                     'term': '\r', 'minlat': 2.2, 'invert': True},
                 'CNTNAP2WT': {'dir': 'Tessa/CNTNAP2_WT', 'term': '\r', 'minlat': 2.2, 'invert': True},
                 'GP43Norm': {'dir': 'Tessa/GP4.3-Thy1-Normal', 'term': '\r', 'minlat': 2.2, 'invert': True},
                 'Yong': {'dir': "Yong\'s ABRs", 'invert': True,
                     'term':'\n', 'minlat': 0.6},
                 'Jun': {'dir': "JUN\'s ABRs", 'term': '\r', 'minlat': 2.2, 'invert': False},
                 'Ruili': {'dir': "Ruilis ABRs", 'invert': True, 'nameselect': 'CBA',
                     'term':'\n', 'minlat': 0.6},
                 'RuiliCBAP40': {'dir': 'RuiliABRData_2010-2015/CBA-P21-P40',
                     'invert': True, 'nameselect': 'CBA',
                     'term':'\n', 'minlat': 0.6},
                 'RuiliCBAP20': {'dir': 'RuiliABRData_2010-2015/CBA-P10-P20',
                     'invert': True, 'nameselect': 'CBA',
                     'term':'\n', 'minlat': 0.6},
                 'Eveleen': {'dir': "Eveleen\'s ABRs", 'term': '\r', 'minlat': 2.2, 'invert': False},
                 'Amber': {'dir': "Amber_ABR_data", 'term': '\r', 'minlat': 2.2, 'invert': True}
                 }
