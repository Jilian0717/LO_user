""""
Dictionaries of defaults to be used for plotting.

"""

vlims_dict = {'salt': (27, 30),
        'temp': (7, 18),
        'dye_01': (0,1),
        'NO3': (0, 44),
        'NH4': (0, 20),
        'phytoplankton': (0,0.5),
        'zooplankton': (0, 4),
        'oxygen': (0, 50),
        'TIC': (2000, 2400),
        'alkalinity': (2000,2400),
        'PH': (7, 8.5),
        'ARAG': (.2, 2.2),
        'SdetritusN': (0,2),
        'LdetritusN': (0,0.1)}

vlims_dict_diff = {'salt': (-5, 5),
        'temp': (-5, 5),
        'dye_01': (-1,1),
        'NO3': (-10, 10),
        'NH4': (-10, 10),
        'phytoplankton': (-1,1),
        'zooplankton': (-1, 1),
        'oxygen': (-10, 10),
        'TIC': (-50, 50),
        'alkalinity': (-50,50),
        'SdetritusN': (-1,1),
        'LdetritusN': (-1,1)}

units_dict = {'salt': '$(g\ kg^{-1})$',
             'temp': ' $(^{\circ}C)$',
             'NO3': ' $(\mu mol\ L^{-1})$',
             'phytoplankton': ' $(mg\ Chl\ m^{-3})$',
             'zooplankton': ' $(\mu mol\ N\ L^{-1})$',
             'oxygen': ' $(mg\ L^{-1})$',
             'TIC': ' $(\mu mol\ L^{-1})$',
             'alkalinity': ' $(\mu\ equivalents\ L^{-1})$',
             'PH': '',
             'ARAG': '',
             'detritus': ' $(\mu mol\ L^{-1})$',
             'Ldetritus': ' $(\mu mol\ L^{-1})$',
             'w': ' $(m\ s^{-1})$',
             'u': ' $(m\ s^{-1})$',
             'v': ' $(m\ s^{-1})$',
             'ubar': ' $(m\ s^{-1})$',
             'vbar': ' $(m\ s^{-1})$',
             'zeta': ' (m)'}
