import numpy as np

constituentType = np.dtype([('name', '|S8'),
                            ('eta', '2f8'),
                            ('u', '2f8'),
                            ('v', '2f8'),
                            ('r', 'f8'),
                            ('zeta', 'f8'),
                            ('konans', 'i4')])

stationType = np.dtype([('name', '|S32'),
                        ('id', 'i4'),
                        ('xy', '2f8'),
                        ('proj', '|S16'),
                        ('constituents', constituentType, 200)])
                        
def createStation(name,id,xy,proj,constituents=None):
    station = np.zeros(1, dtype=stationType)
    station['name'][0] = name
    station['id'][0]   = id
    station['xy'][0]   = xy
    station['proj'][0] =proj
    if constituents==None:
        station['constituents']['name'][0][0:73]=[
            'Z0',
 'O1',
 'K1',
 'M2',
 'S2',
 'M4',
 'MS4',
 'N2',
 'K2',
 'Q1',
 'J1',
 'MU2',
 'P1',
 'MM',
 'OO1',
 'L2',
 'NO1',
 'MN4',
 'M8',
 '2Q1',
 '2N2',
 'NU2',
 'M3',
 'MO3',
 'MK3',
 'SK3',
 'SN4',
 'S4',
 '2MN6',
 'M6',
 '2MS6',
 '2SM6',
 'MSF',
 'ALP1',
 'UPS1',
 'EPS2',
 'ETA2',
 '2MK5',
 '2SK5',
 '3MK7',
 'M10',
 'SIG1',
 'OQ2',
 'MSN2',
 '2SM2',
 'SK4',
 'MF',
 'SSA',
 'PHI1',
 'SO1',
 'MKS2',
 'SO3',
 'MK4',
 '2MK6',
 'MSK6',
 'TAU1',
 'BET1',
 'SA',
 'RHO1',
 'PI1',
 'S1',
 'LDA2',
 'T2',
 'OP2',
 'MSN6',
 'CHI1',
 'PSI1',
 'THE1',
 'R2',
 'MSM',
 'H1',
 'H2',
 'GAM2']
    return station
        