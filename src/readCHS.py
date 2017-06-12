
import numpy as np
import timeit, functools,csv

from multiprocessing import Pool
from ios import *

# ------------------------------------------------------------------------------
# Type
constituentType = np.dtype([('name', '|S8'),
                            ('eta', 'f8', (1, 2)),
                            ('u', 'f8', (1, 2)),
                            ('v', 'f8', (1, 2)),
                            ('r', 'f8'),
                            ('zeta', 'f8'),
                            ('konans', 'i4')])

stationType = np.dtype([('name', '|S32'),
                        ('id', 'i4'),
                        ('xy', '2f8',),
                        ('constituents', constituentType, 200)])
# ------------------------------------------------------------------------------
# Functions
def getCHS(filePath):
    # type: (object, object) -> object
    with open(filePath, 'rb') as file:
        reader = csv.DictReader(file)
        header = reader.fieldnames
        ncount=0
        for row in reader:
            ncount +=1
        file.seek(0)
        reader = csv.DictReader(file)

        stations = np.zeros(ncount, dtype=stationType)
        istation = 0
        for row in reader:
            id = int(float(row["ID"]))
            longitude = float(row["Longitude"])
            latitude = float(row["Latitude"])
            station_name = row["Station Name"]
            icon = 0
            for i in range(4, len(row), 3):
                constituent_name = header[i].split('_')[0]
                strAmp = constituent_name + '_Ampm'
                strLocalPhase = constituent_name + '_LocalPh'
                strUTMPhase = constituent_name + '_UTCPhas'

                amplitude = float(row[strAmp[0:10]])
                localPhase = float(row[strLocalPhase[0:10]])
                utmPhase = float(row[strUTMPhase[0:10]])

                stations['constituents']['name'][istation, icon] = constituent_name
                stations['constituents']['eta'][istation, icon] = (amplitude, utmPhase)
                icon += 1
            stations['name'][istation] = station_name
            stations['id'][istation] = id
            stations['xy'][istation] = (longitude, latitude)
            istation += 1
        return stations

# ------------------------------------------------------------------------------
# Test Functions
def testCHS():
    csvPath = r'CHS_Constants.csv'
    
    ios = IOS()
    datetimes = np.arange('2000-01-01', '2001-01-01', np.timedelta64(15, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)

    t = timeit.Timer(functools.partial(ios.generateTS,datetimes=datetimes, stations=stations))
    print t.timeit(1)

def testCHSParallel():
    csvPath = r'CHS_Constants.csv'
    
    ios = IOS()
    datetimes = np.arange('2000-01-01', '2001-01-01', np.timedelta64(15, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:1000]

    t = functools.partial(ios.pGenerateTS, datetimes=datetimes, stations=stations)
    p = Pool(15)
    print timeit.timeit(lambda: p.map(t, xrange(1000)), number=1)

def testLatitude():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2001-01-01', np.timedelta64(15, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:4]

    stations[1] = stations[0]
    stations[2] = stations[0]
    stations[3] = stations[0]
    stations['xy'][1, 1] = stations['xy'][0, 1] + 0.1
    stations['xy'][2, 1] = stations['xy'][0, 1] + 1.0
    stations['xy'][3, 1] = stations['xy'][0, 1] + 10.0
    t = timeit.Timer(functools.partial(ios.generateTS,datetimes=datetimes, stations=stations))
    print t.timeit(1)
    
def testSomething():
    csvPath = r'CHS_Constants.csv'
    
    ios = IOS()
    datetimes = np.arange('2000-01-01', '2001-01-01', np.timedelta64(15, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:4]
    
    t = timeit.Timer(functools.partial(ios.generateTS,datetimes=datetimes, stations=stations))
    print t.timeit(1)
    
