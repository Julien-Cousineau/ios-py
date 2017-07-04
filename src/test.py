from readCHS import *

# ----------------------------------------------------------------------------------------------------------------------
# Test Functions
def testCHS1():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2000-02-01', np.timedelta64(15, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)

    t = timeit.Timer(functools.partial(ios.generateTS, datetimes=datetimes, stations=stations))
    print("Testing CHS(# of station: %1, # of steps: %2) - %3".format(len(stations), len(datetimes), t.timeit(1)))

def testCHS2():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2001-01-01', np.timedelta64(15, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:1]

    t = timeit.Timer(functools.partial(ios.generateTS, datetimes=datetimes, stations=stations))
    print("Testing CHS(# of station: %1, # of steps: %2) - %3".format(len(stations),len(datetimes),t.timeit(1)))


def testCHSParallel():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2001-01-01', np.timedelta64(15, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:1000]

    t = functools.partial(ios.pGenerateTS, datetimes=datetimes, stations=stations)
    p = Pool(10)
    print("Testing CHS (parallel) - " + timeit.timeit(lambda: p.map(t, xrange(1000)), number=1))

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
    t = timeit.Timer(functools.partial(ios.generateTS, datetimes=datetimes, stations=stations))
    print("Testing CHS (latitude) - " + t.timeit(1))





if __name__ == "__main__":
    testCHS1()
    testCHS2()
# ----------------------------------------------------------------------------------------------------------------------
