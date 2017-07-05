from readCHS import *
from ios import *
import timeit, functools
# ----------------------------------------------------------------------------------------------------------------------
# Test Functions
def testCHS1():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2003-01-01', np.timedelta64(3, 'h'), dtype='datetime64')
    stations = getCHS(csvPath)
    cons = stations['constituents']['name'][0, :]
    stations = stations[:1000]
    # stations = np.concatenate((stations, stations, stations, stations, stations, stations, stations, stations, stations, stations))
    cons = np.array(['Z0', 'M2', 'S2', 'M4', 'K1'])

    t = timeit.Timer(functools.partial(ios.Run, datetimes=datetimes, stations=stations,cons=cons))
    print("Testing CHS(# of station: {0}, # of steps: {1}) - {2}".format(len(stations), len(datetimes), t.timeit(1)))

def testCHS2():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2003-01-01', np.timedelta64(3, 'h'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:1]
    cons = stations['constituents']['name'][0, :]

    t = timeit.Timer(functools.partial(ios.Run, datetimes=datetimes, stations=stations,cons=cons))
    print("Testing CHS(# of station: {0}, # of steps: {1}) - {2}".format(len(stations),len(datetimes),t.timeit(1)))


def testCHS1Parallel():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2003-01-01', np.timedelta64(3, 'h'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:1000]
    # stations = np.concatenate((stations, stations, stations, stations, stations, stations, stations, stations, stations, stations))
    # cons = stations['constituents']['name'][0, 1]
    cons = np.array(['Z0', 'M2', 'S2', 'M4', 'K1'])

    t = timeit.Timer(functools.partial(ios.RunP,nprocessor=10, datetimes=datetimes, stations=stations,cons=cons))

    print("Testing CHS Parallel(# of station: {0}, # of steps: {1}) - {2}".format(len(stations), len(datetimes), t.timeit(1)))


# def testLatitude():
#     csvPath = r'CHS_Constants.csv'
#
#     ios = IOS()
#     datetimes = np.arange('2000-01-01', '2000-01-02', np.timedelta64(15, 'm'), dtype='datetime64')
#     stations = getCHS(csvPath)
#     stations = stations[:4]
#
#     stations[1] = stations[0]
#     stations[2] = stations[0]
#     stations[3] = stations[0]
#     stations['xy'][1, 1] = stations['xy'][0, 1] + 0.1
#     stations['xy'][2, 1] = stations['xy'][0, 1] + 1.0
#     stations['xy'][3, 1] = stations['xy'][0, 1] + 10.0
#     t = timeit.Timer(functools.partial(ios.generateTS, datetimes=datetimes, stations=stations))
#     print("Testing CHS (latitude) - " + t.timeit(1))
#




if __name__ == "__main__":
    testCHS1()
    # testCHS2()
    #testCHS1Parallel()
# ----------------------------------------------------------------------------------------------------------------------
