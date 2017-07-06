from readCHS import *
from ios import *
import timeit, functools
from datetime import date, datetime, timedelta
# ----------------------------------------------------------------------------------------------------------------------
# Test Functions
def testCHS1():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2000-02-01', np.timedelta64(10, 'm'), dtype='datetime64')
    stations = getCHS(csvPath)
    cons = stations['constituents']['name'][0, :]
    # stations = stations[:1]
    stations=stations[np.where(stations['id'] == 60)]
    # stations = np.concatenate((stations, stations, stations, stations, stations, stations, stations, stations, stations, stations))
    # cons = np.array(['Z0', 'M2'])

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

def testCHSExtract():
    csvPath = r'CHS_Constants.csv'

    ios = IOS()
    datetimes = np.arange('2000-01-01', '2000-02-01', np.timedelta64(3, 'h'), dtype='datetime64')
    stations = getCHS(csvPath)
    stations = stations[:15]
    cons = np.array(['Z0', 'M2', 'S2', 'M4', 'K1'])

    t = timeit.Timer(functools.partial(ios.extractConstituents, datetimes=datetimes,ts=[], stations=stations,cons=cons))
    print("Testing CHS(# of station: {0}, # of steps: {1}) - {2}".format(len(stations), len(datetimes), t.timeit(1)))


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

def testCHSExtractInput():
    csvPath = r'CHS_Constants.csv'
    inputTS = r'timeseries.csv'
    ios = IOS()

    vars = [dict(name="Date", type="date"),
            dict(name="WL", type="float"),
            dict(name="U", type="float"),
            dict(name="V", type="float")]
    data = []
    with open(inputTS, 'rb') as csvfile:
        reader = csv.DictReader(csvfile)
        header = reader.fieldnames

        for row in reader:
            tvars = []
            tdata = []
            for var in vars:
                if var['name'] in header:
                    tvars.append(var['name'])
                    if var['type'] == 'date':
                        tdata.append(datetime.strptime(row[var['name']], "%d/%m/%Y %H:%M"))
                    if var['type'] == 'float':
                        tdata.append(float(row[var['name']]))
            t = dict(zip(tvars, tdata))
            data.append(t)


    start = datetime.strptime('2000-01-0120:00:00', "%Y-%m-%d%H:%M:%S")
    end = datetime.strptime('2000-02-0100:00:00', "%Y-%m-%d%H:%M:%S")

    fdata = [row for row in data if start <= row['Date'] <= end]

    datetimes = []
    WLs = []
    Us = []
    Vs = []
    for row in fdata:
        datetimes.append(np.datetime64(row['Date']))
        WLs.append(row['WL'])
        Us.append(row['U'])
        Vs.append(row['V'])

    datetimes = np.asarray(datetimes)
    WLs=np.asarray(WLs)
    Us=np.asarray(Us)
    Vs=np.asarray(Vs)

    stations = getCHS(csvPath)
    # np.where((stations['xy'] > 45.15) & (stations['xy'] < 45.25))[0]
    stations = stations[18:19]
    cons = stations['constituents']['name'][0, :]
    cons = np.array(['Z0', 'M2'])



    t = timeit.Timer(functools.partial(ios.extractConstituents, datetimes=datetimes, ts=np.asarray([[WLs,Us,Vs]]), stations=stations, cons=cons))
    print("Testing CHS(# of station: {0}, # of steps: {1}) - {2}".format(len(stations), len(datetimes), t.timeit(1)))


if __name__ == "__main__":
    # testCHS1()
    # testCHS2()
    #testCHS1Parallel()
    # testCHSExtract()
    testCHSExtractInput()
# ----------------------------------------------------------------------------------------------------------------------
