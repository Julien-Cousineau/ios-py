import csv
import numpy as np
from dtypes import stationType
import pandas as pd

def to_csv(filePath,stations,type='eta'):
    
    tmpdata = np.column_stack((stations['xy'][:,0],stations['xy'][:,1],stations['name'],stations['id']))
    tmpdata = pd.DataFrame(data=tmpdata,columns=['Longitude','Latitude','Station Name','ID'])
    
    for con in  stations[0]['constituents']:
        if len(con['name'])==0:continue
        if type=='eta':
            tmpdata[con['name']+"_Ampm"]=con['eta'][0]
            tmpdata[con['name']+"_LocalPh"]=-1.0
            tmpdata[con['name']+"_UTCPhas"]=con['eta'][1]
        else:
            tmpdata[con['name']+"_u"]=con['u'][0]
            tmpdata[con['name']+"_v"]=con['v'][0]
            tmpdata[con['name']+"_uUTC"]=con['u'][1]
            tmpdata[con['name']+"_vUTC"]=con['v'][1]
    tmpdata.to_csv(filePath,index=None, header=True)

def read_csv(filePath):
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
                strAmp = strAmp[0:12]
                # strLocalPhase = constituent_name + '_LocalPh'
                strUTMPhase = constituent_name + '_UTCPhas'
                strUTMPhase = strUTMPhase[:12]

                amplitude = float(row[strAmp])
                # localPhase = float(row[strLocalPhase[0:10]])
                utmPhase = float(row[strUTMPhase])

                stations['constituents']['name'][istation, icon] = constituent_name
                stations['constituents']['eta'][istation, icon] = (amplitude, utmPhase)
                icon += 1
            stations['name'][istation] = station_name
            stations['id'][istation] = id
            stations['xy'][istation] = (longitude, latitude)
            istation += 1
        return stations

def read_npy(filepath):
    fp = np.memmap(filepath, dtype='int32', mode='r',shape=1)
    return np.memmap(filepath, dtype=stationType, mode='r+', shape=fp[0],offset=4)
  
def to_npy(filepath,stations):
    nstations = len(stations)
    fp = np.memmap(filepath, dtype='int32', mode='w+', shape=1)
    fp[0]=nstations
    fp = np.memmap(filepath, dtype=stationType, mode='r+', shape=nstations,offset=4)
    fp[:] = stations
    return fp

