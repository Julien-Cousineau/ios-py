import csv
import numpy as np
from ios import stationType

def readCHS(filePath):
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

