import sys
sys.path.append('../src')
from readCHS import readCHS
from ios import *
import timeit, functools
import pandas as pd
import filecmp


def CHS():
  filePath = r'CHS_Constants.csv'
  outPath = r'35.TS.csv'

  datetimes = np.arange('2000-01-01', '2001-01-01', np.timedelta64(10, 'm'), dtype='datetime64')

  stations = readCHS(filePath)

  stations = stations[np.where(stations['id']==35)]
  cons = stations['constituents']['name'][0, :]

  cons = np.delete(cons, np.where(cons == 'Z0'))
  ios = IOS(cons,stations,debug=True)
  ts = ios.getTS(datetimes)
  ios.to_csv(outPath,type=0, index=None, header=True,date_format="%Y/%m/%d %H:%M:%S")
  
  
if __name__ == "__main__":
  CHS()
  # compare()