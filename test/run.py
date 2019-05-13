import sys
sys.path.append('../src')
from readCHS import readCHS
from ios import *
import timeit, functools
import pandas as pd
import filecmp


def NOAA():
  filePath = r'noaas.csv'
  outPath = r'8410140.Julien.txt'

  datetimes = np.arange('2001-06-21', '2001-09-22', np.timedelta64(10, 'm'), dtype='datetime64')

  stations = readCHS(filePath)

  stations = stations[np.where(stations['id']==8410140)]
  cons = stations['constituents']['name'][0, :]

  cons = np.delete(cons, np.where(cons == 'Z0'))
  ios = IOS(cons,stations,debug=True)
  ts = ios.getTS(datetimes)
  ios.to_csv(outPath,type=0, index=None, header=True,date_format="%Y/%m/%d %H:%M:%S")
  
def compare():
  print filecmp.cmp('8410140.Julien.template.txt','8410140.Julien.txt')
  
  
if __name__ == "__main__":
  NOAA()
  compare()