import sys
sys.path.append('../src')
from readCHS import read_csv
from ios import *
import timeit, functools
import pandas as pd
import filecmp


def CHS():
  
  filePath = r'35.TS.csv'
  CHSPath =  r'CHS_Constants.csv'

  # datetimes = np.arange('2001-06-21', '2001-09-22', np.timedelta64(10, 'm'), dtype='datetime64')

  stations = read_csv(CHSPath)
  stations = stations[np.where(stations['id']==35)]
  # print stations
  # cons = stations['constituents']['name'][0, :]
  cons = stations[0]['constituents'][np.argsort(stations['constituents']['eta'][0, :, 0])[::-1][:10]]['name']
  # cons = np.delete(cons, np.where(cons == 'Z0'))
  # cons = np.array(['Z0', 'M2','S2','MS4'])
  # print cons
  
  ios = IOS(cons,stations,debug=True)
  
  ts = pd.read_csv(filePath,parse_dates=['Datetime'])
  datetime = ts['Datetime'].values
  WL =ts['35'].values
  U =np.zeros(len(WL))
  V =np.zeros(len(WL))
  
 
  ios.extractConstituents(datetime,np.asarray([[WL,U,V]]))
  ios.stations_to_csv('35.Julien.csv')
  
  # ios.to_csv(outPath,type=0, index=None, header=True,date_format="%Y/%m/%d %H:%M:%S")
  
def compare():
  print filecmp.cmp('35.Julien.template.csv','35.Julien.csv')
  
  
if __name__ == "__main__":
  CHS()
  compare()