import numpy as np
import os,sys
import math,time,copy_reg,types
from multiprocessing import Pool
import scipy.linalg as linalg
import pandas as pd


def _pickle_method(m):
  if m.im_self is None:
    return getattr, (m.im_class, m.im_func.func_name)
  else:
    return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)


from modelTemporal import SWAVE
from readCHS import to_csv




class IOS(object):
  def __init__(self,cons,stations,proj='epsg:3857',folder='',npy="dummy.npy",ngroup=1000,debug=False):
    # self.reset()
    # Properties
    self.debug = debug
    self.stations=stations
    self.cons = cons
    self.proj = proj
    self.folder = folder
    self.npy = npy
    self.ngroup = ngroup
    self.npyPath = os.path.join(self.folder,self.npy)
    # self._refstations = None
    
    # self._datetimes=None
    
    self.waveS = SWAVE(self.cons)
    self.waveS.setStations(stations,proj)
    

  
  def getMemoryFile(self):
    npyPath = self.npyPath
    waveS = self.waveS
    shape = (waveS.nstations,3,waveS.ndatetime)
    return np.memmap(npyPath, dtype='float32', mode='r+', shape=shape, order='F')
  
  def createMemoryFile(self):
    npyPath = self.npyPath
    waveS = self.waveS
    shape = (waveS.nstations,3,waveS.ndatetime)
    return np.memmap(npyPath, dtype='float32', mode='w+', shape=shape, order='F')
  
  def getTS(self,datetimes):
    self.datetimes = datetimes
    self.waveS.setDatetime(datetimes)
    self.file = self.createMemoryFile()
    
    if self.debug:
      self._getTS(0)
    else:
      nstations = self.waveS.nstations
      ngroup = self.ngroup
      pool = Pool(processes=1)
      pool.map(self._getTS, range(int(math.ceil(nstations / ngroup))+1))
      pool.close()
      pool.join()
  
    return self.file
    
    
    # return results[0]
  def to_csv(self,outPath,type=0,**kwargs):
    mf = self.getMemoryFile()
    df = pd.DataFrame(data=mf.T[:, type],columns=self.stations['id'])
    df.insert(0,"Datetime",self.datetimes)
    df.to_csv(outPath,**kwargs)
    
  def to_ccsv(self,outPath,type=0,**kwargs):
    stations=self.stations
    to_csv(outPath,stations)    
  
  def getWaveD(self,i):
    istart = self.ngroup * i
    iend = self.ngroup * (i+1)
    waveD = self.waveS.waveD[self.waveS.refslatIndex[istart:iend]] 
    return waveD,waveD.shape[0]
  
  def _getTS(self,i):
    istart = self.ngroup * i
    iend = self.ngroup * (i+1)
    waveS =self.waveS
    ndatetime = waveS.ndatetime
    waveD,size = self.getWaveD(i)
    res=np.zeros((size,3,ndatetime), dtype='float32')
    for j,var in enumerate(['eta','u','v']):
      [M_A,M_P] = self.getRes(i,type='M',value=var)
      [m_A,m_P] = self.getRes(i,type='m',value=var)
      M_res = np.sum(M_A, axis=-1)
      m_res = np.sum(m_A, axis=-1)
      res[:, j, :] = M_res + m_res
    mf = self.getMemoryFile()
    mf[istart:iend]=res
    del mf
    
  def getRes(self,i,type='M',value='eta'):
    waveS = self.waveS
    waveD,size= self.getWaveD(i)
    dthr  = waveS.model['dthr']
    names = waveS.conNames if type == 'M' else waveS.shallowNames
    
    freq  = waveS.model[type + '_freq']
    u     = waveD[type+'_u']
    v     = waveS.model[type+'_v']
    rc    = np.where(self.stations['constituents']['name'][0,:,np.newaxis]== names)[0]
    
    cts   = self.stations['constituents'][value][:, rc]
    amp   = cts[:, :, 0]
    theta = cts[:, :, 1]
    _f     = waveD[type + '_f']
    
    phase  = freq * dthr[:,np.newaxis] + u + v
    revgmt = phase - theta / 360.0
    radgmt = 2.0 * np.pi * np.modf(revgmt)[0]
    
    f      = _f  * amp[:, np.newaxis] 
    A      = f * np.cos(radgmt)
    P      = f * np.sin(radgmt)
    
    return [A,P]
    
  
  def __generateTSConstituents(self,i,res,type="M"):
    [A,P] = res
    waveS = self.waveS
    waveD,size = self.getWaveD(i)
    ndatetime = waveS.ndatetime
    
    if type =="M":
      A = np.delete(A, 0, axis=2) # Remove Zo
      P = np.delete(P, 0, axis=2) # Remove Zo
      all = np.zeros([size,ndatetime, ((waveS.nrc-1) * 2) + 1], dtype=np.float64)
      all[:, :, 0] = 1.0
      all[:, :, 1::2] = A
      all[:, :, 2::2] = P
    else:
      all = np.zeros([size,ndatetime, waveS.nsc * 2], dtype=np.float64)
      all[:, :, ::2] = A
      all[:, :, 1::2] = P
    return all
  
  def generateTSConstituents(self,i):
    # M_A,M_P = self.getRes(i)
    # m_A,m_P = self.getRes(i,type='m')
    M_all = self.__generateTSConstituents(i,self.getRes(i))
    m_all = self.__generateTSConstituents(i,self.getRes(i,type='m'),type='m')
    
    return np.append(M_all,m_all,axis=2)
    
  def extractConstituents(self,datetimes,values):
    
    datetimes = self.datetimes =  datetimes
    self.waveS.setDatetime(datetimes)
    for k, type in enumerate(['eta', 'u', 'v']):
      self.stations['constituents'][type][:,:,  0]= 1.0
      self.stations['constituents'][type][:,:,  1]= 0.0
    stations=self.stations

    M_res = self.generateTSConstituents(0)

    waveS = self.waveS
    conNames = waveS.conNames
    shallowNames = waveS.shallowNames    
    rc = np.where(stations['constituents']['name'][0,:]==conNames[:,np.newaxis])[1]
    sc = np.where(stations['constituents']['name'][0,:]== shallowNames[:,np.newaxis])[1]
    nrc1 = len(rc)

    for j,station in enumerate(stations):
      U, s, Vh = linalg.svd(M_res[j], full_matrices=False)
      M_inv = np.dot(np.dot(Vh.T, linalg.inv(np.diag(s))), U.T)
      
      M_constants = np.einsum('ij,kj->ki', M_inv, values[j])
      tmp_zo=M_constants[:, 0]
      u =np.insert(M_constants[:, 1::2], 0, tmp_zo, axis=1)
      tmp_zo[:]=0.0
      v = np.insert(M_constants[:, 2::2], 0, tmp_zo, axis=1)
      A = np.sqrt(u**2+v**2)
      P = np.arctan2(v , u) * 180 / np.pi
      P[P < 0] += 360.0 # [-180,180] ~> [0,360]
      

      for k,type in enumerate(['eta','u','v']):
        station['constituents'][type][rc, 0]=A[k][:nrc1]
        station['constituents'][type][rc, 1]=P[k][:nrc1]
        station['constituents'][type][sc, 0]=A[k][nrc1:]
        station['constituents'][type][sc, 1]=P[k][nrc1:]
   
    
    # Remove default amplitude of 1.0
    index =  np.where(self.stations['constituents']['name'][0,:] == np.append(conNames,shallowNames)[:,np.newaxis])[1]
    index =  np.delete(np.arange(len(self.stations['constituents']['name'][0,:])),index)
    for k, type in enumerate(['eta', 'u', 'v']):
      stations['constituents'][type][:,index,  0]= 0.0
