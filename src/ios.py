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
    
  

  def getPSType(self,nMc, nmc):
    return np.dtype([('M', 'f8', (nMc, 2)),
                     ('m', 'f8', (nmc, 2)),
                     ])
  

  
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
  
  def _getTS(self,i):
    istart = self.ngroup * i
    iend = self.ngroup * (i+1)
    waveS =self.waveS
    stations = self.stations
    ndatetime = waveS.ndatetime
    M_freq = waveS.model['M_freq']
    m_freq = waveS.model['m_freq']
    dthr = waveS.model['dthr']
    M_v = waveS.model['M_v']
    m_v = waveS.model['m_v']
    
    waveD = waveS.waveD[self.waveS.refslatIndex[istart:iend]]
    conNames = waveS.conNames
    shallowNames = waveS.shallowNames
    
    rc = np.where(self.stations['constituents']['name'][0,:,np.newaxis]== conNames)[0]
    sc = np.where(self.stations['constituents']['name'][0, :, np.newaxis] == shallowNames)[0]
    
    size = waveD.shape[0] 
    waveTS = np.zeros(size, dtype=self.getPSType(waveS.nrc, waveS.nsc))
    
    res_1=np.zeros((size,3,ndatetime), dtype='float32')
    
    MPhase = M_freq *dthr[:,np.newaxis] + waveD['M_u'] + M_v
    mPhase = m_freq * dthr[:,np.newaxis] + waveD['m_u']  + m_v
    
    for i,var in enumerate(['eta','u','v']):
      waveTS['M'] = stations['constituents'][var][:, rc]
      waveTS['m'] = stations['constituents'][var][:, sc]
      
      M_revgmt = MPhase - waveTS['M'][:, :, 1][:,np.newaxis] / 360.0
      M_radgmt = 2.0 * np.pi * np.modf(M_revgmt)[0]
      M_res = np.sum(waveD['M_f'] * waveTS['M'][:, :, 0][:, np.newaxis] * np.cos(M_radgmt), axis=-1)
      
      m_revgmt = mPhase  - waveTS['m'][:, :, 1][:,np.newaxis] / 360.0
      m_radgmt = 2.0 * np.pi * np.modf(m_revgmt)[0]
      m_res= np.sum(waveD['m_f'] * waveTS['m'][:, :, 0][:, np.newaxis] * np.cos(m_radgmt), axis=-1)
      
      res_1[:, i, :] = M_res + m_res
        
    mf = self.getMemoryFile()
    mf[istart:iend]=res_1
    del mf
  
  
  def getRes(self,i,type='M'):
    istart = self.ngroup * i
    iend = self.ngroup * (i+1)    
    waveS = self.waveS
    names = waveS.conNames if type == 'M' else waveS.shallowNames
    msize  = waveS.nrc if type == 'M' else waveS.nsc
    rc = np.where(self.stations['constituents']['name'][0,:]==names[:,np.newaxis])[0]
    waveD = waveS.waveD[self.waveS.refslatIndex[istart:iend]]
    size = waveD.shape[0] 
    cts = self.stations['constituents']['eta'][:, rc]
    MPhase = M_freq * dthr[:,np.newaxis] + waveD[type+'_u'] + waveS.model[type+'_v']
    
  
  def generateTSConstituents(self,i):
    istart = self.ngroup * i
    iend = self.ngroup * (i+1)
    waveS = self.waveS
    M_freq = waveS.model['M_freq']
    m_freq = waveS.model['m_freq']
    dthr = waveS.model['dthr']
    M_v = waveS.model['M_v']
    m_v = waveS.model['m_v']
    
    ndatetime = waveS.ndatetime
    
    conNames = waveS.conNames
    shallowNames = waveS.shallowNames

    rc = np.where(self.stations['constituents']['name'][0,:]==conNames[:,np.newaxis])[0]
    sc = np.where(self.stations['constituents']['name'][0,:]==shallowNames[:,np.newaxis])[0]

    waveD = waveS.waveD[self.waveS.refslatIndex[istart:iend]]
    size = waveD.shape[0] 
    waveTS = np.zeros(size, dtype=self.getPSType(waveS.nrc, waveS.nsc))

    waveTS['M'] = self.stations['constituents']['eta'][:, rc]
    waveTS['m'] = self.stations['constituents']['eta'][:, sc]
    
    MPhase = M_freq * dthr[:,np.newaxis] + waveD['M_u'] + M_v
    mPhase = m_freq * dthr[:,np.newaxis] + waveD['m_u'] + m_v
    
    M_revgmt = MPhase - waveTS['M'][:, :, 1][:,np.newaxis] / 360.0
    M_radgmt = 2.0 * np.pi * np.modf(M_revgmt)[0]
    
    M_resA = waveD['M_f'] * waveTS['M'][:, :, 0][:, np.newaxis] * np.cos(M_radgmt)
    M_resP = waveD['M_f'] * waveTS['M'][:, :, 0][:, np.newaxis] * np.sin(M_radgmt)
    M_resA = np.delete(M_resA, 0, axis=2) # Remove Zo
    M_resP = np.delete(M_resP, 0, axis=2) # Remove Zo
  
    
    M_all = np.zeros([size,ndatetime, ((waveS.nrc-1) * 2) + 1], dtype=np.float64)
    M_all[:, :, 0] = 1.0
    M_all[:, :, 1::2] = M_resA
    M_all[:, :, 2::2] = M_resP
    
    if(waveS.nsc>0):
      m_revgmt = mPhase  - waveTS['m'][:, :, 1][:,np.newaxis] / 360.0
      m_radgmt = 2.0 * np.pi * np.modf(m_revgmt)[0]
      m_resA = waveD['m_f'] * waveTS['m'][:, :, 0][:, np.newaxis] * np.cos(m_radgmt)
      m_resP = waveD['m_f'] * waveTS['m'][:, :, 0][:, np.newaxis] * np.sin(m_radgmt)
      
      m_all = np.zeros([size,ndatetime, waveS.nsc * 2], dtype=np.float64)
      m_all[:, :, ::2] = m_resA
      m_all[:, :, 1::2] = m_resP
      M_all=np.append(M_all,m_all,axis=2)
    
    
    # print M_all
    return M_all
  
  
  
  
  def resetConstants(self):
    for k, type in enumerate(['eta', 'u', 'v']):
      self.stations['constituents'][type][:,:,  0]= 1.0
      self.stations['constituents'][type][:,:,  1]= 0.0

  def extractConstituents(self,datetimes,values):
    
    datetimes = self.datetimes =  datetimes
    self.waveS.setDatetime(datetimes)
    self.resetConstants()
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
      M_u =np.insert(M_constants[:, 1::2], 0, tmp_zo, axis=1)
      tmp_zo[:]=0.0
      M_v = np.insert(M_constants[:, 2::2], 0, tmp_zo, axis=1)
      M_A = np.sqrt(M_u**2+M_v**2)
      M_P = np.arctan2(M_v , M_u) * 180 / np.pi
      M_P[M_P < 0] += 360.0 # [-180,180] ~> [0,360]
      

      for k,type in enumerate(['eta','u','v']):
        # print station['constituents'][type][rc, 0].shape
        
        station['constituents'][type][rc, 0]=M_A[k][:nrc1]
        station['constituents'][type][rc, 1]=M_P[k][:nrc1]
        station['constituents'][type][sc, 0]=M_A[k][nrc1:]
        station['constituents'][type][sc, 1]=M_P[k][nrc1:]
   
    
    rrr = np.where(self.stations['constituents']['name'][0,:] == np.append(conNames,shallowNames)[:,np.newaxis])[1]
    www =  np.delete(np.arange(len(self.stations['constituents']['name'][0,:])),rrr)
    for k, type in enumerate(['eta', 'u', 'v']):
      stations['constituents'][type][:,www,  0]= 0.0
   
