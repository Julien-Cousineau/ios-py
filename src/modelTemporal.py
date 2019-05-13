import numpy as np
import time
from iosdatabase import getIOSDatabase
from astro import getAstro

def transform(xy,proj):
  import pyproj
  p1 = pyproj.Proj(init=proj)
  p2 = pyproj.Proj(init='epsg:3857')
  fx, fy = pyproj.transform(p1, p2, xy[:, 0], xy[:, 1])
  return np.dstack([fx, fy])[0]

class SWAVE:
  def __init__(self, cons):
    self.newcons,self.newshalls = getIOSDatabase()
    self.setCons(cons)
    self.datetimes = None
    self.stations = None
    self.nstations = None
    self.refslat = None
    
    
  
  def getDWaveType(self):
    return np.dtype([('M_u', 'f8', (self.ndatetime, self.nrc)),
                     ('M_f', 'f8', (self.ndatetime, self.nrc)),
                     ('m_u', 'f8', (self.ndatetime, self.nsc)),
                     ('m_f', 'f8', (self.ndatetime, self.nsc))
                     ])
    
  
  def setCons(self,cons):
    self.cons = cons
    self.rc  = np.where(self.newcons['name'] == self.cons[:, None])[1] # getNewconsIndex
    self.sc  = np.where(self.newshalls['name'] == self.cons[:, None])[1] # getNewshallsIndex
    self.rc1 = np.where(self.newcons['name'] == self.cons[:, None])[0] # getConIndex
    self.sc1 = np.where(self.newshalls['name'] == self.cons[:, None])[0] # getShallowIndex
    self.conNames = self.newcons['name'][np.where(self.newcons['name'] == self.cons[:, None])[1]] # get name of constituents
    self.shallowNames = self.newshalls['name'][np.where(self.newshalls['name'] == self.cons[:, None])[1]] # get name of shallow constituents
    
    self.nrc = nrc = len(self.rc)
    self.nsc = nsc = len(self.sc)
    self.modelDType = [
      ('datetime', 'datetime64[m]'),
      ('dthr', 'f8'),
      ('M_freq', 'f8', (nrc,)),
      ('M_v', 'f8', (nrc,)),
      ('M_uu', 'f8', (nrc, 10)),
      ('m_freq', 'f8', (nsc,)),
      ('m_v', 'f8', (nsc,)),
    ]
  
  def setDatetime(self,datetimes):
    start = time.clock()
    self.datetimes = datetimes
    
    
    rc = self.rc
    sc = self.sc
    newcons = self.newcons[rc]
    newshalls = self.newshalls[sc]
    ndatetime = self.ndatetime = len(datetimes)
    slat = self.refslat
    
    
    model = self.model = np.zeros(ndatetime, dtype=self.modelDType)
    
    model['dthr'] = np.asarray([dt.hour + dt.minute / 60.0 + dt.second / 3600.0 for dt in datetimes.astype(object)])
    
    doods = newcons['dood']
    phase = newcons['phase']
    satDeld = newcons['sats']['deld']
    satPhase = newcons['sats']['phase']
    astro = getAstro(datetimes)

    M_freq = model['M_freq'] = (doods[:, 0] * astro['dtau'][:, np.newaxis] +
                                  doods[:, 1] * astro['ds'][:, np.newaxis] +
                                  doods[:, 2] * astro['dh'][:, np.newaxis] +
                                  doods[:, 3] * astro['dp'][:, np.newaxis] +
                                  doods[:, 4] * astro['dnp'][:, np.newaxis] +
                                  doods[:, 5] * astro['dpp'][:, np.newaxis]) / (24.0 * 365.0)

    M_v = model['M_v'] = np.modf(doods[:, 0] * astro['tau'][:, None] +
                                   doods[:, 1] * astro['s'][:, None] +
                                   doods[:, 2] * astro['h'][:, None] +
                                   doods[:, 3] * astro['p'][:, None] +
                                   doods[:, 4] * astro['np'][:, None] +
                                   doods[:, 5] * astro['pp'][:, None] +
                                   phase[:])[0]

    M_uu = model['M_uu'] = np.modf(satDeld[:, :, 0] * astro['p'][:, None, None] +
                              satDeld[:, :, 1] * astro['np'][:, None, None] +
                              satDeld[:, :, 2] * astro['pp'][:, None, None] +
                              satPhase[:])[0]

    satCorr  = self.satCorr = newcons['sats']['corr']
    satRatio = self.satRatio = newcons['sats']['ratio']

    # Shallow Constituents
    icons = self.icons = newshalls['shallcons']['icon']
    for value in np.arange(0, 45):
      if len(np.where(rc == value)[0]) > 0:
        icons[icons == value] = np.where(rc == value)[0][0]
      else:
        icons[icons == value] = -1
    
    factors = self.factors = newshalls['shallcons']['factor']
    
    model['m_freq'] = np.sum(M_freq[:, icons] * factors, axis=-1)
    model['m_v'] = np.sum(M_v[:, icons] * factors, axis=-1)
    
    
    latLimit = self.latLimit
    waveD = self.waveD = np.zeros(latLimit, dtype=self.getDWaveType())
    p = slat[:, None]
    satShape = (M_freq.shape[0], slat.shape[0]) + satCorr.shape
  
    satAdj = np.zeros(satShape, dtype=np.float64)
  
    satAdj[:, :, satCorr == 0] = satRatio[satCorr == 0]
    satAdj[:, :, satCorr == 1] = satRatio[satCorr == 1] * 0.36309 * (1.0 - 5.0 * p * p) / p
    satAdj[:, :, satCorr == 2] = satRatio[satCorr == 2] * 2.59808 * p
  
    sumc = 1.0 + np.sum(np.swapaxes(satAdj, 0, 1) * np.cos(M_uu * 2.0 * np.pi)[:, :], axis=-1)  # [:,:,:]
    sums = np.sum(np.swapaxes(satAdj, 0, 1) * np.sin(M_uu * 2.0 * np.pi)[:, :], axis=-1)
  
    waveD['M_f'] = np.sqrt((sumc * sumc) + (sums * sums))
    waveD['M_u'] = np.arctan2(sums, sumc) / (2.0 * np.pi)
  
    m_f = np.power(waveD['M_f'][:, :, icons], np.abs(factors))
    m_f[m_f == 0] = 1.0
    waveD['m_f'] = np.prod(m_f, axis=-1)
    waveD['m_u'] = np.sum(waveD['M_u'][:, :, icons] * factors, axis=-1)
    
    print("Created Model in {0}".format(time.clock() - start))

  def setStations(self,stations,proj,latLimit=1000):
    
    nstations = self.nstations = len(stations)
    slat = self.getLatitudes(stations,proj)
    
    if nstations > latLimit:
      self.latLimit=latLimit
      maxY = np.max(slat)
      minY = np.min(slat)
      step = (maxY - minY) / float(latLimit)
      refslat = self.refslat = np.arange(minY, maxY, step)
      self.refslatIndex = np.argmin(np.abs(refslat - slat), axis=1)
    else:
      self.latLimit = nstations
      self.refslat = slat
      self.refslatIndex = np.arange(nstations)
    

  def getLatitudes(self,stations,proj):
    if (proj == 'epsg:3857'):
      latitudes = stations['xy'][:, 1]
    else:
      latitudes = transform(stations['xy'],proj)[:,1]
    return  np.sin((np.pi) * (latitudes / 180.0))
