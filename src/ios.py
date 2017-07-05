import numpy as np

import math,time,functools,copy_reg,types
from multiprocessing import Pool




def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)


class IOS:
    # ------------------------------------------------------------------------------------------------------------------
    # Types
    def __init__(self):
        self.astroType = np.dtype([('d1', 'f8'),
                              ('d12', 'f8'),
                              ('d2', 'f8'),
                              ('d22', 'f8'),
                              ('d23', 'f8'),
                              ('f', 'f8'),
                              ('f2', 'f8'),
                              ('h', 'f8'),
                              ('pp', 'f8'),
                              ('s', 'f8'),
                              ('p', 'f8'),
                              ('np', 'f8'),
                              ('dh', 'f8'),
                              ('dpp', 'f8'),
                              ('ds', 'f8'),
                              ('dp', 'f8'),
                              ('dnp', 'f8'),
                              ('hh', 'f8'),
                              ('tau', 'f8'),
                              ('dtau', 'f8')])
        self.constituentType = np.dtype([('name', '|S8'),
                                    ('eta', 'f8', (1, 2)),
                                    ('u', 'f8', (1, 2)),
                                    ('v', 'f8', (1, 2)),
                                    ('r', 'f8'),
                                    ('zeta', 'f8'),
                                    ('konans', 'i4')])

        self.satType = np.dtype([('deld', '3i4'),
                            ('phase', 'f8'),
                            ('ratio', 'f8'),
                            ('corr', 'i4'),
                            ('adj', 'f8')])

        self.shallconType = np.dtype([('icon', 'i4'),
                                 ('factor', 'f8')])

        self.newconType = np.dtype([('name', '|S8'),
                               ('dood', '6f8'),
                               ('nsats', 'i4'),
                               ('sats', self.satType, 10),
                               ('phase', 'f8'),
                               ('f', 'f8'),
                               ('u', 'f8'),
                               ('v', 'f8'),
                               ('freq', 'f8')])

        self.conRType = np.dtype([('phase', 'f8'),
                             ('uu', 'f8', 5),
                             ('sumc', 'f8'),
                             ('sums', 'f8'),
                             ('f', 'f8'),
                             ('u', 'f8'),
                             ('v', 'f8'),
                             ('freq', 'f8')])
        
        self.tempType = np.dtype([('ajd', 'f8')])

        self.stationType = np.dtype([('name', '|S32'),
                                ('id', 'i4'),
                                ('xy', '2f8',),
                                ('constituents', self.constituentType, 200)])

        self.newcons, self.newshalls = self.getIOS()
    # ------------------------------------------------------------------------------------------------------------------
    # Functions for types (Initialization)
    def getNewshallType(self):
        return np.dtype([('name', '|S8'),
                         ('ncon', 'i4'),
                         ('shallcons', self.shallconType, 5),
                         ('f', 'f8'),
                         ('u', 'f8'),
                         ('v', 'f8'),
                         ('freq', 'f8')])

    def getSWaveType(self,ndatetime, nMc, nmc):
        return np.dtype([('ndatetime', 'i4'),
                         ('datetime', 'datetime64[m]', ndatetime),
                         ('dthr', 'f8', ndatetime),
                         ('rc', 'i4', (nMc, 1)),
                         ('sc', 'i4', (nmc, 1)),
                         ('rc1', 'i4', (nMc, 1)),
                         ('sc1', 'i4', (nmc, 1)),
                         ('M_freq', 'f8', (ndatetime, nMc)),
                         ('M_v', 'f8', (ndatetime, nMc)),
                         ('M_uu', 'f8', (ndatetime, nMc, 10)),
                         ('satCorr', 'i4', (nMc, 10)),
                         ('satRatio', 'f8', (nMc, 10)),
                         ('factors', 'f8', (nmc, 5)),
                         ('icons', 'i4', (nmc, 5)),
                         ('m_freq', 'f8', (ndatetime, nmc)),
                         ('m_v', 'f8', (ndatetime, nmc)),
                         ])

    def getDWaveType(self,ndatetime, nMc, nmc):
        return np.dtype([('M_u', 'f8', (ndatetime, nMc)),
                         ('M_f', 'f8', (ndatetime, nMc)),
                         ('m_u', 'f8', (ndatetime, nmc)),
                         ('m_f', 'f8', (ndatetime, nmc))
                         ])

    def getPSType(self,nMc, nmc):
        return np.dtype([('M_A', 'f8', (nMc, 1)),
                         ('M_phi', 'f8', (nMc, 1)),
                         ('m_A', 'f8', (nmc, 1)),
                         ('m_phi', 'f8', (nmc, 1))
                         ])

    # ------------------------------------------------------------------------------------------------------------------
    # Functions
    def getAstro(self,datetimes):
        # Description here...

    
        _d1 = datetimes - np.datetime64('1899-12-31T12:00')
        d1 = _d1.astype('timedelta64[s]') / np.timedelta64(1, 'D') # d1 [1-D Array of Datetimes]
    
        _h1 = datetimes - np.datetime64('1899-12-31T00:00')
        h1 = _h1.astype('timedelta64[s]') / np.timedelta64(1, 'D')


        nDatetimes = len(d1)
        astro = np.zeros(nDatetimes, dtype=self.astroType)
    
        astro['d1'] = d1
        astro['d12'] = d1 * d1
        astro['d2'] = d1 * 1.0E-04
        astro['d22'] = astro['d2'] * astro['d2']
        astro['d23'] = astro['d2'] ** 3.0
        astro['f'] = 360.0
        astro['f2'] = astro['f'] #/ 365.0
        astro['h'] = np.modf((2.79696678E+02 + d1 * 9.856473354E-01 + astro['d22'] * 2.267E-05) / astro['f'])[0]
        astro['pp'] = \
        np.modf((2.81220844E+02 + d1 * 4.70684E-05 + astro['d22'] * 3.39E-05 + astro['d23'] * 7.0E-08) / astro['f'])[0]
        astro['s'] = \
        np.modf((2.70434164E+02 + d1 * 1.31763965268E+01 - astro['d22'] * 8.5E-05 + astro['d23'] * 3.9E-08) / astro['f'])[0]
        astro['p'] = \
        np.modf((3.34329556E+02 + d1 * 1.114040803E-01 - astro['d22'] * 7.739E-04 - astro['d23'] * 2.6E-07) / astro['f'])[0]
        astro['np'] = \
        np.modf((-2.59183275E+02 + d1 * 5.29539222E-02 - astro['d22'] * 1.557E-04 - astro['d23'] * 5.0E-08) / astro['f'])[0]
        astro['dh'] = (9.856473354E-01 + d1 * 2.267E-05 * 2.0E-08) / astro['f2']
        astro['dpp'] = (4.70684E-05 + d1 * 3.39E-05 * 2.0E-08 + astro['d12'] * 7.0E-08 * 3.0E-12) / astro['f2']
        astro['ds'] = (1.31763965268E+01 - d1 * 8.5E-05 * 2.0E-08 + astro['d12'] * 3.9E-08 * 3.0E-12) / astro['f2']
        astro['dp'] = (1.114040803E-01 - d1 * 7.739E-04 * 2.0E-08 - astro['d12'] * 2.6E-07 * 3.0E-12) / astro['f2']
        astro['dnp'] = (5.29539222E-02 - d1 * 1.557E-04 * 2.0E-08 - astro['d12'] * 5.0E-08 * 3.0E-12) / astro['f2']
    
        # _ktmp = kh.replace(year=kh.year + 1) - datetime(1, 1, 1)
        # ktmp = float(_ktmp.days) + float(_ktmp.seconds) / 86400.0
        astro['hh'] = np.modf(h1)[0]
        astro['tau'] = astro['hh'] / 24.0 + astro['h'] - astro['s']
        astro['dtau'] = 1 + astro['dh'] - astro['ds']
    
        return astro


    def getIOS(self):
        # Description here...

        filePath = r'IOS_tidetbl.txt'
        with open(filePath, 'r') as file:
            text = []
            for line in file:
                text.append(line)
    
            newcons = []
            newshalls = []
            con = True
            icount = 0
            for i in range(0, len(text)):
                line = text[i]
                readline = filter(None, line.replace('\r', ' ').
                           replace('\n', '').split(' '))
                if not len(readline) == 0:
                    if con:
                        newcons, icount = self.getNewcon(readline, newcons, icount)
                    else:
                        newshalls = self.getNewshall(readline, newcons, newshalls)
                else:
                    con = False
    
            newcons = np.asarray(newcons)
            newshalls = np.asarray(newshalls)
    
            return newcons, newshalls
    
    def getNewcon(self,readline, newcons, icount):
        # Description here...

        temp = []
        for s in readline:
            if s != "":
                temp.append(s)
        if len(temp) != 0:
            name = temp[0]
            _newcons = [newcon for newcon in newcons if newcon['name'] == name]
            if len(_newcons) == 0:
                nsats = int(readline[8])
                newcon = np.zeros(1, dtype=self.newconType)
                newcon['name'] = name
                newcon['dood'] = [float(temp[1]), float(temp[2]), float(temp[3]), 
                                  float(temp[4]), float(temp[5]), float(temp[6])]
                newcon['phase'] = float(temp[7])
                newcon['nsats'] = nsats
                newcons.append(newcon)
                icount = 0
            else:
    
                for i in range(1, len(temp), 5):
                    sat = np.zeros(1, dtype=self.satType)
                    sat['deld'] = [int(float(temp[i])), int(float(temp[i + 1])), 
                                                        int(float(temp[i + 2]))]
                    sat['phase'] = float(temp[i + 3])
    
                    if "R1" in temp[i + 4]:
                        sat['corr'] = 1
                        sat['ratio'] = float(temp[i + 4].split("R1")[0])
                    elif "R2" in temp[i + 4]:
                        sat['corr'] = 2
                        sat['ratio'] = float(temp[i + 4].split("R2")[0])
                    else:
                        sat['corr'] = 0
                        sat['ratio'] = float(temp[i + 4])
    
                    newcon = _newcons[0]
                    newcon['sats'][0][icount] = sat[0]
                    icount += 1
        return newcons, icount
    
    def getNewshall(self,readline, newcons, newshalls):
        # Description here...

        name = readline[0]
        ncon = int(readline[1])
    
        newshall = np.zeros(1, dtype=self.getNewshallType())
        newshall['name'] = name
        newshall['ncon'] = ncon
    
        icount = 0
        for i in range(2, len(readline), 2):
            factor = readline[i]
            conname = readline[i + 1]
            if factor != "" and conname != "":
                inewcon = [newcon['name'][0] for newcon in newcons].index(conname)
                shallcon = np.zeros(1, dtype=self.shallconType)
                shallcon['factor'] = float(factor)
                shallcon['icon'] = inewcon
                newshall['shallcons'][0, icount] = shallcon[0]
                icount += 1
    
        newshalls.append(newshall)
    
        return newshalls


    def getRefStations(self,minLat, maxLat, nstep, cons):
        # Description here...
        # Function to reduce # of station

        step = (maxLat - minLat) / float(nstep)
        latitudes = np.arange(minLat, maxLat, step)
        refstations = np.zeros(nstep, dtype=self.stationType)
        for istation, lat in enumerate(latitudes):
            refstations['name'][istation] = str(lat)
            refstations['id'][istation] = istation
            refstations['xy'][istation] = (0, lat)
            for icon, con_name in enumerate(cons):
                refstations['constituents']['name'][istation, icon] = con_name
                refstations['constituents']['eta'][istation, icon] = (1.0, 0.0)
        return refstations

    def getClosestRefStation(self,refstations, stations):
        idx = np.argmin(np.abs(refstations['xy'][:, 1] - stations['xy'][:, 1, np.newaxis]), axis=1)
        return idx

    # ------------------------------------------------------------------------------------------------------------------
    # Functions (Generate Time-Series)
    def Run(self,datetimes, stations,cons):
        waveS,refstations,idx = self._generateTimeV(datetimes,stations,cons)

        start = time.clock()
        waveD = self.generateSpaceV(waveS, refstations)
        print("generateSpaceV {0}".format(time.clock() - start))

        start = time.clock()
        stationids = stations['id']
        values =  self.generateTimeSpace(waveS, waveD, stations, idx)
        print("generateTimeSpace {0}".format(time.clock() - start))
        with open(r'temp_0.csv', 'w') as f:
            self.writeCSV(f,stationids,values)

    def RunP(self, nprocessor,datetimes, stations,cons):
        waveS, refstations, idx = self._generateTimeV(datetimes, stations,cons)

        start = time.clock()
        p = Pool(nprocessor)
        t = functools.partial(self._pGenerateSpaceV, waveS=waveS, stations=refstations)
        nrefstations = len(refstations)
        waveD = np.asarray(p.map(t, xrange(nrefstations))).ravel()
        print("generateSpaceV {0}".format(time.clock() - start))
        p.close()
        p.join()


        p = Pool(nprocessor)

        tt = functools.partial(self._pGenerateTimeSpace, waveS=waveS, waveD=waveD, stations=stations, idx=idx)

        nrows = int(math.ceil(len(stations) / 100))
        divider=10

        nrowspergroup=int(math.ceil(nrows/divider))
        with open(r'temp_10.csv', 'w') as f:
            for i in range(divider):
                start = time.clock()
                istart=i*nrowspergroup
                iend =istart + nrowspergroup
                results = p.map(tt, xrange(istart,iend))
                print("generateTimeSpace {0}".format(time.clock() - start))
                for result in results:
                    stationids=result[0]
                    values = result[1]
                    self.writeCSV(f, stationids, values)


    def writeCSV(self,f,stationids,values):
        for i in range(len(stationids)):
            id = stationids[i]
            ts = values[i]
            f.write(str(id) + ",")
            f.write(",".join(map("{:.3f}".format, ts)))
            f.write("\n")

    def _generateTimeV(self,datetimes,stations,cons):
        start = time.clock()
        waveS = self.generateTimeV(datetimes, cons)
        print("generateTimeV {0}".format(time.clock() - start))

        nstations = len(stations)
        if nstations > 10:
            maxY = np.max(stations['xy'][:, 1])
            minY = np.min(stations['xy'][:, 1])
            refstations = self.getRefStations(minY, maxY, 10, cons)
        else:
            refstations = stations

        idx = self.getClosestRefStation(refstations, stations)
        return waveS,refstations,idx

    def _pGenerateSpaceV(self,x, waveS, stations):
        pstations = stations[x:x + 1]
        return self.generateSpaceV(waveS, pstations)

    def _pGenerateTimeSpace(self,x, waveS, waveD, stations, idx):
        i = 100 * x
        j = 100 * (x + 1)
        if j > len(stations):
            j = len(stations)
        pstations = stations[i:j]
        pidx = idx[i:j]

        values = self.generateTimeSpace(waveS, waveD, pstations, pidx)



        return [pstations['id'],values]

    def generateTimeV(self, datetimes, cons):
        newcons = self.newcons
        newshalls = self.newshalls
        astro = self.getAstro(datetimes)


        b = newcons['name'][:, 0]
        c = newshalls['name'][:, 0]
        rc = np.where(b == cons[:, None])[1]
        sc = np.where(c == cons[:, None])[1]
        rc1 = np.where(b == cons[:, None])[0]
        sc1 = np.where(c == cons[:, None])[0]

        ndatetime = len(datetimes)
        nMc = len(rc)
        nmc = len(sc)

        waveS = np.zeros(1, dtype=self.getSWaveType(ndatetime, nMc, nmc))
        waveS['ndatetime'][0] = ndatetime
        waveS['datetime'][0] = datetimes
        waveS['dthr'][0] = np.asarray(
            [dt.hour + dt.minute / 60.0 + dt.second / 3600.0 for dt in datetimes.astype(object)])
        waveS['rc'][0][:, 0] = rc
        waveS['sc'][0][:, 0] = sc
        waveS['rc1'][0][:, 0] = rc1
        waveS['sc1'][0][:, 0] = sc1

        # Major Constituents
        doods = newcons['dood'][rc, 0]
        phase = newcons['phase'][rc, 0]
        satDeld = newcons['sats']['deld'][rc, 0]
        satPhase = newcons['sats']['phase'][rc, 0]

        waveS['M_freq'][0] = (doods[:, 0] * astro['dtau'][:, None] +
                              doods[:, 1] * astro['ds'][:, None] +
                              doods[:, 2] * astro['dh'][:, None] +
                              doods[:, 3] * astro['dp'][:, None] +
                              doods[:, 4] * astro['dnp'][:, None] +
                              doods[:, 5] * astro['dpp'][:, None]) / 24.0

        waveS['M_v'][0] = np.modf(doods[:, 0] * astro['tau'][:, None] +
                                  doods[:, 1] * astro['s'][:, None] +
                                  doods[:, 2] * astro['h'][:, None] +
                                  doods[:, 3] * astro['p'][:, None] +
                                  doods[:, 4] * astro['np'][:, None] +
                                  doods[:, 5] * astro['pp'][:, None] +
                                  phase[:])[0]

        waveS['M_uu'][0] = np.modf(satDeld[:, :, 0] * astro['p'][:, None, None] +
                                   satDeld[:, :, 1] * astro['np'][:, None, None] +
                                   satDeld[:, :, 2] * astro['pp'][:, None, None] +
                                   satPhase[:])[0]

        waveS['satCorr'] = newcons['sats']['corr'][rc, 0]
        waveS['satRatio'] = newcons['sats']['ratio'][rc, 0]

        # Minor Constituents
        waveS['icons'] = newshalls['shallcons']['icon'][sc, 0]
        for value in np.arange(0, 45):
            if len(np.where(rc == value)[0]) > 0:
                waveS['icons'][waveS['icons'] == value] = np.where(rc == value)[0][0]
            else:
                waveS['icons'][waveS['icons'] == value] = -1
        waveS['factors'] = newshalls['shallcons']['factor'][sc, 0]
        waveS['m_freq'][0] = np.sum(waveS['M_freq'][0][:, waveS['icons'][0]] * waveS['factors'][0], axis=-1)
        waveS['m_v'][0] = np.sum(waveS['M_v'][0][:, waveS['icons'][0]] * waveS['factors'][0], axis=-1)
        return waveS

    def generateSpaceV(self,waveS, refstations):
        ndatetime = waveS['ndatetime'][0]
        rc1 = waveS['rc1'][0][:, 0]
        sc1 = waveS['sc1'][0][:, 0]
        nMc = len(rc1)
        nmc = len(sc1)
        npoint = len(refstations)
        waveD = np.zeros(npoint, dtype=self.getDWaveType(ndatetime, nMc, nmc))

        slat = np.sin((np.pi) * (refstations['xy'][:, 1] / 180.0))
        p = slat[:, None]

        satShape = (waveS['M_freq'][0].shape[0], slat.shape[0],) + waveS['satCorr'][0].shape
        satAdj = np.zeros(satShape, dtype=np.float64)
        satAdj[:, :, waveS['satCorr'][0] == 0] = waveS['satRatio'][0][waveS['satCorr'][0] == 0]
        satAdj[:, :, waveS['satCorr'][0] == 1] = waveS['satRatio'][0][waveS['satCorr'][0] == 1] * 0.36309 * (
        1.0 - 5.0 * p * p) / p
        satAdj[:, :, waveS['satCorr'][0] == 2] = waveS['satRatio'][0][waveS['satCorr'][0] == 2] * 2.59808 * p

        sumc = 1.0 + np.sum(np.swapaxes(satAdj, 0, 1) * np.cos(waveS['M_uu'] * 2.0 * np.pi)[:, :], axis=-1)  # [:,:,:]
        sums = np.sum(np.swapaxes(satAdj, 0, 1) * np.sin(waveS['M_uu'] * 2.0 * np.pi)[:, :], axis=-1)

        waveD['M_f'] = np.sqrt((sumc * sumc) + (sums * sums))
        waveD['M_u'] = np.arctan2(sums, sumc) / (2.0 * np.pi)

        m_f = np.power(waveD['M_f'][:, :, waveS['icons'][0]], np.abs(waveS['factors'][0]))
        m_f[m_f == 0] = 1.0
        waveD['m_f'] = np.prod(m_f, axis=-1)
        waveD['m_u'] = np.sum(waveD['M_u'][:, :, waveS['icons'][0]] * waveS['factors'][0], axis=-1)
        return waveD

    def generateTimeSpace(self,waveS, waveD, stations, idx):

        rc1 = waveS['rc1'][0][:, 0]
        sc1 = waveS['sc1'][0][:, 0]
        nMc = len(rc1)
        nmc = len(sc1)
        nstations = len(stations)

        waveTS = np.zeros(nstations, dtype=self.getPSType(nMc, nmc))
        waveTS['M_A'][:, :, 0] = stations['constituents']['eta'][:, rc1, 0, 0]
        waveTS['M_phi'][:, :, 0] = stations['constituents']['eta'][:, rc1, 0, 1]
        waveTS['m_A'][:, :, 0] = stations['constituents']['eta'][:, sc1, 0, 0]
        waveTS['m_phi'][:, :, 0] = stations['constituents']['eta'][:, sc1, 0, 1]

        M_revgmt = waveS['M_freq'][0] * waveS['dthr'][0][:, None] + waveD['M_u'][idx] + waveS['M_v'][0] - waveTS['M_phi'][
                                                                                                          :, :, 0][:,
                                                                                                          np.newaxis] / 360.0
        M_radgmt = 2.0 * np.pi * np.modf(M_revgmt)[0]
        M_res = np.sum(waveD['M_f'][idx] * waveTS['M_A'][:, :, 0][:, np.newaxis] * np.cos(M_radgmt), axis=-1)
        M_res2 = np.sum(waveD['M_f'][idx] * waveTS['M_A'][:, :, 0][:, np.newaxis] * np.sin(M_radgmt), axis=-1)

        m_revgmt = waveS['m_freq'][0] * waveS['dthr'][0][:, None] + waveD['m_u'][idx] + waveS['m_v'][0] - waveTS[
                                                                                                              'm_phi'][
                                                                                                          :, :, 0][:,
                                                                                                          np.newaxis] / 360.0
        m_radgmt = 2.0 * np.pi * np.modf(m_revgmt)[0]
        m_res = np.sum(waveD['m_f'][idx] * waveTS['m_A'][:, :, 0][:, np.newaxis] * np.cos(m_radgmt), axis=-1)
        m_res2 = np.sum(waveD['m_f'][idx] * waveTS['m_A'][:, :, 0][:, np.newaxis] * np.sin(m_radgmt), axis=-1)

        ress = M_res + m_res - waveTS['M_A'][:, 0]
        res2s = M_res2 + m_res2 - waveTS['M_A'][:, 0]
        return ress


    def extractConstituents(self,datetimes,stations,cons,):
        waveS, refstations, idx = self._generateTimeV(datetimes, stations, cons)
        waveD = self.generateSpaceV(waveS, refstations)
        values = self.generateTimeSpace(waveS, waveD, stations, idx)
