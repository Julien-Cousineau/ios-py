import numpy as np

class IOS:
    # --------------------------------------------------------------------------
    # Types
    def __init__(self):
        self.satType = np.dtype([('deld', '3i4'),
                            ('phase', 'f8'),
                            ('ratio', 'f8'),
                            ('corr', 'i4'),
                            ('adj', 'f8')])
        
        self.newconType = np.dtype([('name', '|S8'),
                               ('dood', '6f8'),
                               ('nsats', 'i4'),
                               ('sats', self.satType, 10),
                               ('phase', 'f8'),
                               ('f', 'f8'),
                               ('u', 'f8'),
                               ('v', 'f8'),
                               ('freq', 'f8')])
        
        self.shallconType = np.dtype([('icon', 'i4'),
                                 ('factor', 'f8')])            
        
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
        
        self.conRType = np.dtype([('phase', 'f8'),
                             ('uu', 'f8', 5),
                             ('sumc', 'f8'),
                             ('sums', 'f8'),
                             ('f', 'f8'),
                             ('u', 'f8'),
                             ('v', 'f8'),
                             ('freq', 'f8')])
        
        self.tempType = np.dtype([('ajd', 'f8')])
        
        self.newcons, self.newshalls = self.getIOS()
    # --------------------------------------------------------------------------
    # Functions (Initialization)
    def getNewshallType(self,ncon):
        return np.dtype([('name', '|S8'),
                         ('ncon', 'i4'),
                         ('shallcons', self.shallconType, 5),
                         ('f', 'f8'),
                         ('u', 'f8'),
                         ('v', 'f8'),
                         ('freq', 'f8')])
                         
    def getWaveType(self,npoint,ndatetime,nMc):
        return np.dtype([('freq', 'f8',(ndatetime,nMc)),
                            ('u', 'f8',(npoint,ndatetime,nMc)),
                            ('v', 'f8',(ndatetime,nMc)),
                            ('f', 'f8', (npoint,ndatetime, nMc)),
                            ('A', 'f8', nMc),
                            ('phi', 'f8', nMc)])
    
    def getMinorWaveType(self,npoint,ndatetime,nMc):
        return np.dtype([('m_freq', 'f8',(ndatetime,nMc)),
                            ('m_u', 'f8',(npoint,ndatetime,nMc)),
                            ('m_v', 'f8',(ndatetime,nMc)),
                            ('m_f', 'f8', (npoint,ndatetime, nMc)),
                            ('m_A', 'f8', nMc),
                            ('m_phi', 'f8', nMc)])

    def getAstro(self,datetimes):
        # d1 [1-D Array of Datetimes]
    
        _d1 = datetimes - np.datetime64('1899-12-31T12:00')
        d1 = _d1.astype('timedelta64[s]') / np.timedelta64(1, 'D')
    
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
        name = readline[0]
        ncon = int(readline[1])
    
        newshall = np.zeros(1, dtype=self.getNewshallType(ncon))
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
    
    # Functions (Generate Time-Series)
    def pGenerateTS(self,x,datetimes, stations):
        pstations = stations[x:x + 1]
        self.generateTS(datetimes, pstations)

    def generateTS(self,datetimes, stations):
        newcons = self.newcons
        newshalls = self.newshalls
        astro = self.getAstro(datetimes)
        dthr = np.asarray([dt.hour + dt.minute / 60.0 + dt.second / 3600.0 for dt in datetimes.astype(object)])
    
        slat = np.sin((np.pi) * (stations['xy'][:,1] / 180.0))
    
        a = stations['constituents']['name'][0, :]
        b = newcons['name'][:, 0]
        c = newshalls['name'][:, 0]
        rc = np.where(b == a[:, None])[1]
        sc = np.where(c == a[:, None])[1]
        rc1 = np.where(b == a[:, None])[0]
        sc1 = np.where(c == a[:, None])[0]
    
        M_A = stations['constituents']['eta'][:,rc1,0,0][0]
        M_phi = stations['constituents']['eta'][:,rc1,0,1][0]
        m_A = stations['constituents']['eta'][:,sc1,0,0][0]
        m_phi = stations['constituents']['eta'][:,sc1,0,1][0]
    
    
        npoint = len(slat)
        ndatetime = len(datetimes)
        nMc = len(rc)
        nmc = len(rc)
        
        print nMc
        print ndatetime
        print npoint
        # M_wave = np.zeros(1, dtype=self.getWaveType(npoint,ndatetime,nMc))
        # m_wave = np.zeros(1, dtype=self.getMinorWaveType(npoint, ndatetime, nmc))
        
        
        # Major Constituents
        doods = newcons['dood'][rc, 0]
        phase = newcons['phase'][rc, 0]
        satDeld = newcons['sats']['deld'][rc, 0]
        satPhase = newcons['sats']['phase'][rc, 0]
    
        #M_freq = M_wave['freq'][0]
        M_freq = (doods[:, 0] * astro['dtau'][:, None] +
                  doods[:, 1] * astro['ds'][:, None] +
                  doods[:, 2] * astro['dh'][:, None] +
                  doods[:, 3] * astro['dp'][:, None] +
                  doods[:, 4] * astro['dnp'][:, None] +
                  doods[:, 5] * astro['dpp'][:, None]) / 24.0
    
        M_v = np.modf(doods[:, 0] * astro['tau'][:, None] +
                      doods[:, 1] * astro['s'][:, None] +
                      doods[:, 2] * astro['h'][:, None] +
                      doods[:, 3] * astro['p'][:, None] +
                      doods[:, 4] * astro['np'][:, None] +
                      doods[:, 5] * astro['pp'][:, None] +
                      phase[:])[0]
    
        uu = np.modf(satDeld[:, :, 0] * astro['p'][:, None, None] +
                     satDeld[:, :, 1] * astro['np'][:, None, None] +
                     satDeld[:, :, 2] * astro['pp'][:, None, None] +
                     satPhase[:])[0]
    
        satCorr = newcons['sats']['corr'][rc, 0]
        # satAdj = newcons['sats']['adj'][rc,0]
        satRatio = newcons['sats']['ratio'][rc, 0]
    
        satShape = (M_freq.shape[0], slat.shape[0],) + satCorr.shape
        satAdj = np.zeros(satShape, dtype=np.float64)
    
        p = slat[:, None]
        satAdj[:, :, satCorr == 0] = satRatio[satCorr == 0]
        satAdj[:, :, satCorr == 1] = satRatio[satCorr == 1] * 0.36309 * (1.0 - 5.0 * p * p) / p
        satAdj[:, :, satCorr == 2] = satRatio[satCorr == 2] * 2.59808 * p
    
        sumc = 1.0 + np.sum(np.swapaxes(satAdj, 0, 1) * np.cos(uu * 2.0 * np.pi)[:, :], axis=-1)  # [:,:,:]
        sums = np.sum(np.swapaxes(satAdj, 0, 1) * np.sin(uu * 2.0 * np.pi)[:, :], axis=-1)
    
        M_f = np.sqrt((sumc * sumc) + (sums * sums))
        M_u = np.arctan2(sums, sumc) / (2.0 * np.pi)
    
    
    
        # Minor Constituents
        icons = newshalls['shallcons']['icon'][sc, 0]
        for value in np.arange(0,45):
            icons[icons==value]=np.where(rc==value)[0][0]
        factors = newshalls['shallcons']['factor'][sc, 0]
        m_f = np.power(M_f[:, :, icons], np.abs(factors))
        m_f[m_f == 0] = 1.0
        m_f = np.prod(m_f, axis=-1)
        m_u = np.sum(M_u[:, :, icons] * factors, axis=-1)
        m_v = np.sum(M_v[:, icons] * factors, axis=-1)
        m_freq = np.sum(M_freq[:, icons] * factors, axis=-1)
    
    
    
    
        M_revgmt = M_freq * dthr[:, None] + M_u + M_v - M_phi / 360.0
        M_radgmt = 2.0 * np.pi * np.modf(M_revgmt)[0]
        M_res = np.sum(M_f * M_A * np.cos(M_radgmt), axis=-1)
        M_res2 = np.sum(M_f * M_A * np.sin(M_radgmt), axis=-1)
    
        m_revgmt = m_freq * dthr[:, None] + m_u + m_v - m_phi / 360.0
        m_radgmt = 2.0 * np.pi * np.modf(m_revgmt)[0]
        m_res = np.sum(m_f * m_A * np.cos(m_radgmt), axis=-1)
        m_res2 = np.sum(m_f * m_A * np.sin(m_radgmt), axis=-1)
    
        ress = M_res + m_res -M_A[0]
        res2s = M_res2 + m_res2 - M_A[0]