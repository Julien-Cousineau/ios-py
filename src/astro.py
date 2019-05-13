import numpy as np

astroType = np.dtype([('d1', 'f8'),
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


def getAstro(datetimes):
  astro = np.zeros(len(datetimes), dtype=astroType)
  
  kb = np.datetime64('1899-12-31T00:00')
  kh = (datetimes - kb).astype('timedelta64[D]').astype('f8')
  d1 = kh - 0.5
  
  astro['d1'] = d1
  astro['d12'] = d1 * d1
  astro['d2'] = d1 * 1.0E-04
  astro['d22'] = astro['d2'] * astro['d2']
  astro['d23'] = astro['d2'] ** 3.0
  astro['f'] = 360.0
  astro['f2'] = astro['f'] / 365.0
  astro['h'] = np.modf((2.79696678E+02 + d1 * 9.856473354E-01 + astro['d22'] * 2.267E-05) / astro['f'])[0]
  astro['pp'] = np.modf((2.81220844E+02 + d1 * 4.70684E-05 + astro['d22'] * 3.39E-05 + astro['d23'] * 7.0E-08) / astro['f'])[0]
  astro['s'] = np.modf((2.70434164E+02 + d1 * 1.31763965268E+01 - astro['d22'] * 8.5E-05 + astro['d23'] * 3.9E-08) / astro['f'])[0]
  astro['p'] = np.modf((3.34329556E+02 + d1 * 1.114040803E-01 - astro['d22'] * 7.739E-04 - astro['d23'] * 2.6E-07) / astro['f'])[0]
  astro['np'] = np.modf((-2.59183275E+02 + d1 * 5.29539222E-02 - astro['d22'] * 1.557E-04 - astro['d23'] * 5.0E-08) / astro['f'])[0]
  astro['dh'] = (9.856473354E-01 + d1 * 2.267E-05 * 2.0E-08) / astro['f2']
  astro['dpp'] = (4.70684E-05 + d1 * 3.39E-05 * 2.0E-08 + astro['d12'] * 7.0E-08 * 3.0E-12) / astro['f2']
  astro['ds'] = (1.31763965268E+01 - d1 * 8.5E-05 * 2.0E-08 + astro['d12'] * 3.9E-08 * 3.0E-12) / astro['f2']
  astro['dp'] = (1.114040803E-01 - d1 * 7.739E-04 * 2.0E-08 - astro['d12'] * 2.6E-07 * 3.0E-12) / astro['f2']
  astro['dnp'] = (5.29539222E-02 - d1 * 1.557E-04 * 2.0E-08 - astro['d12'] * 5.0E-08 * 3.0E-12) / astro['f2']
  
  ktmp = kh * 24.0
  astro['hh'] = ktmp - (np.floor(ktmp / 24.0) * 24.0)
  astro['tau'] = astro['hh'] / 24.0 + astro['h'] - astro['s']
  astro['dtau'] = 365.0 + astro['dh'] - astro['ds']
  
  return astro