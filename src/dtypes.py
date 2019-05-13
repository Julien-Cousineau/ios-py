import numpy as np

constituentType = np.dtype([('name', '|S8'),
                            ('eta', '2f8'),
                            ('u', '2f8'),
                            ('v', '2f8'),
                            ('r', 'f8'),
                            ('zeta', 'f8'),
                            ('konans', 'i4')])

stationType = np.dtype([('name', '|S32'),
                        ('id', 'i4'),
                        ('xy', '2f8'),
                        ('proj', '|S16'),
                        ('constituents', constituentType, 200)])