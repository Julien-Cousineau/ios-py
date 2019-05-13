import os
import numpy as np

filePath = r'IOS_tidetbl.txt'

satType = np.dtype([('deld', '3i4'),
                             ('phase', 'f8'),
                             ('ratio', 'f8'),
                             ('corr', 'i4'),
                             ('adj', 'f8')])
shallconType = np.dtype([('icon', 'i4'),
                                  ('factor', 'f8')])

newconType = np.dtype([('name', '|S8'),
                                ('dood', '6f8'),
                                ('nsats', 'i4'),
                                ('sats', satType, 10),
                                ('phase', 'f8'),
                                ('f', 'f8'),
                                ('u', 'f8'),
                                ('v', 'f8'),
                                ('freq', 'f8')])

newshallType = np.dtype([('name', '|S8'),
                     ('ncon', 'i4'),
                     ('shallcons', shallconType, 5),
                     ('f', 'f8'),
                     ('u', 'f8'),
                     ('v', 'f8'),
                     ('freq', 'f8')])


def getNewshall(readline, newcons, newshalls):
  # Description here...
  
  name = readline[0]
  ncon = int(readline[1])
  
  newshall = np.zeros(1, dtype=newshallType)
  newshall['name'] = name
  newshall['ncon'] = ncon
  
  icount = 0
  for i in range(2, len(readline), 2):
    factor = readline[i]
    conname = readline[i + 1]
    if factor != "" and conname != "":
      inewcon = [newcon['name'][0] for newcon in newcons].index(conname)
      shallcon = np.zeros(1, dtype=shallconType)
      shallcon['factor'] = float(factor)
      shallcon['icon'] = inewcon
      newshall['shallcons'][0, icount] = shallcon[0]
      icount += 1
  
  newshalls.append(newshall)
  
  return newshalls

def getNewcon(readline, newcons, icount):
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
      newcon = np.zeros(1, dtype=newconType)
      newcon['name'] = name
      newcon['dood'] = [float(temp[1]), float(temp[2]), float(temp[3]),
                        float(temp[4]), float(temp[5]), float(temp[6])]
      newcon['phase'] = float(temp[7])
      newcon['nsats'] = nsats
      newcons.append(newcon)
      icount = 0
    else:
      
      for i in range(1, len(temp), 5):
        sat = np.zeros(1, dtype=satType)
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

def getIOSDatabase():
  # Description here...
  print os.path.abspath(filePath)
  folder = os.path.dirname(os.path.realpath(__file__))
  
  with open(os.path.join(folder,filePath), 'r') as file:
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
          newcons, icount = getNewcon(readline, newcons, icount)
        else:
          newshalls = getNewshall(readline, newcons, newshalls)
      else:
        con = False
    
    newcons = np.asarray(newcons)
    newshalls = np.asarray(newshalls)
    
    return newcons[:,0], newshalls[:,0]