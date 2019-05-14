import os,sys,argparse
from ios import IOS
import numpy as np
import pandas as pd
from dtypes import stationType

def parseArgs():
    command_parser  = argparse.ArgumentParser(description='IOS')
    
    subparsers = command_parser.add_subparsers(help='Choose a command')
    
    decompose_parser = subparsers.add_parser('decompose', help='"decompose" help')
    decompose_parser.add_argument('--folder',required=True, help='Folder path')
    decompose_parser.add_argument('--timeSeriesPath',required=True,help='Timeseries path(.csv)')
    decompose_parser.add_argument('--datetimeColumn',required=True,help='Datetime column')
    decompose_parser.add_argument('--stationName',required=True,help='Station Name')
    decompose_parser.add_argument('--stationID',type=int,required=True,help='Station ID')
    decompose_parser.add_argument('--longitude',type=float,required=True,help='Station longitude')
    decompose_parser.add_argument('--latitude',type=float,required=True,help='Station latitude')
    decompose_parser.add_argument('--constituents',default=None, help='Select constituents')
    decompose_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.csv)')
    decompose_parser.set_defaults(action='decompose')
    
    decomposeA_parser = subparsers.add_parser('decomposeArray', help='"decomposeArray" help')
    decomposeA_parser.add_argument('--folder',required=True, help='Folder path')
    decomposeA_parser.add_argument('--timeSeriesPath',required=True,help='Timeseries path(.csv)')
    decomposeA_parser.add_argument('--stationFilePath',required=True,help='Station ID')
    decomposeA_parser.add_argument('--constituents',default=None, help='Select constituents')
    decomposeA_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.csv)')
    decomposeA_parser.set_defaults(action='decomposeArray')
    
    decomposeSLF_parser = subparsers.add_parser('decomposeSLF', help='"decomposeSLF" help')
    decomposeSLF_parser.add_argument('--folder',required=True, help='Folder path')
    decomposeSLF_parser.add_argument('--slfPath',required=True,help='SLF path(.slf)')
    decomposeSLF_parser.add_argument('--constituents',default=None, help='Select constituents')
    decomposeSLF_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.npy)')
    decomposeSLF_parser.set_defaults(action='decomposeSLF')    
    
    compose_parser = subparsers.add_parser('compose', help='"compose" help')
    compose_parser.add_argument('--folder',required=True, help='Folder path')
    compose_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.csv or .npy)')
    compose_parser.add_argument('--constituents',default=None, help='Select constituents')
    compose_parser.add_argument('--stations',default=None, help='Station IDs')
    compose_parser.add_argument('--startDate',required=True, help='Start Date')
    compose_parser.add_argument('--endDate',required=True, help='End Date')
    compose_parser.add_argument('--step',type=int,default=10, help='Step in minutes')
    compose_parser.add_argument('--etaOnly',type=bool,default=True, help='Eta only (exclude u and v)')
    compose_parser.add_argument('--outPath',required=True, help='Output timeseries Path (.csv, .slf)')
    compose_parser.set_defaults(action='compose')
    
    return command_parser.parse_args()

def main(props):
    if props['action']=='decompose':return decompose(props)
    if props['action']=='decomposeArray':return decomposeArray(props)
    if props['action']=='decomposeSLF':return decomposeSLF(props)
    if props['action']=='compose':return compose(props)

def decompose(props):
    folder = props['folder']
    timeSeriesPath = os.path.join(folder,props['timeSeriesPath'])
    stationName = props['stationName']
    stationID = props['stationID']
    longitude = props['longitude']
    latitude = props['latitude']
    constituents = props['constituents'].split(",")
    conFilePath = os.path.join(folder,props['conFilePath'])
    npyPath = os.path.splitext(conFilePath)[0] + ".npy"
    
    stations = np.zeros(1, dtype=stationType)
    stations['name'][0] = stationName
    stations['id'][0]   = stationID
    stations['xy'][0]   = (longitude, latitude)
    stations['proj'][0] ='epsg:3857'
    
    ios = IOS(constituents,stations,npy=npyPath)
    
    ts = pd.read_csv(timeSeriesPath,parse_dates=['Datetime'])
    if not 'Datetime' in ts.columns:sys.exit("Datetime column does not exist")
    if not 'eta' in ts.columns:sys.exit("eta column does not exist")
    
    datetime = ts['Datetime'].values
    eta      = ts['eta'].values
    u        = ts['u'].values if 'u' in ts.columns else np.zeros(len(eta))
    v        = ts['v'].values if 'v' in ts.columns else np.zeros(len(eta))
    
    ios.extractConstituents(datetime,np.asarray([[eta,u,v]]))
    ios.to_ccsv(conFilePath)    

def decomposeArray(props):
    folder = props['folder']
    timeSeriesPath = os.path.join(folder,props['timeSeriesPath'])
    stationFilePath = os.path.join(folder,props['stationFilePath'])
    conFilePath = os.path.join(folder,props['conFilePath'])
    npyPath = os.path.splitext(conFilePath)[0] + ".npy"
    
    constituents = props['constituents'].split(",")
    
    df = pd.read_csv(stationFilePath)
    
    stations = np.zeros(len(df), dtype=stationType)
    stations['name'] = df['StationName'].values
    stations['id']   = df['StationName'].values
    stations['xy']   = np.columnstack(df['Longitude'].values,  df['Latitude'].values)
    stations['proj'][:] ='epsg:3857'
    
    ios = IOS(constituents,stations,npy=npyPath)
    
    ts = pd.read_csv(timeSeriesPath,parse_dates=['Datetime'])
    if not 'Datetime' in ts.columns:sys.exit("Datetime column does not exist")
    if not 'eta' in ts.columns:sys.exit("eta column does not exist")
    
    datetime = ts['Datetime'].values
    eta      = ts['eta'].values
    u        = ts['u'].values if 'u' in ts.columns else np.zeros(len(eta))
    v        = ts['v'].values if 'v' in ts.columns else np.zeros(len(eta))
    
    ios.extractConstituents(datetime,np.asarray([[eta,u,v]]))
    ios.to_ccsv(conFilePath)

def decomposeSLF(props):
    return

def compose(props):
    return

def generate(props):
    modelPath = os.path.join(props['folder'],props['modelPath'])
    outPath = os.path.join(props['folder'],props['outPath'])
    start = props['startDate']
    end = props['endDate']
    step = props['step']
    type = props['type']
    
    datetimes = np.arange(start, end, np.timedelta64(step, 'm'), dtype='datetime64')
    
    return None

def extract(props):
    
    return None

if __name__== "__main__":
    main(vars(parseArgs()))