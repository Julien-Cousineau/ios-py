import os,sys,argparse
from ios import IOS
import numpy as np
import pandas as pd
from dtypes import stationType,createStation
from readCHS import read_csv,read_npy
import filecmp

def parseArgs():
    command_parser  = argparse.ArgumentParser(description='IOS')
    
    subparsers = command_parser.add_subparsers(help='Choose a command')
    
    decompose_parser = subparsers.add_parser('decompose', help='"decompose" help')
    decompose_parser.add_argument('--folder',required=True, help='Folder path')
    decompose_parser.add_argument('--timeSeriesPath',required=True,help='Timeseries path(.csv)')
    decompose_parser.add_argument('--stationName',required=True,help='Station Name')
    decompose_parser.add_argument('--stationID',type=int,required=True,help='Station ID')
    decompose_parser.add_argument('--longitude',type=float,required=True,help='Station longitude')
    decompose_parser.add_argument('--latitude',type=float,required=True,help='Station latitude')
    decompose_parser.add_argument('--constituents',default=None, help='Select constituents')
    decompose_parser.add_argument('--consRemove',default=None, help='Remove constituents')
    decompose_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.csv)')
    decompose_parser.set_defaults(action='decompose')
    
    decomposeA_parser = subparsers.add_parser('decomposeArray', help='"decomposeArray" help')
    decomposeA_parser.add_argument('--folder',required=True, help='Folder path')
    decomposeA_parser.add_argument('--timeSeriesPath',required=True,help='Timeseries path(.csv)')
    decomposeA_parser.add_argument('--stationFilePath',required=True,help='Station ID')
    decomposeA_parser.add_argument('--constituents',default=None, help='Select constituents')
    decomposeA_parser.add_argument('--consRemove',default=None, help='Remove constituents')
    decomposeA_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.csv)')
    decomposeA_parser.set_defaults(action='decomposeArray')
    
    decomposeSLF_parser = subparsers.add_parser('decomposeSLF', help='"decomposeSLF" help')
    decomposeSLF_parser.add_argument('--folder',required=True, help='Folder path')
    decomposeSLF_parser.add_argument('--slfPath',required=True,help='SLF path(.slf)')
    decomposeSLF_parser.add_argument('--constituents',default=None, help='Select constituents')
    decomposeSLF_parser.add_argument('--consRemove',default=None, help='Remove constituents')
    decomposeSLF_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.npy)')
    decomposeSLF_parser.set_defaults(action='decomposeSLF')    
    
    compose_parser = subparsers.add_parser('compose', help='"compose" help')
    compose_parser.add_argument('--folder',required=True, help='Folder path')
    compose_parser.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.csv or .npy)')
    compose_parser.add_argument('--constituents',default=None, help='Select constituents')
    compose_parser.add_argument('--consRemove',default=None, help='Remove constituents')
    compose_parser.add_argument('--stationIDS',required=True, help='Station IDs')
    compose_parser.add_argument('--startDate',required=True, help='Start Date')
    compose_parser.add_argument('--endDate',required=True, help='End Date')
    compose_parser.add_argument('--step',type=int,default=10, help='Step in minutes')
    compose_parser.add_argument('--vars',default=None, help='eta,u,v')
    compose_parser.add_argument('--outPath',required=True, help='Output timeseries Path (.csv, .slf)')
    compose_parser.set_defaults(action='compose')

    validate_parser = subparsers.add_parser('validate', help='"validate" help')
    validate_parser.add_argument('--folder',required=True, help='Folder path')
    validate_parser.add_argument('--input',required=True, help='Input path')
    validate_parser.add_argument('--output',required=True, help='Output path')
    validate_parser.set_defaults(action='validate')
    
    getTopConstants = subparsers.add_parser('getTopConstants', help='"getTopConstants" help')
    getTopConstants.add_argument('--folder',required=True, help='Folder path')
    getTopConstants.add_argument('--conFilePath',required=True, help='Harmonic constants filepath(.csv or .npy)')
    getTopConstants.add_argument('--stationIDS',required=True, help='Station IDS')
    getTopConstants.add_argument('--nTop',type=int,required=True, help='N top constants')
    getTopConstants.set_defaults(action='getTopConstants')    
    
    return command_parser.parse_args()

def main(props):
    if props['action']=='decompose':return decompose(props)
    if props['action']=='decomposeArray':return decomposeArray(props)
    if props['action']=='decomposeSLF':return decomposeSLF(props)
    if props['action']=='compose':return compose(props)
    if props['action']=='validate':return validate(props)
    if props['action']=='getTopConstants':return getTopConstants(props)

def getTopConstants(props):
    nTop = props['nTop']
    stations = getStations(props)

    index =np.argsort(-stations['constituents']['eta'][:,:,0])[:,:nTop]
    
    cons =  stations['constituents']['name'][0,index]
    
    for i in range(len(stations)):
        print stations[i]['id']
        print cons[i]

def getDatetimes(props):
    start = props['startDate']
    end = props['endDate']
    step = props['step']
    return np.arange(start, end, np.timedelta64(step, 'm'), dtype='datetime64')

def getStations(props):
    folder = props['folder']
    stationIDS = np.asarray(props['stationIDS'].split(","),dtype=np.int32)
    conFilePath = os.path.join(folder,props['conFilePath'])
    ext =  os.path.splitext(conFilePath)[1]
    if ext != '.csv' and ext != '.npy':sys.exit("Con filepath should be .csv or .npy")
    if ext=='.csv':stations = read_csv(conFilePath)
    if ext=='.npy':stations = read_npy(conFilePath)
   
    stations = stations[np.isin(stations['id'],stationIDS)]
    if len(stations)==0:sys.exit("Stations ids does not exist")
    return stations

def getCons(props):
    constituents = props['constituents']
    consRemove = props['consRemove']
    
    if constituents is None:
        stations = getStations(props)
        cons = stations['constituents']['name'][0, :]
    else:
        cons = np.asarray(constituents.split(","))
    
    cons=cons[cons!='']
    cons = cons if consRemove is None else cons[np.isin(cons, consRemove.split(','),invert=True)]
    return cons

def getConsStations(props):
    return getCons(props),getStations(props)

def compose(props):
    outPath = os.path.join(props['folder'],props['outPath'])
    
    datetimes = getDatetimes(props)
    cons,stations = getConsStations(props)
    
    vars = None if props['vars'] is None else props['vars'].split(',')
    ios = IOS(cons,stations,debug=True)
    ios.getTS(datetimes)
    ios.to_csv(outPath,vars=vars, index=None, header=True,date_format="%Y/%m/%d %H:%M:%S")


def validate(props):
    folder = props['folder']
    input = os.path.join(folder,props['input'])
    output = os.path.join(folder,props['output'])
    print filecmp.cmp(input,output)

def decompose(props):
    folder = props['folder']
    timeSeriesPath = os.path.join(folder,props['timeSeriesPath'])
    stationName = props['stationName']
    stationID = props['stationID']
    longitude = props['longitude']
    latitude = props['latitude']
    
    conFilePath = os.path.join(folder,props['conFilePath'])
    fileType = '.csv' if os.path.splitext(conFilePath)[1]=='.csv' else '.npy'
    
    stations = createStation(stationName,stationID,(longitude, latitude),'epsg:3857')
    
    cons = getCons(props)
    ios = IOS(cons,stations,npy='temp.npy')
    
    ts = pd.read_csv(timeSeriesPath,parse_dates=['Datetime'])
    if not 'Datetime' in ts.columns:sys.exit("Datetime column does not exist")
    if not 'eta' in ts.columns:sys.exit("eta column does not exist")
    
    datetime = ts['Datetime'].values
    eta      = ts['eta'].values
    u        = ts['u'].values if 'u' in ts.columns else np.zeros(len(eta))
    v        = ts['v'].values if 'v' in ts.columns else np.zeros(len(eta))
   
    ios.extractConstituents(datetime,np.asarray([[eta,u,v]]))
    if fileType=='.csv': ios.stations_to_csv(conFilePath)
    else:ios.stations_to_npy(conFilePath)

def decomposeArray(props):
    folder = props['folder']
    timeSeriesPath = os.path.join(folder,props['timeSeriesPath'])
    stationFilePath = os.path.join(folder,props['stationFilePath'])
    constituents = props['constituents'].split(",")
    conFilePath = os.path.join(folder,props['conFilePath'])
    fileType = '.csv' if os.path.splitext(conFilePath)[1]=='.csv' else '.npy'
    
    df = pd.read_csv(stationFilePath)
    
    stations = np.zeros(len(df), dtype=stationType)
    stations['name'] = df['StationName'].values
    stations['id']   = df['StationName'].values
    stations['xy']   = np.columnstack(df['Longitude'].values,  df['Latitude'].values)
    stations['proj'][:] ='epsg:3857'
    
    ios = IOS(constituents,stations,npy='temp.npy')
    
    ts = pd.read_csv(timeSeriesPath,parse_dates=['Datetime'])
    if not 'Datetime' in ts.columns:sys.exit("Datetime column does not exist")
    if not 'eta' in ts.columns:sys.exit("eta column does not exist")
    
    datetime = ts['Datetime'].values
    eta      = ts['eta'].values
    u        = ts['u'].values if 'u' in ts.columns else np.zeros(len(eta))
    v        = ts['v'].values if 'v' in ts.columns else np.zeros(len(eta))
    
    ios.extractConstituents(datetime,np.asarray([[eta,u,v]]))
    if fileType=='.csv': ios.stations_to_csv(conFilePath)
    else:ios.stations_to_npy(conFilePath)

def decomposeSLF(props):
    return




if __name__== "__main__":
    main(vars(parseArgs()))