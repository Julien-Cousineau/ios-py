import os,sys,argparse
from xxxx import main

def parseArgs():
    command_parser  = argparse.ArgumentParser(description='Script XXX')
    
    subparsers = command_parser.add_subparsers(help='Choose a command')
    
    command1_parser = subparsers.add_parser('command1', help='"Command 1" help')
    command1_parser.add_argument('--argument1',required=True, help='argument for Command 1')
    command1_parser.set_defaults(action='command1')
    
    command2_parser = subparsers.add_parser('command2', help='"Command 2" help')
    command2_parser.add_argument('--argument1',required=True, help='argument for Command 2')
    command2_parser.set_defaults(action='command2')
    
    return command_parser.parse_args()

if __name__== "__main__":
    main(vars(parseArgs()))