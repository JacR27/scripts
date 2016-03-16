#import time
import datetime
from pprint import pprint as pp
import re
import yaml
import numpy as np


def open_starling_stdout(filename):
    is_start = 1
    with open(filename) as fh:
        starttime = get_start_time(fh)
        preamble = parse_preamble(fh,starttime)
        data = parse_run(fh,starttime)
        dump_data(data)
        
        


def dump_data(data,fmat='csv'):
    if fmat == 'csv':
        with open("data.csv",'w') as oh:
            oh.write("task\tadded\tlaunched\tstarting\tcomplete\n")
            for i in sorted(list(data.keys())):
                d = data[i]
                if len(d) != 4:
                    pp([i,d])
                    continue
                index = i
                added = d[0][0]
                launched = d[1][0]
                starting = d[2][0]
                complete = d[3][0]
                oh.write("{}\t{}\t{}\t{}\t{}\n".format(index,added,launched,starting,complete))
                
    elif fmat == 'yaml':
        yaml.dump(data, out.yaml, default_flow_style=False)
    
        

#FINALRUN = "[WorkflowRunner] Starling workflow successfully completed."
FINALPREAMBLE = "[WorkflowRunner] [RunParameters] mailTo"
FINALRUN = "[WorkflowRunner] Manta workflow successfully completed."
def get_start_time(fh):
    starttime = convert_date_time(fh.readline().strip().split()[0])
    fh.seek(0)
    return starttime

def parse_preamble(fh,starttime):
    preamble = []
    for line in fh:
        line = process_line(line,starttime)
        preamble.append(line)
        if FINALPREAMBLE in line[1]:
            return preamble
        
def parse_run(fh,starttime):
    data = {}
    COMMANDADDER = '[WorkflowRunner] Adding'
    COMMANDLAUNCHER = '[TaskManager] Launching'
    COMMANDRUNNER = '[TaskRunner'
    COMMANDCOMPLETE = '[TaskManager] Completed'
    for line in fh:
        line = process_line(line,base = starttime)
        if FINALRUN in line[1]:
            return data
        if COMMANDADDER in line[1] or  COMMANDLAUNCHER in line[1]:
            task = re.sub(r"CallGenome\+callGenomeSegment_|CallGenome\+","",re.findall(r"'([^']*)'",line[1])[0],count=1)
            if data.get(task) == None:
                data[task] = [line]
            else:
                data[task].append(line)
        
        elif  COMMANDCOMPLETE in line[1]:
            task = re.sub(r"CallGenome\+callGenomeSegment_|CallGenome\+","",re.findall(r"'([^']*)'",line[1])[0],count=1)
            if data.get(task) == None:
                data[task] = [line]
            else:
                data[task].append(line)
        
        elif COMMANDRUNNER in line[1]:
            task = re.sub(r"CallGenome\+callGenomeSegment_|CallGenome\+","",line[1].split()[0][12:-1],count=1)
            if data.get(task) == None:
                data[task] = [line]
            else:
                data[task].append(line)
        
        else:
            print(line)
            
        
def process_line(line,base):
    line = line.strip().split(None, 3)
    dateTime = (convert_date_time(line[0]) - base).total_seconds()
    job = line[3]
    return [dateTime, job]

def convert_date_time(dateTimeString):
    return datetime.datetime.strptime(dateTimeString, '[%Y-%m-%dT%H:%M:%S.%f]')

if __name__ == '__main__':
    open_starling_stdout('Manta.Run1.stderr.txt')
