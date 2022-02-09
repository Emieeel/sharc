#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:24:40 2022

@author: emielkoridon
"""

import os
import shutil
# External Calls to MOLCAS
import subprocess as sp
# Command line arguments
import sys
# Regular expressions
import re
# debug print for dicts and arrays
import pprint
# sqrt and other math
import math
import cmath
# runtime measurement
import datetime
# copy of arrays of arrays
from copy import deepcopy
# parallel calculations
from multiprocessing import Pool
import time
from socket import gethostname
import itertools
# write debug traces when in pool threads
import traceback

# =========================================================0
# compatibility stuff

if sys.version_info[0]!=3:
    print('This is a script for Python 3!')
    sys.exit(0)

version='2.1'
versiondate=datetime.date(2019,9,1)



changelogstring='''
02.07.2022:
- Writing SA-OO-VQE / SHARC interface...
'''

# ======================================================================= #
# holds the system time when the script was started
starttime=datetime.datetime.now()

# global variables for printing (PRINT gives formatted output, DEBUG gives raw output)
DEBUG=False
PRINT=True

# conversion factors
au2a=0.529177211
rcm_to_Eh=4.556335e-6



# =============================================================================================== #
# =============================================================================================== #
# =========================================== general routines ================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def readfile(filename):
  try:
    f=open(filename)
    out=f.readlines()
    f.close()
  except IOError:
    print('File %s does not exist!' % (filename))
    sys.exit(12)
  return out

# ======================================================================= #
def writefile(filename,content):
  # content can be either a string or a list of strings
  try:
    f=open(filename,'w')
    if isinstance(content,list):
      for line in content:
        f.write(line)
    elif isinstance(content,str):
      f.write(content)
    else:
      print('Content %s cannot be written to file!' % (content))
    f.close()
  except IOError:
    print('Could not write to file %s!' % (filename))
    sys.exit(13)


# ======================================================================= #
def isbinary(path):
  return (re.search(r':.* text',sp.Popen(["file", '-L', path], stdout=sp.PIPE).stdout.read())is None)

# ======================================================================= #
def eformat(f, prec, exp_digits):
    '''Formats a float f into scientific notation with prec number of decimals and exp_digits number of exponent digits.

    String looks like:
    [ -][0-9]\.[0-9]*E[+-][0-9]*

    Arguments:
    1 float: Number to format
    2 integer: Number of decimals
    3 integer: Number of exponent digits

    Returns:
    1 string: formatted number'''

    s = "% .*e"%(prec, f)
    mantissa, exp = s.split('e')
    return "%sE%+0*d"%(mantissa, exp_digits+1, int(exp))

# ======================================================================= #
def measuretime():
    '''Calculates the time difference between global variable starttime and the time of the call of measuretime.

    Prints the Runtime, if PRINT or DEBUG are enabled.

    Arguments:
    none

    Returns:
    1 float: runtime in seconds'''

    endtime=datetime.datetime.now()
    runtime=endtime-starttime
    if PRINT or DEBUG:
        hours=runtime.seconds/3600
        minutes=runtime.seconds/60-hours*60
        seconds=runtime.seconds%60
        print('==> Runtime:\n%i Days\t%i Hours\t%i Minutes\t%i Seconds\n\n' % (runtime.days,hours,minutes,seconds))
    total_seconds=runtime.days*24*3600+runtime.seconds+runtime.microseconds/1.e6
    return total_seconds

# ======================================================================= #
def removekey(d,key):
    '''Removes an entry from a dictionary and returns the dictionary.

    Arguments:
    1 dictionary
    2 anything which can be a dictionary keyword

    Returns:
    1 dictionary'''

    if key in d:
        r = dict(d)
        del r[key]
        return r
    return d

# ======================================================================= #         OK
def containsstring(string,line):
    '''Takes a string (regular expression) and another string. Returns True if the first string is contained in the second string.

    Arguments:
    1 string: Look for this string
    2 string: within this string

    Returns:
    1 boolean'''

    a=re.search(string,line)
    if a:
        return True
    else:
        return False


# =============================================================================================== #
# =============================================================================================== #
# ============================= iterator routines  ============================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def itmult(states):

    for i in range(len(states)):
        if states[i]<1:
            continue
        yield i+1
    return

# ======================================================================= #
def itnmstates(states):

    for i in range(len(states)):
        if states[i]<1:
            continue
        for k in range(i+1):
            for j in range(states[i]):
                yield i+1,j+1,k-i/2.
    return


# =============================================================================================== #
# =============================================================================================== #
# =========================================== output extraction ================================= #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def getversion(out,MOLCAS):
    # allowedrange=[ (18.0,18.999), (8.29999,8.30001) ]
    # # first try to find $MOLCAS/.molcasversion
    # molcasversion=os.path.join(MOLCAS,'.molcasversion')
    # if os.path.isfile(molcasversion):
    #     vf=open(molcasversion)
    #     string=vf.readline()
    #     vf.close()
    # # otherwise try to read this from the output file
    # else:
    #     string=''
    #     for i in range(50):
    #         line=out[i]
    #         s=line.split()
    #         for j,el in enumerate(s):
    #             if 'version' in el:
    #                 string=s[j+1]
    #                 break
    #         if string!='':
    #             break
    # a=re.search('[0-9]+\.[0-9]+',string)
    # if a==None:
    #     print('No MOLCAS version found.\nCheck whether MOLCAS path is set correctly in MOLCAS.resources\nand whether $MOLCAS/.molcasversion exists.')
    #     sys.exit(17)
    # v=float(a.group())
    # if not any( [ i[0]<=v<=i[1] for i in allowedrange ] ):
    #     #allowedrange[0]<=v<=allowedrange[1]:
    #     print('MOLCAS version %3.1f not supported! ' % (v))
    #     sys.exit(18)
    # if DEBUG:
    #     print('Found MOLCAS version %3.1f\n' % (v))
    return 0.01

# =============================================================================================== #
# =============================================================================================== #
# =========================================== print routines ==================================== #
# =============================================================================================== #
# =============================================================================================== #

# ======================================================================= #
def printheader():
    '''Prints the formatted header of the log file. Prints version number and version date

    Takes nothing, returns nothing.'''

    print(starttime,gethostname(),os.getcwd())
    if not PRINT:
        return
    string='\n'
    string+='  '+'='*80+'\n'
    string+='||'+' '*80+'||\n'
    string+='||'+' '*27+'SHARC - SA-OO-VQE - Interface'+' '*27+'||\n'
    string+='||'+' '*80+'||\n'
    string+='||'+' '*19+'Authors: Emiel Koridon'+' '*20+'||\n'
    string+='||'+' '*80+'||\n'
    string+='||'+' '*int(36-(len(version)+1)//2)+'Version: %s' % (version)+' '*int(35-(len(version))//2)+'||\n'
    lens=len(versiondate.strftime("%d.%m.%y"))
    string+='||'+' '*int(37-lens//2)+'Date: %s' % (versiondate.strftime("%d.%m.%y"))+' '*int(37-(lens+1)//2)+'||\n'
    string+='||'+' '*80+'||\n'
    string+='  '+'='*80+'\n\n'
    print(string)
    if DEBUG:
        print(changelogstring)


# =============================================================================================== #
# =============================================================================================== #
# =========================================== SUBROUTINES TO readQMin =========================== #
# =============================================================================================== #
# =============================================================================================== #


# ======================================================================= #
def checkscratch(SCRATCHDIR):
    '''Checks whether SCRATCHDIR is a file or directory. If a file, it quits with exit code 1, if its a directory, it passes. If SCRATCHDIR does not exist, tries to create it.

    Arguments:
    1 string: path to SCRATCHDIR'''

    exist=os.path.exists(SCRATCHDIR)
    if exist:
        isfile=os.path.isfile(SCRATCHDIR)
        if isfile:
            print('$SCRATCHDIR=%s exists and is a file!' % (SCRATCHDIR))
            sys.exit(28)
    else:
        try:
            os.makedirs(SCRATCHDIR)
        except OSError:
            print('Can not create SCRATCHDIR=%s\n' % (SCRATCHDIR))
            sys.exit(29)

# ======================================================================= #
def removequotes(string):
  if string.startswith("'") and string.endswith("'"):
    return string[1:-1]
  elif string.startswith('"') and string.endswith('"'):
    return string[1:-1]
  else:
    return string

# ======================================================================= #
def getsh2caskey(sh2cas,key):
  i=-1
  while True:
    i+=1
    try:
      line=re.sub('#.*$','',sh2cas[i])
    except IndexError:
      break
    line=line.split(None,1)
    if line==[]:
      continue
    if key.lower() in line[0].lower():
      return line
  return ['','']

# ======================================================================= #
def get_sh2cas_environ(sh2cas,key,environ=True,crucial=True):
  line=getsh2caskey(sh2cas,key)
  if line[0]:
    LINE=line[1]
    LINE=removequotes(LINE).strip()
  else:
    if environ:
      LINE=os.getenv(key.upper())
      if not LINE:
        if crucial:
          print('Either set $%s or give path to %s in MOLCAS.resources!' % (key.upper(),key.upper()))
          sys.exit(30)
        else:
          return ''
    else:
      if crucial:
        print('Give path to %s in MOLCAS.resources!' % (key.upper()))
        sys.exit(31)
      else:
        return ''
  LINE=os.path.expandvars(LINE)
  LINE=os.path.expanduser(LINE)
  if containsstring(';',LINE):
    print("$%s contains a semicolon. Do you probably want to execute another command after %s? I can't do that for you..." % (key.upper(),key.upper()))
    sys.exit(32)
  return LINE

# ======================================================================= #
def get_pairs(QMinlines,i):
  nacpairs=[]
  while True:
    i+=1
    try:
      line=QMinlines[i].lower()
    except IndexError:
      print('"keyword select" has to be completed with an "end" on another line!')
      sys.exit(33)
    if 'end' in line:
      break
    fields=line.split()
    try:
      nacpairs.append([int(fields[0]),int(fields[1])])
    except ValueError:
      print('"nacdr select" is followed by pairs of state indices, each pair on a new line!')
      sys.exit(34)
  return nacpairs,i

# ======================================================================= #         OK
def readQMin(QMinfilename):
    '''Reads the time-step dependent information from QMinfilename. This file contains all information from the current SHARC job: geometry, velocity, number of states, requested quantities along with additional information. The routine also checks this input and obtains a number of environment variables necessary to run MOLCAS.

    Steps are:
    - open and read QMinfilename
    - Obtain natom, comment, geometry (, velocity)
    - parse remaining keywords from QMinfile
    - check keywords for consistency, calculate nstates, nmstates
    - obtain environment variables for path to MOLCAS and scratch directory, and for error handling

    Arguments:
    1 string: name of the QMin file

    Returns:
    1 dictionary: QMin'''

    # read QMinfile
    QMinlines=readfile(QMinfilename)
    QMin={}

    # Get natom
    try:
        natom=int(QMinlines[0])
    except ValueError:
        print('first line must contain the number of atoms!')
        sys.exit(35)
    QMin['natom']=natom
    if len(QMinlines)<natom+4:
        print('Input file must contain at least:\nnatom\ncomment\ngeometry\nkeyword "states"\nat least one task')
        sys.exit(36)

    # Save Comment line
    QMin['comment']=QMinlines[1]

    # Get geometry and possibly velocity (for backup-analytical non-adiabatic couplings)
    QMin['geo']=[]
    QMin['veloc']=[]
    hasveloc=True
    for i in range(2,natom+2):
        if not containsstring('[a-zA-Z][a-zA-Z]?[0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*.*[-]?[0-9]+[.][0-9]*', QMinlines[i]):
            print('Input file does not comply to xyz file format! Maybe natom is just wrong.')
            sys.exit(37)
        fields=QMinlines[i].split()
        for j in range(1,4):
            fields[j]=float(fields[j])
        QMin['geo'].append(fields[0:4])
        if len(fields)>=7:
            for j in range(4,7):
                fields[j]=float(fields[j])
            QMin['veloc'].append(fields[4:7])
        else:
            hasveloc=False
    if not hasveloc:
        QMin=removekey(QMin,'veloc')


    # Parse remaining file
    i=natom+1
    while i+1<len(QMinlines):
        i+=1
        line=QMinlines[i]
        line=re.sub('#.*$','',line)
        if len(line.split())==0:
            continue
        key=line.lower().split()[0]
        if 'savedir' in key:
            args=line.split()[1:]
        else:
            args=line.lower().split()[1:]
        if key in QMin:
            print('Repeated keyword %s in line %i in input file! Check your input!' % (key,i+1))
            continue  # only first instance of key in QM.in takes effect
        if len(args)>=1 and 'select' in args[0]:
            pairs,i=get_pairs(QMinlines,i)
            QMin[key]=pairs
        else:
            QMin[key]=args

    if 'unit' in QMin:
        if QMin['unit'][0]=='angstrom':
            factor=1./au2a
        elif QMin['unit'][0]=='bohr':
            factor=1.
        else:
            print('Dont know input unit %s!' % (QMin['unit'][0]))
            sys.exit(38)
    else:
        factor=1./au2a

    for iatom in range(len(QMin['geo'])):
        for ixyz in range(3):
            QMin['geo'][iatom][ixyz+1]*=factor


    if not 'states' in QMin:
        print('Keyword "states" not given!')
        sys.exit(39)
    # Calculate states, nstates, nmstates
    for i in range(len(QMin['states'])):
        QMin['states'][i]=int(QMin['states'][i])
    reduc=0
    for i in reversed(QMin['states']):
        if i==0:
            reduc+=1
        else:
            break
    for i in range(reduc):
        del QMin['states'][-1]
    nstates=0
    nmstates=0
    for i in range(len(QMin['states'])):
        nstates+=QMin['states'][i]
        nmstates+=QMin['states'][i]*(i+1)
    QMin['nstates']=nstates
    QMin['nmstates']=nmstates


    # Various logical checks
    if not 'states' in QMin:
        print('Number of states not given in QM input file %s!' % (QMinfilename))
        sys.exit(40)

    possibletasks=['h','soc','dm','grad','overlap','dmdr','socdr','ion','phases']
    if not any([i in QMin for i in possibletasks]):
        print('No tasks found! Tasks are "h", "soc", "dm", "grad","dmdr", "socdr", "overlap" and "ion".')
        sys.exit(41)

    if 'samestep' in QMin and 'init' in QMin:
        print('"Init" and "Samestep" cannot be both present in QM.in!')
        sys.exit(42)

    if 'phases' in QMin:
        QMin['overlap']=[]

    if 'overlap' in QMin:# and 'init' in QMin:
        # print('"overlap" and "phases" cannot be calculated in the first timestep! Delete either "overlap" or "init"')
        print("Overlap not implemented for SAOOVQE")
        sys.exit(43)

    if not 'init' in QMin and not 'samestep' in QMin:
        QMin['newstep']=[]

    if not any([i in QMin for i in ['h','soc','dm','grad']]) and 'overlap' in QMin:
        QMin['h']=[]

    if len(QMin['states'])>8:
        print('Higher multiplicities than octets are not supported!')
        sys.exit(44)

    if 'h' in QMin and 'soc' in QMin:
        QMin=removekey(QMin,'h')

    if 'nacdt' in QMin:
        print('Within the SHARC-MOLCAS interface, "nacdt" is not supported.')
        sys.exit(45)

    if 'molden' in QMin:
        os.environ['MOLCAS_MOLDEN']='ON'
        if 'samestep' in QMin:
            print('HINT: Not producing Molden files in "samestep" mode!')
            del QMin['molden']

    #if 'ion' in QMin:
        #print 'Ionization probabilities not implemented!'
        #sys.exit(46)

    # Check for correct gradient list
    if 'grad' in QMin:
        if len(QMin['grad'])==0 or QMin['grad'][0]=='all':
            QMin['grad']=[ i+1 for i in range(nmstates)]
            #pass
        else:
            for i in range(len(QMin['grad'])):
                try:
                    QMin['grad'][i]=int(QMin['grad'][i])
                except ValueError:
                    print('Arguments to keyword "grad" must be "all" or a list of integers!')
                    sys.exit(47)
                if QMin['grad'][i]>nmstates:
                    print('State for requested gradient does not correspond to any state in QM input file state list!')
                    sys.exit(48)

    # Process the overlap requests
    # identically to the nac requests
    if 'overlap' in QMin:
        if len(QMin['overlap'])>=1:
            nacpairs=QMin['overlap']
            for i in range(len(nacpairs)):
                if nacpairs[i][0]>nmstates or nacpairs[i][1]>nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(49)
        else:
            QMin['overlap']=[ [j+1,i+1] for i in range(nmstates) for j in range(i+1)]

    # Process the non-adiabatic coupling requests
    # type conversion has already been done
    if 'nacdr' in QMin:
        if len(QMin['nacdr'])>=1:
            nacpairs=QMin['nacdr']
            for i in range(len(nacpairs)):
                if nacpairs[i][0]>nmstates or nacpairs[i][1]>nmstates:
                    print('State for requested non-adiabatic couplings does not correspond to any state in QM input file state list!')
                    sys.exit(50)
        else:
            QMin['nacdr']=[ [j+1,i+1] for i in range(nmstates) for j in range(i)]


    # obtain the statemap 
    statemap={}
    i=1
    for imult,istate,ims in itnmstates(QMin['states']):
        statemap[i]=[imult,istate,ims]
        i+=1
    QMin['statemap']=statemap

    # get the set of states for which gradients actually need to be calculated
    gradmap=set()
    if 'grad' in QMin:
        for i in QMin['grad']:
            gradmap.add( tuple(statemap[i][0:2]) )
    gradmap=list(gradmap)
    gradmap.sort()
    QMin['gradmap']=gradmap

    # get the list of statepairs for NACdr calculation
    nacmap=set()
    if 'nacdr' in QMin:
        for i in QMin['nacdr']:
            s1=statemap[i[0]][0:2]
            s2=statemap[i[1]][0:2]
            if s1[0]!=s2[0] or s1==s2:
                continue
            if s1[1]>s2[1]:
                continue
            nacmap.add(tuple(s1+s2))
    nacmap=list(nacmap)
    nacmap.sort()
    QMin['nacmap']=nacmap







    # open MOLCAS.resources
    filename='MOLCAS.resources'
    if os.path.isfile(filename):
        sh2cas=readfile(filename)
    else:
        print('Warning: No MOLCAS.resources found!')
        print('Reading resources from SH2CAS.inp')
        sh2cas=readfile('SH2CAS.inp')

    QMin['pwd']=os.getcwd()

    QMin['molcas']=get_sh2cas_environ(sh2cas,'molcas')
    os.environ['MOLCAS']=QMin['molcas']

    QMin['tinker']=get_sh2cas_environ(sh2cas,'tinker',crucial=False)
    if QMin['tinker']=='':
        QMin=removekey(QMin,'tinker')
    else:
        os.environ['TINKER']=QMin['tinker']

    if 'ion' in QMin:
        QMin['wfoverlap']=get_sh2cas_environ(sh2cas,'wfoverlap',crucial=False)
        if not QMin['wfoverlap']:
            ciopath=os.path.join(os.path.expandvars(os.path.expanduser('$SHARC')),'wfoverlap.x')
            if os.path.isfile(ciopath):
                QMin['wfoverlap']=ciopath
            else:
                print('Give path to wfoverlap.x in MOLCAS.resources!')
                sys.exit(51)
        # read ncore and ndocc from resources
        line=getsh2caskey(sh2cas,'numfrozcore')
        if line[0]:
            try:
                QMin['ncore']=max(0,int(line[1]))
            except ValueError:
                print('numfrozcore does not evaluate to numerical value!')
                sys.exit(52)
        line=getsh2caskey(sh2cas,'numocc')
        if line[0]:
            try:
                QMin['ndocc']=int(line[1])
            except ValueError:
                print('numocc does not evaluate to numerical value!')
                sys.exit(53)


    # Set up scratchdir
    line=get_sh2cas_environ(sh2cas,'scratchdir',False,False)
    if line==None:
        line=QMin['pwd']+'/SCRATCHDIR/'
    line=os.path.expandvars(line)
    line=os.path.expanduser(line)
    line=os.path.abspath(line)
    #checkscratch(line)
    QMin['scratchdir']=line


    # Set up savedir
    if 'savedir' in QMin:
        # savedir may be read from QM.in file
        line=QMin['savedir'][0]
    else:
        line=get_sh2cas_environ(sh2cas,'savedir',False,False)
        if line==None or line=='':
            line=os.path.join(QMin['pwd'],'SAVEDIR')
    line=os.path.expandvars(line)
    line=os.path.expanduser(line)
    line=os.path.abspath(line)
    if 'init' in QMin:
        checkscratch(line)
    QMin['savedir']=line


    line=getsh2caskey(sh2cas,'debug')
    if line[0]:
        if len(line)<=1 or 'true' in line[1].lower():
            global DEBUG
            DEBUG=True

    line=getsh2caskey(sh2cas,'no_print')
    if line[0]:
        if len(line)<=1 or 'true' in line[1].lower():
            global PRINT
            PRINT=False

    QMin['memory']=500
    line=getsh2caskey(sh2cas,'memory')
    if line[0]:
        try:
            QMin['memory']=int(line[1])
        except ValueError:
            print('MOLCAS memory does not evaluate to numerical value!')
            sys.exit(54)
    else:
        print('WARNING: Please set memory for MOLCAS in MOLCAS.resources (in MB)! Using 500 MB default value!')
    os.environ['MOLCASMEM']=str(QMin['memory'])
    os.environ['MOLCAS_MEM']=str(QMin['memory'])

    QMin['ncpu']=1
    line=getsh2caskey(sh2cas,'ncpu')
    if line[0]:
        try:
            QMin['ncpu']=int(line[1])
        except ValueError:
            print('Number of CPUs does not evaluate to numerical value!')
            sys.exit(55)

    QMin['mpi_parallel']=False
    line=getsh2caskey(sh2cas,'mpi_parallel')
    if line[0]:
        QMin['mpi_parallel']=True


    QMin['schedule_scaling']=0.6
    line=getsh2caskey(sh2cas,'schedule_scaling')
    if line[0]:
        try:
            x=float(line[1])
            if 0<x<=2.:
                QMin['schedule_scaling']=x
        except ValueError:
            print('Schedule scaling does not evaluate to numerical value!')
            sys.exit(56)

    QMin['Project']='MOLCAS'
    os.environ['Project']=QMin['Project']

    QMin['delay']=0.0
    line=getsh2caskey(sh2cas,'delay')
    if line[0]:
        try:
            QMin['delay']=float(line[1])
        except ValueError:
            print('Submit delay does not evaluate to numerical value!')
            sys.exit(57)

    QMin['Project']='MOLCAS'
    os.environ['Project']=QMin['Project']
    os.environ['MOLCAS_OUTPUT']='PWD'

    line=getsh2caskey(sh2cas,'always_orb_init')
    if line[0]:
        QMin['always_orb_init']=[]
    line=getsh2caskey(sh2cas,'always_guess')
    if line[0]:
        QMin['always_guess']=[]
    if 'always_orb_init' in QMin and 'always_guess' in QMin:
        print('Keywords "always_orb_init" and "always_guess" cannot be used together!')
        sys.exit(58)

    # open template
    template=readfile('MOLCAS.template')

    QMin['template']={}
    integers=['nactel','inactive','ras2','frozen']
    strings =['basis','method','baslib']
    floats=['ipea','imaginary','gradaccumax','gradaccudefault','displ', 'rasscf_thrs_e', 'rasscf_thrs_rot', 'rasscf_thrs_egrd','cholesky_accu']
    booleans=['cholesky','no-douglas-kroll','qmmm','cholesky_analytical']
    for i in booleans:
        QMin['template'][i]=False
    QMin['template']['roots'] = [0 for i in range(8)]
    QMin['template']['rootpad'] = [0 for i in range(8)]
    QMin['template']['method']='casscf'
    QMin['template']['baslib']=''
    QMin['template']['ipea']=0.25
    QMin['template']['imaginary']=0.00
    QMin['template']['frozen']=-1
    QMin['template']['gradaccumax']=1.e-2
    QMin['template']['gradaccudefault']=1.e-4
    QMin['template']['displ']=0.005
    QMin['template']['cholesky_accu']=1e-4
    QMin['template']['rasscf_thrs_e']=1e-8
    QMin['template']['rasscf_thrs_rot']=1e-4            # TODO: apparent default in MOLCAS is 0.1
    QMin['template']['rasscf_thrs_egrd']=1e-4
    QMin['template']['pcmset']={'solvent':'water', 'aare':0.4,'r-min':1.0,'on':False}
    QMin['template']['pcmstate']=(QMin['statemap'][1][0],QMin['statemap'][1][1])


    for line in template:
        orig=re.sub('#.*$','',line).split(None,1)
        line=re.sub('#.*$','',line).lower().split()
        if len(line)==0:
            continue
        if 'spin' in line[0]:
            QMin['template']['roots'][int(line[1])-1]=int(line[3])
        elif 'roots' in line[0]:
            for i,n in enumerate(line[1:]):
                QMin['template']['roots'][i]=int(n)
        elif 'rootpad' in line[0]:
            for i,n in enumerate(line[1:]):
                QMin['template']['rootpad'][i]=int(n)
        elif 'baslib' in line[0]:
            QMin['template']['baslib']=os.path.abspath(orig[1])
        elif line[0] in integers:
            QMin['template'][line[0]]=int(line[1])
        elif line[0] in booleans:
            QMin['template'][line[0]]=True
        elif line[0] in strings:
            QMin['template'][line[0]]=line[1]
        elif line[0] in floats:
            QMin['template'][line[0]]=float(line[1])
        elif 'pcmset' in line[0]:
            # order: solvent, aare, r-min
            QMin['template']['pcmset']['on']=True
            QMin['template']['pcmset']['solvent']=line[1]
            if len(line)>=3:
                QMin['template']['pcmset']['aare']=float(line[2])
            if len(line)>=4:
                QMin['template']['pcmset']['r-min']=float(line[3])
        elif 'pcmstate' in line[0]:
            QMin['template']['pcmstate']=(int(line[1]),int(line[2]))

    # roots must be larger or equal to states
    for i,n in enumerate(QMin['template']['roots']):
        if i==len(QMin['states']):
            break
        if not n>=QMin['states'][i]:
            print('Too few states in state-averaging in multiplicity %i! %i requested, but only %i given' % (i+1,QMin['states'][i],n))
            sys.exit(59)

    # check rootpad
    for i,n in enumerate(QMin['template']['rootpad']):
        if i==len(QMin['states']):
            break
        if not n>=0:
            print('Rootpad must not be negative!')
            sys.exit(60)

    # condense roots list
    for i in range(len(QMin['template']['roots'])-1,0,-1):
        if QMin['template']['roots'][i]==0:
            QMin['template']['roots'].pop(i)
        else:
            break
    QMin['template']['rootpad']=QMin['template']['rootpad'][:len(QMin['template']['roots'])]

    # check roots versus number of electrons
    #nelec=QMin['template']['inactive']*2+QMin['template']['nactel']
    #for i,n in enumerate(QMin['states']):
        #if n>0:
            #if not (QMin['template']['nactel']+i)%2==0:
                #print 'Number of electrons is %i, but states of multiplicity %i are requested.' % (nelec,i+1)
                #sys.exit(61)

    necessary=['basis','nactel','ras2','inactive']
    for i in necessary:
        if not i in QMin['template']:
            print('Key %s missing in template file!' % (i))
            sys.exit(62)


    # logic checks:
    if QMin['template']['pcmset']['on']:
        if QMin['template']['qmmm']:
            print('PCM and QM/MM cannot be used together!')

    # QM/MM mode needs Tinker
    if QMin['template']['qmmm'] and not 'tinker' in QMin:
        print('Please set $TINKER or give path to tinker in MOLCAS.resources!')
        sys.exit(63)

    # find method
    allowed_methods=['casscf','caspt2','ms-caspt2']
    for i,m in enumerate(allowed_methods):
        if QMin['template']['method']==m:
            QMin['method']=i
            break
    else:
        print('Unknown method "%s" given in MOLCAS.template' % (QMin['template']['method']))
        sys.exit(64)

    # decide which type of gradients to do:
    # 0 = analytical CASSCF gradients in one MOLCAS input file (less overhead, but recommended only under certain circumstances)
    # 1 = analytical CASSCF gradients in separate MOLCAS inputs, possibly distributed over several CPUs (DEFAULT)
    # 2 = numerical gradients (CASPT2, MS-CASPT2, Cholesky-CASSCF; or for dmdr and socdr), possibly distributed over several CPUs
    if 'dmdr' in QMin or 'socdr' in QMin or 'grad' in QMin or 'nacdr' in QMin:
        if 'dmdr' in QMin or 'socdr' in QMin:
            QMin['gradmode']=2
        elif QMin['template']['cholesky'] and not QMin['template']['cholesky_analytical']:
            QMin['gradmode']=2
        elif QMin['template']['pcmset']['on']:
            QMin['gradmode']=2
        elif QMin['method']>0:
            QMin['gradmode']=2
        else:
            if QMin['ncpu']>0:
                QMin['gradmode']=1
            else:
                # check if any gradient to be calculated is a SS-CASSCF gradient
                ss_grads=False
                for i in QMin['gradmap']:
                    if QMin['template']['roots'][i[0]-1]==1:
                        ss_grads=True
                if ss_grads:
                    QMin['gradmode']=1
                else:
                    QMin['gradmode']=0
        if QMin['gradmode']==2:
            QMin['displ']=QMin['template']['displ']/au2a
    else:
        QMin['gradmode']=0
    QMin['ncpu']=max(1,QMin['ncpu'])

    # currently, QM/MM is only allowed with CASSCF and analytical gradients on one CPU
    # will be available in the future
    if QMin['template']['qmmm'] and QMin['gradmode']==2:
            print('QM/MM is only possible currently with analytical gradients.')
            sys.exit(65)
    if 'nacdr' in QMin and QMin['gradmode']==2:
            print('Nonadiabatic coupling vectors are only possible currently with analytical gradients.')
            sys.exit(66)
    if QMin['template']['qmmm'] and 'nacdr' in QMin:
        print('Nonadiabatic coupling vectors are currently not available with QM/MM.')
        sys.exit(67)


    # QM/MM setup
    if QMin['template']['qmmm']:
        QMin['active_qmmm_atoms'] = []
        keycontent=readfile('MOLCAS.qmmm.key')
        for line in keycontent:
            if 'qmmm' in line.lower() and not 'electrostatics' in line.lower():
                QMin['total_qmmm_natom'] = int(line.split()[1])
            elif containsstring('^qm ',line.lower()) or containsstring('^mm ',line.lower()):
                line=line.split()
                if len(line)==3 and int(line[1])<0:     # range definition (e.g. QM -1 4)
                    start=-int(line[1])
                    end=int(line[2])
                    QMin['active_qmmm_atoms']+=[i for i in range(start,end+1)]
                else:                                   # list definition (e.g. MM 5 6 7 8)
                    QMin['active_qmmm_atoms']+=[int(i) for i in line[1:]]
        #print 'total amount of qmmm atoms given in key file:', QMin['total_qmmm_natom']
        #print 'number of indices given in key file:', len(QMin['active_qmmm_atoms'])


    # Check the save directory
    try:
        ls=os.listdir(QMin['savedir'])
        err=0
    except OSError:
        print('Problems reading SCRADIR=%s' % (QMin['savedir']))
        sys.exit(68)
    if 'init' in QMin:
        err=0
    elif 'samestep' in QMin:
        for imult,nstates in enumerate(QMin['states']):
            if nstates<1:
                continue
            if not 'MOLCAS.%i.JobIph' % (imult+1) in ls:
                print('File "MOLCAS.%i.JobIph" missing in SAVEDIR!' % (imult+1))
                err+=1
        if 'overlap' in QMin:
            for imult,nstates in enumerate(QMin['states']):
                if nstates<1:
                    continue
                if not 'MOLCAS.%i.JobIph.old' % (imult+1) in ls:
                    print('File "MOLCAS.%i.JobIph.old" missing in SAVEDIR!' % (imult+1))
                    err+=1
    elif 'overlap' in QMin:
        for imult,nstates in enumerate(QMin['states']):
            if nstates<1:
                continue
            if not 'MOLCAS.%i.JobIph' % (imult+1) in ls:
                print('File "MOLCAS.%i.JobIph" missing in SAVEDIR!' % (imult+1))
                err+=1
    if err>0:
        print('%i files missing in SAVEDIR=%s' % (err,QMin['savedir']))
        sys.exit(69)

    QMin['version']=getversion( ['']*50 ,QMin['molcas'] )

    # if PRINT:
    #     printQMin(QMin)

    return QMin



# ========================== Main Code =============================== #
# def main():

# Retrieve PRINT and DEBUG
# try:
#     envPRINT=os.getenv('SH2CAS_PRINT')
#     if envPRINT and envPRINT.lower()=='false':
#         global PRINT
#         PRINT=False
#     envDEBUG=os.getenv('SH2CAS_DEBUG')
#     if envDEBUG and envDEBUG.lower()=='true':
#         global DEBUG
#         DEBUG=True
# except ValueError:
#     print('PRINT or DEBUG environment variables do not evaluate to numerical values!')
#     sys.exit(90)

# # Process Command line arguments
# if len(sys.argv)!=2:
#     print('Usage:\n./SHARC_MOLCAS.py <QMin>\n')
#     print('version:',version)
#     print('date:',versiondate)
#     print('changelog:\n',changelogstring)
#     sys.exit(91)
# QMinfilename=sys.argv[1]
qm_dir = '/home/emielkoridon/PhD_stuff/sharc_saoovqe/formaldimine_test_2/scratch/TRAJ/QM'
os.chdir(qm_dir)
QMinfilename=qm_dir + '/QM.in'

# Print header
printheader()

# Read QMinfile
QMin=readQMin(QMinfilename)

    # # make list of jobs
    # QMin,joblist=generate_joblist(QMin)

    # # run all MOLCAS jobs
    # errorcodes=runjobs(joblist,QMin)

    # # get output
    # QMoutall=collectOutputs(joblist,QMin,errorcodes)

    # # extract data, perform Dyson calculations
    # if 'ion' in QMin:
    #     QMoutDyson=do_Dyson(QMin)
    # else:
    #     QMoutDyson=None

    # # format final output
    # QMout=arrangeQMout(QMin,QMoutall,QMoutDyson)

    # # Measure time
    # runtime=measuretime()
    # QMout['runtime']=runtime

    # # Write QMout
    # writeQMout(QMin,QMout,QMinfilename)

    # # Remove Scratchfiles from SCRATCHDIR
    # if not DEBUG:
    #     cleanupSCRATCH(QMin['scratchdir'])
    #     if 'cleanup' in QMin:
    #         cleanupSCRATCH(QMin['savedir'])
    # if PRINT or DEBUG:
    #     print('#================ END ================#')

# if __name__ == '__main__':
#     main('/home/emielkoridon/PhD_stuff/sharc_saoovqe/QM_test/QM.in')
