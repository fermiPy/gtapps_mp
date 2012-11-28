#!/usr/bin/env python

from multiprocessing import Pool
import numpy as np
import pyfits
import tempfile
import os
import subprocess

from gt_apps import filter,expCube

def ltcube(times):

    '''This is the atomic function that actually runs in the seperate
    threads.  It takes a list as input where the first element is
    tmin, second is tmax, third is spacecraft file, fourth is the
    event file and fifth is the zmax parameter.  It first uses
    gtselect with wide open cuts to divide up the event file then it
    runs gtltcube on that event file.  The temporary event file is
    deleted automatically.  The function returns the name of the
    created ltcube file which can be combined with other files and/or
    deleted later.'''

    print "Starting calculation on interval {} to {}".format(times[0],times[1])
    evfile = tempfile.NamedTemporaryFile(suffix=".fits")
    filter['rad'] = "INDEF"
    filter['evclass'] = 0
    filter['evclsmin'] = 0
    filter['evclsmax'] = 10
    filter['infile'] = times[3]
    filter['outfile'] = evfile.name
    filter['ra'] = "INDEF"
    filter['dec'] = "INDEF"
    filter['tmin'] = times[0]
    filter['tmax'] = times[1]
    filter['emin'] = 0
    filter['emax'] = 400000
    filter['zmax'] = 180
    filter['convtype'] = -1
    filter['chatter'] = 0
    filter.run(print_command=False)

    osfilehandle,outfilename = tempfile.mkstemp(suffix=".fits")
    expCube['evfile'] = evfile.name
    expCube['scfile'] =  times[2]
    expCube['outfile'] = outfilename
    expCube['dcostheta'] = 0.025
    expCube['binsz'] = 1
    expCube['phibins'] = 0
    expCube['zmax'] = times[4]
    expCube['chatter'] = 0
    expCube.run(print_command=False)
    print "Completed calculation on interval {} to {}".format(times[0],times[1])
    return outfilename


def ltsum(filenames, Outfile, SaveTemp):

    '''This function takes a list of livetime cubes and sums them up using
    gtltsum.  It first checks to see if there's only one temporary
    file.  If so, it just copies that to the output file.  If not, it
    creates a temporary file that lists the individual ltcube files
    and operates gtltsum on them.'''

    if len(filenames) <= 1:
        subprocess.call(["cp", filenames[0], Outfile])
    else:
        fileListfile = tempfile.NamedTemporaryFile()
        for filename in filenames:
            fileListfile.file.write(filename + "\n")
        fileListfile.flush()
        subprocess.call(["gtltsum", 
                         "infile1=@"+fileListfile.name, 
                         "outfile="+Outfile])

    if SaveTemp:
        print "Did not delete the following temporary files:"
        print filenames
    else:
        print "Deleting temporary files..."
        for filename in filenames:
            os.remove(filename)


def gtltcube_mp(bins, SCFile, EVFile, OutFile, SaveTemp, zmax):

    '''This functions looks at a spacecraft file and splits the time into
    chunks that match the bin edges in the spacecraft file.  It then
    submits jobs based upon those start and stop times.  This is to
    make the resulting files as close to the original as possible.
    Note that this assumes you are using the full time period in your
    spacecraft file.'''

    print "Opening SC file to determine break points..."
    hdulist = pyfits.open(SCFile)
    scdata = hdulist[1].data
    hdulist.close()
    scstart = scdata.field('START')
    scstop = scdata.field('STOP')
    scstartssplit = np.array_split(scstart,int(bins))
    scstopssplit = np.array_split(scstop,bins) 
    starts = [st[0] for st in scstartssplit]
    stops = [st[-1] for st in scstopssplit]
    scfiles = [SCFile for st in scstartssplit]
    evfiles = [EVFile for st in scstartssplit]
    zmaxes =  [zmax for st in scstartssplit]

    pool = Pool(processes=bins)      
    times = np.array([starts,stops,scfiles,evfiles,zmaxes])
    print "Spawning {} jobs...".format(bins)
    tempfilenames = pool.map(ltcube,times.transpose())
    print "Combining temporary files..."
    ltsum(tempfilenames, OutFile, SaveTemp)

def cli():

    helpString = "Submits the gtltcube program as sperate threads via python and\
                  joins up the resulting temporary exposure cubes at the end\
                  resulting in a single exposure cube for the input event file.\
                  This greatly reduces the running time. For more details on \
                  gtltcube see the gtltcube help file."

    import argparse

    parser = argparse.ArgumentParser(description=helpString)
    parser.add_argument("jobs", type=int, help="The number of jobs you wish to spawn (usually the number of cores on your machine).")
    parser.add_argument("sfile", help="The spacecraft data file. See gtltcube help for more information.")
    parser.add_argument("evfile", help="Input event file.  See gtltcube help for more information.")
    parser.add_argument("outfile", help="Output file name.")

    parser.add_argument("--savetmp", default = False, help="Save the temporary files (default is False).")
    parser.add_argument("--zmax", type=int, default = 180, help="zmax parameter for gtltcube (default is 180)")
    
    args = parser.parse_args()

    gtltcube_mp(args.jobs, args.sfile, args.evfile, args.outfile, args.savetmp, args.zmax)


if __name__ == '__main__': cli()

    
