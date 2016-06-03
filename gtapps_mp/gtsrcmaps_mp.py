#!/usr/bin/env python

from multiprocessing import Pool
import numpy as np
import tempfile
import os
import subprocess
from gtapps_mp.utils import pyfits

from gt_apps import evtbin,srcMaps

def run_gtsrcmaps(square):

    '''This is the atomic function that actually runs in the seperate
    threads.  The input is an array with two elements, the first describes
    the energy binning and the second lists all of the progromatic
    options.  This function runs evtbin to generate a counts cube for the
    sub-energy range and then gtsrcmaps.  It returns the bin number and
    the output filename of the srcmaps file for the sub-energy range.'''


    idx = square[0][0][0]
    emin = square[0][0][1]*0.001
    emax = square[0][-1][2]*0.001
    enumbins = len(square[0])
    options = square[1]

    ccubehandle,ccubefilename = tempfile.mkstemp(suffix=".fits")
    evtbin['evfile'] = options['evfile']
    evtbin['outfile'] = ccubefilename
    evtbin['algorithm'] = 'CCUBE'
    evtbin['nxpix'] = options['nxpix']
    evtbin['nypix'] = options['nypix']
    evtbin['binsz'] = options['binsz']
    evtbin['coordsys'] = options['coordsys']
    evtbin['xref'] = options['xref']
    evtbin['yref'] = options['yref']
    evtbin['axisrot'] = options['axisrot']
    evtbin['proj'] = options['proj']
    evtbin['ebinalg'] = options['ebinalg']
    evtbin['emin'] = emin
    evtbin['emax'] = emax
    evtbin['enumbins'] = enumbins
    evtbin['scfile'] = options['scfile']
    evtbin.run(print_command=False)

    srcmaphandle,srcmapfilename = tempfile.mkstemp(suffix=".fits")
    srcMaps['scfile'] = options['scfile']
    srcMaps['expcube'] = options['expcube']
    srcMaps['cmap'] = ccubefilename
    srcMaps['srcmdl'] = options['srcmdl']
    srcMaps['bexpmap'] = options['bexpmap']
    srcMaps['outfile'] = srcmapfilename
    srcMaps['irfs'] = options['irfs']
    srcMaps['chatter'] = 2
    srcMaps['emapbnds'] = "no"
    srcMaps.run(print_command=False)

    return idx,srcmapfilename

def srcmapssum(results, ref_hdu, outfile, savetemp):

    '''This function takes a list of source maps and joins them together.
    If there is only one file to be joined, it just copies it to the
    outfile.  If there is more than one, it uses pyfits to open them
    all up, add the energy slices together and then outputs this.  The
    first element in the list is the index.'''

    results.sort()

    if len(results) <= 1:
        subprocess.call(["cp", results[0][1], outfile])
    else:
        SrcMaps = [pyfits.open(result[1]) for result in results]
        for idx,sourceHDU in enumerate(SrcMaps[0][3:]):
            newHDU = SrcMaps[0][idx+3].copy()
            newHDU.data = np.delete(SrcMaps[0][idx+3].data,-1,0)
            for SrcMap in SrcMaps[1:]:
                newHDU.data = np.append(newHDU.data,SrcMap[idx+3].data[:-1],axis=0)
            newHDU.data = np.append(newHDU.data,[SrcMaps[-1][idx+3].data[-1]],axis=0)
            #newHDU.add_checksum()
            ref_hdu.append(newHDU)
        ref_hdu.writeto(outfile,clobber='yes')
        
    if savetemp:
        print "Did not delete the following temporary files:"
        print results
    else:
        print "Deleting temporary files..."
        for result in results:
            os.remove(result[1])

def gtsrcmaps_mp(nxpix, nypix, binsz, scfile, evfile, expcube, 
                 ccube, coordsys, proj, xref, yref, axisrot, ebinalg,
                 emin, emax, enumbins,jobs, srcmdl, bexpmap, irfs, 
                 outfile, savetmp):

    '''This function opens up the input counts map and divides up the
    energy bins into seperate jobs which are sent to the pool to be
    run.'''

    ref_hdu = pyfits.open(ccube)
    energies = np.array(ref_hdu[1].data)
    
    #There's an error when you only use one bin so you need to
    #make sure and compute more than one bin per job.
    while len(energies)/float(jobs) < 3:
        print "Too many jobs ({}), reducing by 1".format(jobs)
        jobs -= 1
        print "Jobs is now {}".format(jobs)
    energy_arrays = np.array_split(energies,jobs)

    options = {'nxpix': nxpix,
               'nypix': nypix,
               'binsz': binsz,
               'scfile': scfile,
               'evfile': evfile,
               'coordsys': coordsys,
               'xref': xref,
               'yref': yref,
               'axisrot': axisrot,
               'proj': proj,
               'ebinalg': ebinalg,
               'emin': emin,
               'emax': emax,
               'enumbins': enumbins,
               'expcube': expcube,
               'srcmdl': srcmdl,
               'bexpmap': bexpmap,
               'irfs': irfs}

    SQ = [(array,options) for array in energy_arrays]

    pool = Pool(processes=jobs)
    print "Spawning {} jobs...".format(jobs)
    results = pool.map(run_gtsrcmaps,SQ)
    print "Combining temporary files..."
    srcmapssum(results, ref_hdu, outfile, savetmp)

def cli():

    helpString = "Submits the gtsrcmaps program as sperate threads via python and\
                  joins up the resulting temporary source maps at the end\
                  resulting in a single source map file.\
                  This greatly reduces the running time. For more details on \
                  gtsrcmaps see the gtsrcmaps help file.  Note that the checksum \
                  and datasum are incorrect for the final file."

    import argparse

    parser = argparse.ArgumentParser(description=helpString)
    
    parser.add_argument('jobs', type=int, help="Number of concurrent jobs (threads) to run.")
    parser.add_argument('nxpix', type=int, help="Number of horizontal pixels.  See gtbin for more information.")
    parser.add_argument('nypix', type=int, help="Number of vertical pixels.  See gtbin for more information.")
    parser.add_argument('binsz', type=float, help="Degrees per pixel.  See gtbin for more information.")
    parser.add_argument('scfile', help="Spacecraft file.  See gtsrcmaps for more information.")
    parser.add_argument('evfile', help="Event file.  See gtbin for more information.")
    parser.add_argument('coordsys', help="Coordinate system (CEL or GAL).  See gtbin for more information.")
    parser.add_argument('xref', type=float, help="Horizontal position. See gtbin for more information.")
    parser.add_argument('yref', type=float, help="Vertical position. See gtbin for more information.")
    parser.add_argument('axisrot', type=float, help="Rotation angle.  See gtbin for more information.")
    parser.add_argument('proj', help="Coordintate projection. See gtbin for more information.")
    parser.add_argument('ebinalg', help="Energy bin specification.  See gtbin for more information.")
    parser.add_argument('emin', help="Minimum energy (in MeV). See gtbin for more information.")
    parser.add_argument('emax', help="Maximum energy (in MeV). See gtbin for more information.")
    parser.add_argument('enumbins', help="Number of energy bins. See gtbin for more information.")
    parser.add_argument('expcube', help="Output of gtltcube. See gtsrcmaps for more information.")
    parser.add_argument('srcmdl', help="XML model file.   See gtsrcmaps for more information.")
    parser.add_argument('bexpmap', help="Binned exposure map from gtexpcube2.   See gtsrcmaps for more information.")
    parser.add_argument('ccube', help="Counts cube produced by gtbin.    See gtsrcmaps for more information.")
    parser.add_argument('irfs', help="Name of the IRFs.   See gtsrcmaps for more information.")
    parser.add_argument('outfile', help="Name of the output FITS file")
    
    parser.add_argument("--savetmp", default = False, help="Save the temporary files (default is False).")
    
    args = parser.parse_args()

    gtsrcmaps_mp(args.nxpix, args.nypix, args.binsz, args.scfile, args.evfile, args.expcube, 
                 args.ccube, args.coordsys, args.proj, args.xref, args.yref, args.axisrot, args.ebinalg,
                 args.emin, args.emax, args.enumbins, args.jobs, args.srcmdl, args.bexpmap, args.irfs, 
                 args.outfile, args.savetmp)

if __name__ == '__main__': cli()

    
