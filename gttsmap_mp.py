#!/usr/bin/env python

from multiprocessing import Pool
import numpy as np
import tempfile
import os
import subprocess
import pyfits
import pywcs

from gt_apps import TsMap

def tsmap(square):

    '''This is the atomic function that actually runs in the seperate
    threads.  The input is a python list.  The first element is a list
    of two arrays giving the x and y points.  The second through last
    elements are the scfile, evfile, expcube, nlong, nlat, irfs,
    srcrad and nenergies.  This function creates a temporary file for
    the output and returns that file's name and the x and y pixel coords.'''

    print "Starting analysis of pixel {}...".format(square[0][2:])
    osfilehandle,outfilename = tempfile.mkstemp(suffix=".fits")
    TsMap['xref'] = square[0][0]
    TsMap['yref'] = square[0][1]
    TsMap['nxpix'] = 1 
    TsMap['nypix'] = 1 
    TsMap['evfile'] = square[1]
    TsMap['scfile'] = square[2]
    TsMap['statistic']  = square[3]
    TsMap['expmap'] = square[4]
    TsMap['expcube'] = square[5]
    TsMap['srcmdl'] = square[6]
    TsMap['irfs'] = square[7]
    TsMap['optimizer'] = square[8]
    TsMap['ftol'] = square[9]
    TsMap['binsz'] = square[10]
    TsMap['coordsys'] = square[11]
    TsMap['proj'] = square[12]
    TsMap['outfile'] = outfilename
    TsMap['chatter'] = 0
    TsMap.run(print_command=False)
    print "Completed analysis of pixel {}".format(square[0][2:])
    return outfilename,square[0][2:]

def tsmapsum(files_info, out_hdu, Outfile, SaveTemp):

    '''This function takes a list of tsmaps maps and mosaics them
    together.  If there is only one file to be added, it just copies
    it to the outfile.  Each input file contains a single pixel and
    the second element of the files_info array is the location of that
    pixel in the larger array.'''
    
    if len(files_info) <= 1:
        subprocess.call(["cp", files_info[0][0], Outfile])
    else:
        for file_info in files_info:
            hdu = pyfits.open(file_info[0])
            out_hdu.data[(int(file_info[1][1] - 1),int(file_info[1][0] - 1))] = hdu[0].data
            hdu.close()
        out_hdu.writeto(Outfile, clobber="True")

    if SaveTemp:
        print "Did not delete the following temporary files:"
        print [file_info[0] for file_info in files_info]
    else:
        print "Deleting temporary files..."
        for file_info in files_info:
            os.remove(file_info[0])

def gttsmap_mp(nxpix, nypix, jobs, evfile, sfile, statistic,
               expmap, expcube, srcmdl, IRFS, optimizer, ftol, binsz, 
               coordsys, xref, yref, proj, outfile, savetmp):

    '''This function figures out how to divide up the image into pixels
    and submits each of the pixels as a job to a pool of gttsmap
    functions.  It creates a dummy WCS to find the sky coordinates of each
    pixel and also creats a new FITS HDU that can store the final
    values.'''

    print "THIS PROGRAM DOES NOT PRODUCE VALID RESULTS"

    wcs_temp = pywcs.WCS(naxis=2)
    wcs_temp.wcs.crpix = [(nxpix+1)/2.,(nypix+1)/2.]
    wcs_temp.wcs.cdelt = np.array([-1*binsz,binsz])
    wcs_temp.wcs.crval = [xref, yref]
    if coordsys == 'CEL':
        wcs_temp.wcs.ctype = ["RA---{}".format(proj), "DEC--{}".format(proj)]
    else:
        wcs_temp.wcs.ctype = ["GLAT-{}".format(proj), "GLONG-{}".format(proj)]
    wcs_temp_header = wcs_temp.to_header()
    out_hdu = pyfits.PrimaryHDU(data=np.zeros((nxpix,nypix)),header=wcs_temp_header)

    submaps_crpix = [(x,y) for x in np.arange(1,nxpix+1) for y in np.arange(1,nypix+1)]
    submaps_crval = wcs_temp.wcs_pix2sky(submaps_crpix,1)
    submaps_geometry = np.column_stack((submaps_crval,submaps_crpix))

    SQ = [(row, evfile, sfile, statistic,
           expmap, expcube, srcmdl, IRFS, optimizer, ftol, binsz, 
           coordsys, proj) for row in submaps_geometry]

    pool = Pool(processes=jobs)      
    print "Spawning {} pixels on {} nodes...".format(nxpix*nypix,jobs)
    files_info = pool.map(tsmap,SQ)
    print "Combining temporary files..."
    tsmapsum(files_info, out_hdu, outfile, savetmp)

def cli():

    helpString = "  -> DO NOT USE THIS PROGRAM. RESULTS ARE NOT VALID. <- \
                  Submits the gtexpmap program as sperate threads via python and\
                  joins up the resulting temporary exposure maps at the end\
                  resulting in a single exposure map for the input event file.\
                  This greatly reduces the running time. For more details on \
                  gtexpmap see the gtexpmap help file.  Note that the checksum \
                  and datasum are incorrect for the final file."

    import argparse

    parser = argparse.ArgumentParser(description=helpString)
    parser.add_argument("nxpix", type=int, help="Number of pixels along x-axis.  See gttsmap help for more information.")
    parser.add_argument("nypix", type=int, help="Number of pixels along y-axis.  See gttsmap help for more information.")
    parser.add_argument("jobs", type=int, help="The number of concurrent jobs.")
    parser.add_argument("evfile", help="Input event file.  See gttsmap help for more information.")
    parser.add_argument("sfile", help="The spacecraft data file. See gttsmap help for more information.")
    parser.add_argument("statistic", help="UNBINNED or BINNED. See gttsmap help for more information.")
    parser.add_argument("expmap", help="Input exposure map.  See gttsmap help for more information.")    
    parser.add_argument("expcube", help="Input livetime cube.  See gttsmap help for more information.")
    parser.add_argument("srcmdl", help="XML source model definition.  See gttsmap help for more information.")
    parser.add_argument("IRFS", help="IRFs to use.  See gttsmap help for more information.")
    parser.add_argument("optimizer", help="The optimizer (e.g. NEWMINUIT). See gttsmap help for more information.")
    parser.add_argument("ftol", type=float, help="Relative Fit tolerance. See gttsmap help for more information.")
    parser.add_argument("binsz", type=float, help="Image scale (deg/pix).  See gttsmap help for more information.")
    parser.add_argument("coordsys", help="CEL or GAL.  See gttsmap help for more information.")
    parser.add_argument("xref", type=float, help="x-coord of center (RA or l).  See gttsmap help for more information.")
    parser.add_argument("yref", type=float, help="y-coord of center (DEC or b).  See gttsmap help for more information.")
    parser.add_argument("proj", help="Coordinate projection. See gttsmap help for more information.")
    parser.add_argument("outfile", help="Output file name.")

    parser.add_argument("--savetmp", default = False, help="Save the temporary files (default is False).")
    
    args = parser.parse_args()

    gttsmap_mp(args.nxpix, args.nypix, args.jobs, args.evfile, args.sfile, args.statistic,
               args.expmap, args.expcube, args.srcmdl, args.IRFS, args.optimizer, args.ftol, args.binsz, 
               args.coordsys, args.xref, args.yref, args.proj, args.outfile, args.savetmp)

if __name__ == '__main__': cli()

    
