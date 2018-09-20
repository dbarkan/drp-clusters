from modeller import *
import modeller.salign
import sys
import tempfile
import os
import subprocess
import argparse
import align_native_overlap
import align_disulfides
import cluster_lib

parser = argparse.ArgumentParser()
parser.add_argument("-q", "--drp_query_file", help="Text file with input set of DRPs. One DRP per line, specified as a DRP code (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)", required=True)
parser.add_argument("-p", "--pdb_directory", required=True, help="Location of PDB files. Expected format is identical to the 'divided' PDB FTP site at ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/")
parser.add_argument("-o", "--output_file", help="Full path of output file. If not set, output will be sent to STDOUT")
parser.add_argument("-m", "--align_mode", choices=['full_drp', 'disulfides'], required=True, help="which type of alignment to make (see documentation for details)")
parser.add_argument("-x", "--max_drp_count", help="For testing purposes; limit the number of DRPs to align to this number", type=int)

if (len(sys.argv) < 2):
    print "Please run with '-h' for full usage"
    sys.exit()
config = parser.parse_args(sys.argv[1:])

    
fh = open(config.drp_query_file, 'r')
drpList = []
for line in fh:
    line = line.rstrip('\n\r')
    if (line.strip() == ''):
        continue
    drpList.append(line)
    
outputFh = open(config.output_file, 'w') #initialize output file; otherwise anything existing will be appended to
outputFh.close()

cutoff = len(drpList)
if (config.max_drp_count):
    cutoff = config.max_drp_count
    
for i in range(0, cutoff):
    print 'aligning DRP %s (%s of %s)' % (drpList[i], i + 1, cutoff)
    for j in range(i+1, cutoff):
        pnoa = None
        if (config.align_mode == 'full_drp'):
            pnoa = align_native_overlap.PairwiseNativeOverlapAligner()
        else:
            pnoa = align_disulfides.PairwiseDisulfideAligner()

        iDrp = drpList[i]
        jDrp = drpList[j]
        
        pnoa.setParams(iDrp, jDrp, config.pdb_directory, config.output_file, True)
        try:
            pnoa.execute()
        except cluster_lib.ModellerException, e: #these are occasional errors inherent in modeller, safe to ignore. Other Exceptions will not be caught
            pass



