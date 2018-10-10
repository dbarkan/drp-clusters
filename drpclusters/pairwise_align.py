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
import logging
logger = logging.getLogger('pairwise_align')
logger.setLevel(logging.DEBUG)


class PairwiseAlignPipeline(cluster_lib.Runner):
    def execute(self):
        drpList = cluster_lib.readDrpCodeFile(self.config.drp_query_file)    

        outputFh = open(self.config.output_file, 'w') #initialize output file; otherwise anything existing will be appended to
        outputFh.close()

        cutoff = len(drpList)
        if (self.config.max_drp_count):
            cutoff = self.config.max_drp_count

        for i in range(0, cutoff):
            logger.info('Aligning DRP %s (%s of %s)' % (drpList[i], i + 1, cutoff))
            for j in range(i+1, cutoff):
                pnoa = None
                iDrp = drpList[i]
                jDrp = drpList[j]
                
                if (self.config.align_mode == 'full_drp'):
                    parser = align_native_overlap.getParser()
                    config = parser.parse_args(['-f', iDrp, '-s', jDrp, '-p', self.config.pdb_directory, '-o', self.config.output_file, '-a'])
                    pnoa = align_native_overlap.PairwiseNativeOverlapAligner(config)
                else:
                    parser = align_disulfides.getParser()
                    config = parser.parse_args(['-f', iDrp, '-s', jDrp, '-p', self.config.pdb_directory, '-o', self.config.output_file, '-a'])
                    pnoa = align_disulfides.PairwiseDisulfideAligner(config)


                try:
                    pnoa.execute()
                except cluster_lib.ModellerException, e: #these are occasional errors inherent in modeller, safe to ignore. Other Exceptions will not be caught
                    pass

    def makeLogPrefix(self):
        return 'pairwiseAlign'

    
def getParser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Iterate over pairs of DRP PDB files and align them with one of two methods")
    
    parser.add_argument("-q", dest="drp_query_file", metavar="<file>", required=True,
                        help="Text file with input set of DRPs. One DRP per line, specified as a DRP code\n"
                        "(5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)\n\n")

    parser.add_argument("-p", dest="pdb_directory", metavar="<directory>", required=True,
                        help="Location of PDB files. Created in setup_pdb.py; each file should be named after its DRP code\n\n")

    parser.add_argument("-r", dest="run_directory", metavar='<dir>', required=True,
                        help="Directory to which log and other output is written. Will be created if does not exist\n\n")

    parser.add_argument("-o", dest="output_file", metavar="<file>",
                        help="Full path of output file.\n\n")
    
    parser.add_argument("-m", dest="align_mode", metavar="<string>", choices=['full_drp', 'disulfides'], required=True,
                        help="which type of alignment to make (see documentation for details)\n\n")
    
    parser.add_argument("-x", dest="max_drp_count", metavar="<int>", type=int,
                        help="For testing purposes; limit the number of DRPs to align to this number\n\n")
    return parser

if __name__ == '__main__':

    parser = getParser()
    config = parser.parse_args(sys.argv[1:])

    runner = PairwiseAlignPipeline(config)
    runner.execute()

