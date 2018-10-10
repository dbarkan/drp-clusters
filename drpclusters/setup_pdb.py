from modeller import *
import sys
import tempfile
import os
import subprocess
import argparse
import gzip
import cluster_lib
import logging

logger = logging.getLogger('setup_pdb')
logger.setLevel(logging.DEBUG)

class MissingPdbFileException(cluster_lib.DrpClusterException):
    pass

class PdbCopy(cluster_lib.Runner):
    def makeMirrorPdb(self, pdbId):
        return os.path.join(self.config.pdb_directory, pdbId[1:3], "pdb%s.ent.gz" % pdbId)

    def makeLogPrefix(self):
        return 'setupPdb'
    
    def execute(self):
        if (not os.path.exists(self.config.pdb_output_directory)):
            os.mkdir(self.config.pdb_output_directory)
        fh = open(self.config.drp_query_file, 'r')
        counter = 0
        lengthOutputFh = open(self.getFullOutputFile(os.path.basename(config.drp_length_output_file)), "w")
        for drpCode in cluster_lib.readDrpCodeFile(self.config.drp_query_file):
            counter += 1
            if (counter % 10 == 0):
                logger.info('Copying DRP %s' % counter)

            [pdbId, chainId] = cluster_lib.readDrpCode(drpCode)
            
            fullPdb = self.makeMirrorPdb(pdbId)
            if (not os.path.exists(fullPdb)):
                fullPdb = os.path.join(self.config.pdb_directory, "%s.pdb" % pdbId)
                if (not os.path.exists(fullPdb)):
                    raise MissingPdbFileException("Did not find expected PDB file for DRP code %s\n"
                                                  "Searched for the following:\n%s\n%s\n"
                                                  "Please ensure your local PDB mirror is set up according to specifications in the documentation\n"
                                                  "Please also make sure the path you specified to the PDB mirror root directory is correct" % (pdbId, fullPdb, self.makeMirrorPdb(pdbId)))

            #use MODELLER to read coordinate file from pdbDir
            log.none()
            env = environ()
            env.io.atom_files_directory = ['.', self.config.pdb_directory]
            firstModel = model(env, file=pdbId, model_segment=('FIRST:'+chainId, 'LAST:'+chainId))

            #write to length output file for cluster pipeline downstream
            lengthOutputFh.write("%s\t%s\n" % (drpCode, len(firstModel.residues)))
            
            #in local dir, write out temp file with *only* the chain representing the DRP
            #(much easier to use MODELLER to do this rather than parse the PDB file manually)
            firstModel.write(file="%s.temp.pdb" % drpCode)

            #However MODELLER does not retain SSBOND information which we need. So read the temp file back in
            drpPdbFh = open("%s.temp.pdb" % drpCode)
            pdbLines = []
            for drpLine in drpPdbFh:
                pdbLines.append(drpLine)
            drpPdbFh.close()

            #Read the original PDB file to get SSBOND
            fullFh = gzip.open(fullPdb, 'r')
            ssBondLines = []
            for fullLine in fullFh:
                if(fullLine.startswith('SSBOND')):
                    ssBondLines.append(fullLine)

            #concatenate temp file and SS BOND info into final PDB file
            finalDrpFh = open(os.path.join(self.config.pdb_output_directory, "%s.pdb" % drpCode), 'w')
            outputLines = [pdbLines[0]] + ssBondLines + pdbLines[1:]
            finalDrpFh.write("".join(outputLines))
            finalDrpFh.close()

            #remove temp file
            os.remove("%s.temp.pdb" % drpCode)
        lengthOutputFh.close()

def getParser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Copy DRP PDB files to local directory (only save coordinates for chain representing DRP, along with SSBOND info)")

    parser.add_argument("-r", dest="run_directory", metavar='<dir>', required=True,
                        help="Directory to which log and other output is written. Will be created if does not exist\n\n")
    
    parser.add_argument("-q", dest="drp_query_file", metavar='<file>', required=True,
                        help="Text file with input set of DRPs. One DRP per line, specified as a DRP code\n"
                        "(5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)\n\n")
    
    parser.add_argument("-p", dest="pdb_directory", metavar='<dir>', required=True,
                        help="Location of PDB files. Expected format is identical to the 'divided' PDB FTP site at\n"
                        "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/\n\n")
    
    parser.add_argument("-o", dest="pdb_output_directory", metavar='<dir>', required=True,
                        help="Output directory to which single chain PDB files will be written\n\n")
    
    parser.add_argument("-l", dest="drp_length_output_file", metavar='<file>', default="drp_lengths.txt",
                        help="Output file to write DRP sequence length annotation. Do not include full path (default drp_lengths.txt)\n\n")
    return parser
        
if __name__ == '__main__':
    parser = getParser()
    config = parser.parse_args(sys.argv[1:])

    pc = PdbCopy(config)

    pc.execute()

