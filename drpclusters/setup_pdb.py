from modeller import *
import sys
import tempfile
import os
import subprocess
import argparse
import gzip
import cluster_lib

class PdbCopy:
    def copyPdbs(self, inputDrpFile, pdbDir):
        fh = open(inputDrpFile, 'r')
        counter = 0
        lengthOutputFh = open("drp_lengths.txt", "w")
        for line in fh:
            counter += 1
            if (counter % 10 == 0):
                print 'copy drp %s' % counter
            line = line.rstrip('\n\r')
            if (line.strip() == ''):
                continue

            [pdbId, chainId] = cluster_lib.readDrpCode(line)
            fullPdb = os.path.join(pdbDir, pdbId[1:3], "pdb%s.ent.gz" % pdbId)
            if (not os.path.exists(fullPdb)):
                msg =  "ERROR: did not find expected PDB file %s for DRP code %s\n" % (fullPdb, line)
                msg += "Please ensure your local PDB mirror is set up according to specifications in the documentation\n"
                msg += "Please also make sure the path you specified to the PDB mirror root directory is correct"
                print msg
                sys.exit()
            
            #use MODELLER to read coordinate file from pdbDir
            log.none()
            env = environ()
            env.io.atom_files_directory = ['.', pdbDir]
            firstModel = model(env, file=pdbId, model_segment=('FIRST:'+chainId, 'LAST:'+chainId))

            #write to length output file for cluster pipeline downstream
            lengthOutputFh.write("%s\t%s\n" % (line, len(firstModel.residues)))
            
            #in local dir, write out temp file with *only* the chain representing the DRP
            #(much easier to use MODELLER to do this rather than parse the PDB file manually)
            firstModel.write(file="%s.temp.pdb" % line)

            #However MODELLER does not retain SSBOND information which we need. So read the temp file back in
            drpPdbFh = open("%s.temp.pdb" % line)
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
            finalDrpFh = open("%s.pdb" % line, 'w')
            outputLines = [pdbLines[0]] + ssBondLines + pdbLines[1:]
            finalDrpFh.write("".join(outputLines))
            finalDrpFh.close()

            #remove temp file
            os.remove("%s.temp.pdb" % line)
        lengthOutputFh.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Copy DRP PDB files to cwd (only save coordinates for chain representing DRP, along with SSBOND info)")
    parser.add_argument("-q", "--drp_query_file", help="Text file with input set of DRPs. One DRP per line, specified as a DRP code (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)", required=True)
    parser.add_argument("-p", "--pdb_directory", required=True, help="Location of PDB files. Expected format is identical to the 'divided' PDB FTP site at ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/")
    
    if (len(sys.argv) < 2):
        print "Please run with '-h' for full usage"
        sys.exit()
    config = parser.parse_args(sys.argv[1:])
    pc = PdbCopy()
    pc.copyPdbs(config.drp_query_file, config.pdb_directory)

