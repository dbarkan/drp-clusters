from modeller import *
import sys
import tempfile
import os
import subprocess
import argparse
import gzip
import cluster_lib

class MissingPdbFileException(cluster_lib.DrpClusterException):
    pass

class PdbCopy:
    def makeMirrorPdb(self, config, pdbId):
        return os.path.join(config.pdb_directory, pdbId[1:3], "pdb%s.ent.gz" % pdbId)

    def copyPdbs(self, config):
        fh = open(config.drp_query_file, 'r')
        counter = 0
        lengthOutputFh = open(config.drp_length_output_file, "w")
        for line in fh:
            counter += 1
            if (counter % 10 == 0):
                print 'copy drp %s' % counter
            line = line.rstrip('\n\r')
            if (line.strip() == ''):
                continue

            [pdbId, chainId] = cluster_lib.readDrpCode(line)
            
            fullPdb = self.makeMirrorPdb(config, pdbId)
            if (not os.path.exists(fullPdb)):
                fullPdb = os.path.join(config.pdb_directory, "%s.pdb" % pdbId)
                if (not os.path.exists(fullPdb)):
                    raise MissingPdbFileException("Did not find expected PDB file for DRP code %s\n"
                                                  "Searched for the following:\n%s\n%s\n"
                                                  "Please ensure your local PDB mirror is set up according to specifications in the documentation\n"
                                                  "Please also make sure the path you specified to the PDB mirror root directory is correct" % (pdbId, fullPdb, self.makeMirrorPdb(config, pdbId)))

            #use MODELLER to read coordinate file from pdbDir
            log.none()
            env = environ()
            env.io.atom_files_directory = ['.', config.pdb_directory]
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
            finalDrpFh = open(os.path.join(config.output_directory, "%s.pdb" % line), 'w')
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
    parser.add_argument("-o", "--output_directory", required=True, help="Output directory to which single chain PDB files will be written")
    parser.add_argument("-l", "--drp_length_output_file", help="Output file to write DRP sequence length annotation", default="drp_lengths.txt")
    if (len(sys.argv) < 2):
        print "Please run with '-h' for full usage"
        sys.exit()
    config = parser.parse_args(sys.argv[1:])
    pc = PdbCopy()
    pc.copyPdbs(config)

