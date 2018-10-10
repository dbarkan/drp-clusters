from modeller import *
import modeller.salign
import sys
import tempfile
import os
import subprocess
import argparse
import cluster_lib

class PairwiseNativeOverlapAligner:
    def __init__(self, config):
        self.config = config
        
    def writeOutput(self, label, value):
        outputList = [self.config.first_pdb_code, self.config.second_pdb_code, label, value]
        self.resultFh.write("%s\n" % '\t'.join(str(x) for x in outputList))

    def writeErrorAndRaise(self, errorCode, modellerError=None):
        outputList = [self.config.first_pdb_code, self.config.second_pdb_code, 'native_overlap', "Error: %s" % errorCode]
        self.resultFh.write("%s\n" % '\t'.join(str(x) for x in outputList))
        self.resultFh.close()
        
        if (modellerError):
            print "Skipping alignment between %s and %s due to Modeller error: %s" % (self.config.first_pdb_code, self.config.second_pdb_code, str(modellerError))
        raise cluster_lib.ModellerException(errorCode)
        
    def execute(self):
        appendString = 'w'
        if (self.config.append_to_output):
            appendString = 'a'

        self.resultFh = open(self.config.output_file, appendString)
        log.none()
        env = environ()
        env.io.atom_files_directory = ['.', '../atom_files', self.config.pdb_directory]

        [firstPdb, firstChain] = cluster_lib.readDrpCode(self.config.first_pdb_code)
        [secondPdb, secondChain] = cluster_lib.readDrpCode(self.config.second_pdb_code)

        aln = alignment(env)
        firstModel = None
        secondModel = None
        pdbFile = os.path.join(self.config.pdb_directory, "%s.pdb" % self.config.first_pdb_code)
        if (not os.path.exists(pdbFile)):
            print "Warning: did not find expected PDB file %s" % pdbFile
            self.writeErrorAndRaise("missing_pdb")
        try:
            firstModel = model(env, file=self.config.first_pdb_code, model_segment=('FIRST:'+firstChain, 'LAST:'+firstChain))        
        except Exception, e:
            self.writeErrorAndRaise("first_model", e)
            raise e
        aln.append_model(firstModel, atom_files=self.config.first_pdb_code, align_codes=firstPdb+firstChain)

        try:
            secondModel = model(env, file=self.config.second_pdb_code, model_segment=('FIRST:'+secondChain, 'LAST:'+secondChain))
        except Exception, e:
            self.writeErrorAndRaise("second_model", e)
            raise e

        aln.append_model(secondModel, atom_files=self.config.second_pdb_code, align_codes=secondPdb+secondChain)
        #Run SALIGN
        try:
            saveStdout = sys.stdout
            sys.stdout = open(os.devnull, "w")
            modeller.salign.iterative_structural_align(aln)
            sys.stdout = saveStdout
        except Exception, e:
            sys.stdout = saveStdout
            self.writeErrorAndRaise("salign", e)
        
        #Superpose second onto first to get similarity metrics
        atmsel = selection(firstModel).only_atom_types('CA')
        r = atmsel.superpose(secondModel, aln)

        #Prepare fraction calculations
        firstModelLength = len(firstModel.residues)
        secondModelLength = len(secondModel.residues)

        #prepare sequence identity calculations
        firstSequence = aln[0]
        secondSequence = aln[1]
        sequenceIdentity = firstSequence.get_sequence_identity(secondSequence)

        firstFraction = (r.num_equiv_pos * 1.0) / (firstModelLength * 1.0)
        secondFraction = (r.num_equiv_pos * 1.0) / (secondModelLength * 1.0)
        firstCutoffFraction = (r.num_equiv_cutoff_pos * 1.0) / (firstModelLength * 1.0)
        secondCutoffFraction = (r.num_equiv_cutoff_pos * 1.0) / (secondModelLength * 1.0)

        shorterFraction = min(firstFraction, secondFraction)
        shorterCutoffFraction = min(firstCutoffFraction, secondCutoffFraction)
        longerFraction = max(secondFraction, firstFraction)
        longerCutoffFraction = max(secondCutoffFraction, firstCutoffFraction)

        shorterSequenceProduct = sequenceIdentity * shorterCutoffFraction
        longerSequenceProduct = sequenceIdentity * longerCutoffFraction

        self.writeOutput("native_overlap", r.num_equiv_pos)
        self.writeOutput("native_overlap_rms_cutoff", r.num_equiv_cutoff_pos)
        
        self.writeOutput("second_fraction", secondFraction)
        self.writeOutput("first_fraction", firstFraction)
        self.writeOutput("shorter_fraction", shorterFraction)
        self.writeOutput("longer_fraction", longerFraction)

        self.writeOutput("second_cutoff_fraction", secondCutoffFraction)
        self.writeOutput("first_cutoff_fraction", firstCutoffFraction)
        self.writeOutput("shorter_cutoff_fraction", shorterCutoffFraction)
        self.writeOutput("longer_cutoff_fraction", longerCutoffFraction)

        self.writeOutput("sequence_identity", sequenceIdentity)

        self.writeOutput("shorter_sequence_product", shorterSequenceProduct)
        self.writeOutput("longer_sequence_product", longerSequenceProduct)

        self.resultFh.close()

def getParser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Structurally Align two DRPs over their full length")

    parser.add_argument("-f", dest="first_pdb_code", metavar="<string>", required=True,
                        help="DRP code for first coordinate set to align (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)\n\n")
    
    parser.add_argument("-s", dest="second_pdb_code", metavar="<string>", required=True,
                        help="DRP code for second coordinate set to align (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)\n\n")
    
    parser.add_argument("-p", dest="pdb_directory", metavar="<directory>", required=True,
                        help="Location of PDB files. Created in setup_pdb.py; each file should be named after its DRP code\n\n")
    
    parser.add_argument("-o", dest="output_file", metavar="<directory>", required=True, help="Full path of output file\n\n")
    
    parser.add_argument("-a", dest="append_to_output", action="store_true",
                        help="If set, output will be appended to <output_file>; if not then existing <output_file> will be overwritten\n\n")
    
    return parser
        
if __name__ == '__main__':
    
    parser = getParser()

    config = parser.parse_args(sys.argv[1:])

    pnoa = PairwiseNativeOverlapAligner(config)
    pnoa.execute()
