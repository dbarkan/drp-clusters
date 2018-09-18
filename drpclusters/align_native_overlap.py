from modeller import *
import modeller.salign
import sys
import tempfile
import os
import subprocess
import argparse
import cluster_lib

class PairwiseNativeOverlapAligner:

    def writeOutput(self, label, value):
        outputList = [self.first_pdb_code, self.second_pdb_code, label, value]
        self.resultFh.write("%s\n" % '\t'.join(str(x) for x in outputList))

    def setParams(self, first_pdb_code, second_pdb_code, pdb_directory, output_file, append_to_output):
        self.first_pdb_code = first_pdb_code
        self.second_pdb_code = second_pdb_code
        self.pdb_directory = pdb_directory
        self.output_file = output_file
        self.append_to_output = append_to_output

    def writeErrorAndRaise(self, errorCode, modellerError=None):
        outputList = [self.first_pdb_code, self.second_pdb_code, 'native_overlap', "Error: %s" % errorCode]
        self.resultFh.write("%s\n" % '\t'.join(str(x) for x in outputList))
        self.resultFh.close()
        
        if (modellerError):
            print "Skipping alignment between %s and %s due to Modeller error: %s" % (self.first_pdb_code, self.second_pdb_code, str(modellerError))
        raise cluster_lib.ModellerException(errorCode)
        
    def execute(self):
        appendString = 'w'
        if (self.append_to_output):
            appendString = 'a'

        self.resultFh = open(self.output_file, appendString)
        log.none()
        env = environ()
        env.io.atom_files_directory = ['.', '../atom_files', self.pdb_directory]

        [firstPdb, firstChain] = cluster_lib.readDrpCode(self.first_pdb_code)
        [secondPdb, secondChain] = cluster_lib.readDrpCode(self.second_pdb_code)

        aln = alignment(env)
        firstModel = None
        secondModel = None
        pdbFile = os.path.join(self.pdb_directory, "%s.pdb" % self.first_pdb_code)
        if (not os.path.exists(pdbFile)):
            print "Warning: did not find expected PDB file %s" % pdbFile
            self.writeErrorAndRaise("missing_pdb")
        try:
            firstModel = model(env, file=self.first_pdb_code, model_segment=('FIRST:'+firstChain, 'LAST:'+firstChain))        
        except Exception, e:
            self.writeErrorAndRaise("first_model", e)
            raise e
        aln.append_model(firstModel, atom_files=self.first_pdb_code, align_codes=firstPdb+firstChain)

        try:
            secondModel = model(env, file=self.second_pdb_code, model_segment=('FIRST:'+secondChain, 'LAST:'+secondChain))
        except Exception, e:
            self.writeErrorAndRaise("second_model", e)
            raise e

        aln.append_model(secondModel, atom_files=self.second_pdb_code, align_codes=secondPdb+secondChain)
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
    
if __name__ == '__main__':

    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--first_pdb_code", required=True, help="DRP code for first coordinate set to align (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)")
    parser.add_argument("-s", "--second_pdb_code", required=True, help="DRP code for second coordinate set to align (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)")
    parser.add_argument("-p", "--pdb_directory", required=True, help="Location of PDB files. Expected format is identical to the 'divided' PDB FTP site at ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/")
    parser.add_argument("-o", "--output_file", help="Full path of output file. If not set, output will be sent to STDOUT")
    parser.add_argument("-a", "--append_to_output", action="store_true", help="If set, output will be appended to <output_file>; if not then existing <output_file> will be overwritten")

    if (len(sys.argv) < 2):
        print "Please run with '-h' for full usage"
        sys.exit()
    config = parser.parse_args(sys.argv[1:])

    pnoa = PairwiseNativeOverlapAligner()
    pnoa.setParams(config.first_pdb_code, config.second_pdb_code, config.pdb_directory, config.output_file, config.append_to_output)
    
    pnoa.execute()
