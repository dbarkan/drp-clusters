from modeller import *
import modeller.salign
import sys
import os
import subprocess
import re
import itertools
import time
import math
import argparse
import cluster_lib

class PairwiseDisulfideAligner:
    def __init__(self, config):
        self.config = config
        
    def execute(self):

        self.prepareRun()
        
        self.prepareDrpCodes()

        self.prepareDisulfides()

        self.getBestRms()

    def prepareRun(self):
        log.none()
        self.env = environ()
        self.env.io.atom_files_directory = ['.', '../atom_files', self.config.pdb_directory]
        appendString = 'w'
        if (self.config.append_to_output):
            appendString = 'a'
        self.resultFh = open(self.config.output_file, appendString)

    def writeErrorAndRaise(self, errorCode, modellerError=None):
        outputList = [self.config.ref_pdb_code, self.config.target_pdb_code, 'cysteine_rms', "Error: %s" % errorCode]
        self.resultFh.write("%s\n" % '\t'.join(str(x) for x in outputList))
        self.resultFh.close()

        if (modellerError):
            print "Skipping alignment between %s and %s due to Modeller error: %s" % (self.config.ref_pdb_code, self.config.target_pdb_code, str(modellerError))
            
        raise cluster_lib.ModellerError(errorCode)

    def prepareDrpCodes(self):
        [self.refPdbId, self.refChainId] = cluster_lib.readDrpCode(self.config.ref_pdb_code)
        [self.targetPdbId, self.targetChainId] = cluster_lib.readDrpCode(self.config.target_pdb_code)
        if (self.targetPdbId == self.refPdbId and self.targetChainId == self.refChainId):
            self.writeErrorAndRaise("pdbs_identical")

    def prepareDisulfides(self):
        refDisulfides = self.findDisulfides(self.refPdbId, self.refChainId, 'ref')
        targetDisulfides = self.findDisulfides(self.targetPdbId, self.targetChainId, 'target') 
        self.disulfidePairs = self.makeDisulfidePairs(refDisulfides, targetDisulfides)

    def getBestRms(self):
        iteration = 1
        bestRms = 100
        bestPairList = None
        for (iteration, disulfidePair) in enumerate(self.disulfidePairs, start=1):
            tempRefPdb = "%s%s_reduced_%s" % (self.refPdbId, self.refChainId, iteration)
            tempTargetPdb = "%s%s_reduced_%s" % (self.targetPdbId, self.targetChainId, iteration)

            #make pdb with only cysteines, ordered according to the order in disulfidePairs
            self.writePdb(tempTargetPdb, self.targetPdbId, self.targetChainId, disulfidePair[1])
            self.writePdb(tempRefPdb, self.refPdbId, self.refChainId, disulfidePair[0])

            #align pdbs, superpose to find least squares distance, get the transformation that we can use to align the full proteins, and apply it
            r =   self.getTransformation(tempRefPdb, self.refChainId, tempTargetPdb, self.targetChainId)
            #applyTransformation(refPdbId, self.refChainId, self.targetPdbId, self.targetChainId, rotation, translation, env, workingDir, statsDict)
            if (r.rms < bestRms):
                bestRms = r.rms
                bestPairList = disulfidePair
            os.remove("%s.pdb" % tempRefPdb)
            os.remove("%s.pdb" % tempTargetPdb)
            iteration += 1
        outputList = [self.config.ref_pdb_code, self.config.target_pdb_code, 'cysteine_rms', bestRms, bestPairList[0], bestPairList[1]]
        self.resultFh.write("%s\n" % '\t'.join(str(x) for x in outputList))
        self.resultFh.close()

    def findDisulfides(self, pdbId, chainId, pdbType):
        pdbFile = os.path.join(self.config.pdb_directory, '%s%s.pdb' % (pdbId, chainId))

        if (not os.path.exists(pdbFile)):
            print "Warning: did not find expected PDB file %s" % pdbFile
            self.writeErrorAndRaise("missing_pdb")

        pdbFh = open(pdbFile, 'r')
        lines = pdbFh.readlines()
        ssRe = re.compile('^SSBOND\s+\d+.*?CYS\s(\w)\s+(\d+)\s+CYS\s(\w)\s+(\d+)')
        disulfideCount = 1
        allDisulfides = {}
        uniqueCysteinePositions = {}

        for line in lines:
            ssBond = ssRe.search(line)

            if (ssBond):
                firstChain = ssBond.group(1)
                firstPosition = ssBond.group(2)
                secondChain = ssBond.group(3)
                secondPosition = ssBond.group(4)

                if (firstChain == chainId and secondChain == chainId):
                    if (firstPosition in uniqueCysteinePositions or secondPosition in uniqueCysteinePositions): #don't process DRPs with cysteines involved in more than one bond
                        self.writeErrorAndRaise("%s_multibonded" % pdbType)
                    uniqueCysteinePositions[firstPosition] = 1
                    uniqueCysteinePositions[secondPosition] = 1

                    allDisulfides[disulfideCount] = {}
                    allDisulfides[disulfideCount]["source"] = firstPosition
                    allDisulfides[disulfideCount]["target"] = secondPosition
                    disulfideCount += 1

        if (len(allDisulfides) < 1):
            self.writeErrorAndRaise("%s_missing_disulfides" % pdbType)

        return allDisulfides

    def makeDisulfidePairs(self, refDisulfides, targetDisulfides):
        #each combination is a list of pairs of cysteine positions.
        #we need all possible ways to align pairs. First we make a list of possible permutations of disulfides. 
        #For each permutation, you can all cysteine pairs 1-2 or 2-1, so create all possible orderings of cysteines for each position.
        finalPairs = []  

        #The shorter list stays static while we get all possible orderings for the longer list

        shorterDisulfides = targetDisulfides
        longerDisulfides = refDisulfides

        if (len(refDisulfides.keys()) < len(targetDisulfides.keys())):
            longerDisulfides= targetDisulfides
            shorterDisulfides = refDisulfides

        longerList = sorted(longerDisulfides.keys())
        shorterList = sorted(shorterDisulfides.keys())
        shorterLength = len(shorterList)

        shorterCysteineList = []
        for position in shorterList:
            shorterCysteineList.append(shorterDisulfides[position]["source"])
            shorterCysteineList.append(shorterDisulfides[position]["target"])

        totalPairCount = 0
        permutations = list(itertools.permutations(longerList, shorterLength))  #get all permutations of positions
        for permutation in permutations:

            longerCysteineCombinations = self.makeCysteineCombinations(longerDisulfides, permutation) #get cysteine orderings (1-2 vs 2-1 for each position)

            for longerCysteineCombination in longerCysteineCombinations:
                if (len(refDisulfides.keys()) < len(targetDisulfides.keys())):
                    nextCombination = [shorterCysteineList, longerCysteineCombination]
                else:
                    nextCombination = [longerCysteineCombination, shorterCysteineList]
                finalPairs.append(nextCombination)
                totalPairCount += 1
        return finalPairs

    def makeCysteineCombinations(self, disulfides, disulfidePositionList):

        #total number of orderings is 2^x where x is number of disulfides. To create all orderings, convert number to binary
        #and for each position in binary string, output source-target or target-source depending on if it's 0 or 1. This way we get all
        #possible orderings
        positionCount = len(disulfidePositionList)
        listCount = 2 ** positionCount
        cysteineCombinations = []

        for i in range(listCount):
            nextList = []
            #get binary string. Begins with 0b so chop those two characters off
            binary = bin(i)
            binary = binary[2:]
            tempBinary = binary

            #Want all 0's to be represented to precision of length of list, so pad with 0's if they're missing
            for k in range (positionCount - len(binary)):
                tempBinary = "0" + tempBinary
            binary = tempBinary
            counter = 0
            for j in range(positionCount):
                #iterate through positions and output source-target or target-source
                binaryValue = binary[j]
                position = disulfidePositionList[counter]
                counter += 1
                sourceCysteine = disulfides[position]["source"]
                targetCysteine = disulfides[position]["target"]
                if (binaryValue == "0"):
                    nextList.append(sourceCysteine)
                    nextList.append(targetCysteine)
                else:
                    nextList.append(targetCysteine)
                    nextList.append(sourceCysteine)
            cysteineCombinations.append(nextList)

        return cysteineCombinations

    def writePdb(self, pdbOutputFile, pdbId, chainId, cysteineOrder):
        #go through PDB, find cysteines in disulfide bonds, renumber their residues according to cysteineOrder, and output only those cysteine atoms
        pdbInputFh = open(os.path.join(self.config.pdb_directory, "%s%s.pdb" % (pdbId, chainId)))
        outputFh = open("%s.pdb" % pdbOutputFile, 'w')

        cysteineOrderDict = {}
        outputLines = {}
        cysteineCounter = 1

        #make map of original postion of cysteines to the final position we want them in the pdb file 
        #(this postion is the same as their index in cysteineOrder)
        for position in cysteineOrder:
            cysteineOrderDict[position] = cysteineCounter
            outputLines[cysteineCounter] = []
            cysteineCounter += 1

        pdbRe = re.compile('(^ATOM.*CYS\s)(\w)(\s*)(\d+)(.*)')
        terRe = re.compile('^TER')
        pdbLines = pdbInputFh.readlines()
        foundChain = 0
        for line in pdbLines:
            terLine = terRe.search(line)
            if (terLine and foundChain == 1):  #account for multiple copies, i.e. NMR
                break
            cysLine = pdbRe.search(line)
            if (cysLine):  #found a cysteine atom with the same chain
                firstPart = cysLine.group(1)
                lineChain = cysLine.group(2)
                space = cysLine.group(3)
                resNumber = cysLine.group(4)
                lastPart = cysLine.group(5)
                if (lineChain == chainId and cysteineOrderDict.has_key(resNumber)):  #in a disulfide bond
                    foundChain = 1
                    finalCysPosition = cysteineOrderDict[resNumber]
                    outputString = "%s%s%s%s%s" % (firstPart, lineChain, space, finalCysPosition, lastPart)
                    outputLines[finalCysPosition].append(outputString)   #add line to list of those to be written for this cysteine

        #check we found all cysteines
        for cysPosition in cysteineOrderDict:
            if (not cysteineOrderDict[cysPosition] in outputLines):
                self.writeErrorAndRaise("cysteine_not_found")

        #check we have same number as expected
        if (len(outputLines) != len(cysteineOrderDict.keys())):
            self.writeErrorAndRaise("cysteine_count_mismatch")

        #output
        for finalCysPosition in sorted(outputLines.keys()):
            outputList = outputLines[finalCysPosition]
            for outputString in outputList:
                outputFh.write("%s\n" % outputString)

    def getTransformation(self, reducedRefPdb, refChain, reducedTargetPdb, targetChain):

        #get the rotation and translation for superposing the given target PDB cysteine combination onto the ref pdb
        #PDBs are only cysteines ('reduced' files)
        aln = alignment(self.env)
        r = None
        try:
            reducedRefModel = model(self.env, file=reducedRefPdb, model_segment=('FIRST:'+refChain, 'LAST:'+refChain))
            aln.append_model(reducedRefModel, atom_files=reducedRefPdb, align_codes=reducedRefPdb+refChain)
            reducedTargetModel = model(self.env, file=reducedTargetPdb, model_segment=('FIRST:'+targetChain, 'LAST:'+targetChain))
            aln.append_model(reducedTargetModel, atom_files=reducedTargetPdb, align_codes=reducedTargetPdb+targetChain)

            aln.salign()

            atmsel = selection(reducedRefModel)
            r = atmsel.superpose(reducedTargetModel, aln)
        except Exception as e:
            self.writeErrorAndRaise('modeller', e)
        trans = r.translation
        rotation = r.rotation

        #get sum of distances across equivalent ca's and sg's in system
        return r


def getParser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Structurally Align two DRPs by closest matching set of disulfide bonds")

    parser.add_argument("-f", dest="ref_pdb_code", metavar="<string>", required=True,
                        help="DRP code for first coordinate set to align (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)\n\n")
    
    parser.add_argument("-s", dest="target_pdb_code", metavar="<string>", required=True,
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

    pnoa = PairwiseDisulfideAligner(config)
    
    pnoa.execute()
