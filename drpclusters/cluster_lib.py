import sys
import os
import re
import subprocess
import logging
from operator import itemgetter
import collections
import math
import time
import shutil

logger = logging.getLogger(__name__)
class DrpClusterException(Exception):
    pass

class ModellerException(DrpClusterException):
    pass

class InvalidDrpCodeError(DrpClusterException):
    pass

SingletonTuple = collections.namedtuple("SingletonTuple", ['singletonDrpCode', 'refDrp',
                                                           'distance', 'singletonType'])
def readDrpCode(drpCode):
    """Return tuple consisting of PDB ID and Chain ID given drpCode.

    param drpCode: 5-character drpCode
    raises InvalidDrpCodeError if drpCode is not expected format
    """
    if (len(drpCode) != 5 and not drpCode.endswith('_fit')):
        raise InvalidDrpCodeError("Expected 5 character DRP code, format PDB ID + Chain ID "
                                  "(e.g. 1zdcA). Instead got DRP code %s." % drpCode)

    #These methods are a little simple but allow for modifying drp code format in a single location
    pdbId = drpCode[0:4]
    chainId = drpCode[4]
    return (pdbId, chainId)

def makeDrpCode(pdbId, chainId):
    """Return DRP code to identify a DRP given its PDB and Chain IDs"""
    #exception if not 5 chars?
    return "%s%s" % (pdbId, chainId)

def readDrpCodeFile(fileName):
    pfr = PtFileReader(fileName)
    drpCodeList = []
    for line in pfr.getLines():
        readDrpCode(line) #just validate format
        drpCodeList.append(line)
    return drpCodeList
        
class Drp(object):

    """Class to handle all aspects of DRPs"""

    def __init__(self, drpCode):
        [self.pdbId, self.chainId] = readDrpCode(drpCode)
        self.drpCode = drpCode
        self.model = None
        self.completeModel = None
        self.singletonTuple = None
        self.isCentroid = False

    def getPdbId(self):
        return self.pdbId

    def getChainId(self):
        return self.chainId

    def getDrpCode(self):
        return makeDrpCode(self.pdbId, self.chainId)
    
    def setCentroid(self):
        self.isCentroid = True
        
    def getScopFold(self):
        return self.scopFold

    def setScopFold(self, scopFold):
        self.scopFold = scopFold

    def getClusterMemberId(self):
        return self.drpCode

    def isSingleton(self):
        return self.singletonTuple is not None

    def setSingletonTuple(self, st):
        self.singletonTuple = st

    def getSingletonTuple(self):
        return self.singletonTuple

    def getModel(self, pdbDirectory):
        if (not self.model):
            from modeller import *
            env = environ()
            env.io.atom_files_directory = [pdbDirectory]
            self.model = model(env, file=self.drpCode, model_segment=('FIRST:'+self.chainId,
                                                                      'LAST:'+self.chainId))
        return self.model

    def setDistanceToCentroid(self, distance):
        self.distanceToCentroid = distance

    def getDistanceToCentroid(self):
        return self.distanceToCentroid

class Runner(object):  
    def __init__(self, config):
        self.config = config
        self.makeOutputDirectory()
        self.initLogging()

    def makeOutputDirectory(self):
        runDirectory = self.config.run_directory
        self.fullOutputDir = runDirectory
        if (not os.path.exists(self.fullOutputDir)):
            os.mkdir(self.fullOutputDir)

    def getFullOutputFile(self, fileName):
        return os.path.join(self.fullOutputDir, fileName)

    def initLogging(self, level=logging.INFO):
        formatString =  '[%(asctime)s] [%(name)s] [%(levelname)s]: %(message)s'
        logging.basicConfig(filename=self.makeLogFileName(), filemode='w', format=formatString, level=level)
        consoleHandler = logging.StreamHandler(sys.stdout)  
        consoleHandler.setLevel(level)
        formatter = logging.Formatter(formatString, "%Y-%m-%d %H:%M:%S")
        consoleHandler.setFormatter(formatter)
        logging.getLogger().addHandler(consoleHandler)
        logging.getLogger().setLevel(level)

    def makeLogFileName(self):
        return self.getFullOutputFile("%s.log" % self.makeLogPrefix())

    
class PtFileReader:

    """Essentially python file object that provides commonly used processing functionality"""

    def __init__(self, fileName, skipBlanks=True, skipHash=True):
        """Read input file and save all lines after processing

        Similar to python File object, except we commonly will want to skip blank lines in the file
        (programs often break when getting an unexpected blank line) and also skip lines beginning
        with '#' as these are often comments as opposed to data.

        @param fileName full path file name to be processed
        @param skipBlanks set to true to ignore blank lines
        @param skipHash set to true to ignore lines beginning with '#'
        """
        fh = open(fileName, 'r')
        lines = fh.readlines()
        finalLines = []
        blankRe = re.compile('^\s*$')
        hashRe = re.compile('^\#')

        for line in lines:
            line = line.rstrip('\r\n')
            blankLine = blankRe.search(line)
            if (blankLine and skipBlanks):
                continue
            hashLine = hashRe.search(line)
            if (hashLine and skipHash):
                continue
            finalLines.append(line)
        self.finalLines = finalLines
        self.fileName = fileName
        fh.close()
        
    def getLines(self):
        return self.finalLines

class PtpShutils:
    def copyFile(self, sourceDir, sourceFile, destinationDir, destinationFile=None):
        if (destinationFile is None):
            destinationFile = sourceFile
        destinationFile = os.path.join(destinationDir, destinationFile)
        """Copy source file from sourceDir to destinationDir; don't return until copy is done"""
        self.tryShutil(shutil.copy, os.path.join(sourceDir, sourceFile), destinationFile)

        #possible error: if file is large, conceivably this could return before copy is done
        self.sleepUntilDone(destinationFile, predicate=self.fileDoesNotExist)
    
    def tryShutil(self, function, *args):
        """Run a python shutil module command"""
        try:
            function(*args)
        except IOError as e:
            raise PtpShutilError(e, function, args)

    def moveFile(self, sourceDir, sourceFile, destinationDir, destinationFile=None):
        """Move source file from sourceDir to destinationDir; don't return until move is done"""
        fullSourceFile = os.path.join(sourceDir, sourceFile)
        if (destinationFile is None):
            destinationFile = sourceFile
        fullDestinationFile = os.path.join(destinationDir, destinationFile)
        
        self.tryShutil(shutil.move, fullSourceFile, fullDestinationFile)

        self.sleepUntilDone(fullDestinationFile, predicate=self.fileDoesNotExist)

    def runSubprocess(self, args, checkStdError=True):
        """Run python subprocess module command; by default, raise exception if anything was written to stderr"""
        process = subprocess.Popen(args, shell=False, stderr=subprocess.PIPE)
        processOutput = process.communicate()
        if (processOutput[1] != "" and checkStdError):
            raise Exception("Got subprocess error.\nRan method args %s\nGot stderr %s" % (args, processOutput[1]))
        return processOutput

    def sleepUntilDone(self, fileName, predicate):
        """Sleep until predicate involving fileName is true; useful for avoiding race conditions in file manipulation"""
        sleepTime = 0
        while (predicate(fileName)):
            print "sleep 1 second"
            time.sleep(1)
            sleepTime += 1
            if (sleepTime > 10):
                raise Exception("Timeout on file %s" % fileName)

    def fileDoesNotExist(self, fileName):
        """Simple predicate that returns true if fileName doesn't exist in file system"""
        return not os.path.exists(fileName)

    def fileExists(self, fileName):
        """Simple predicate that returns true if fileName exists in file system"""
        return os.path.exists(fileName)

class PtpShutilError(Exception):
    def __init__(self,  e, function, *args):
        print "got error %s when running function %s args %s" % (str(e), function, str(args))

class ScopInfo:
    def __init__(self, fileName):
        self.noScopInfoKeyword = "Unannotated"

        self.pdbIdsToScopIds = {}
        self.scopFoldToName = {}
        self.scopFoldToName[self.noScopInfoKeyword] = self.noScopInfoKeyword
        reader = PtFileReader(fileName)
        lines = reader.getLines()
        for line in lines:
            cols = line.split('\t')
            scopType = cols[1]
            if (scopType == "px"):
                self.processScopCodes(cols)
            if (scopType == "cf"):
                self.processScopFoldName(cols)

    def processScopFoldName(self, cols):
        self.scopFoldToName[cols[2]] = cols[4]

    def getScopFoldName(self, scopFoldId):
        return self.scopFoldToName[scopFoldId]

    def processScopCodes(self, cols):
        scopFamily = cols[2]
        pdbInfo = cols[4]
        [pdbId, chainInfo] = pdbInfo.split(' ')
        chainList = chainInfo.split(',')
        for chain in chainList:
            singleChainInfo = chain.split(':')
            chainId = singleChainInfo[0]
            drpCode = makeDrpCode(pdbId, chainId)
            self.pdbIdsToScopIds[drpCode] = scopFamily

    def getScopFamilyId(self, drpCode):
        if (drpCode not in self.pdbIdsToScopIds):
            return self.noScopInfoKeyword
        return self.pdbIdsToScopIds[drpCode]

    def getScopFoldId(self, drpCode):
        familyId = self.getScopFamilyId(drpCode)
        if (familyId == self.noScopInfoKeyword):
            return familyId
        foldCols = familyId.split(".")
        return ".".join(foldCols[0:2])

class DisulfideDrpDistances:
      
    def stringToList(self, inputString):
        inputString = inputString.lstrip("'[")
        inputString = inputString.rstrip("]'")
        inputList = inputString.split("', '")

        return inputList

    def readFile(self, fileName, filterDrps=[]):
        self.drpDistanceTuple = collections.namedtuple("drpDistanceTuple", ['distance', 'firstCysOrder', 'secondCysOrder'])
        self.distances = {}

        reader = PtFileReader(fileName)
      
        lines = reader.getLines()
        for (lineCount, line) in enumerate(lines):

            line = line.rstrip("\n\r")
            if ("Error" in line):
                continue
            [firstPdbCode, secondPdbCode, distanceType, distance, firstCysteineString, secondCysteineString] = line.split('\t')

            #2ny9X   2frbA   cysteiene_rms   2.58908390999   ['35', '16', '3', '30'] ['2', '7', '3', '13']
            firstCysteineOrder = self.stringToList(firstCysteineString)
            secondCysteineOrder = self.stringToList(firstCysteineString)

            self.addDistance(firstPdbCode, secondPdbCode, distance, firstCysteineOrder, secondCysteineOrder)
            self.addDistance(secondPdbCode, firstPdbCode, distance, secondCysteineOrder, firstCysteineOrder)

    def addDistance(self, firstCode, secondCode, distance, firstOrder, secondOrder):
        if (firstCode not in self.distances):
            self.distances[firstCode] = {}
        self.distances[firstCode][secondCode] = self.drpDistanceTuple(distance, firstOrder, secondOrder)
    
    def getDisulfideDistanceTuple(self, firstCode, secondCode):
        if (firstCode not in self.distances):
            print "didn't get first pdb code %s in distances" % firstCode
            return None
        if (secondCode not in self.distances[firstCode]):
            print "didn't get distance betweenn first pdb code %s and second %s" % (firstCode, secondCode)
            return None
        return self.distances[firstCode][secondCode]
            

class DistanceMatrix:

    def __init__(self):
        self.matrix = {}

    def getDistance(self, firstId, secondId):
        if (firstId not in self.matrix):
            return None
        elif (secondId not in self.matrix[firstId]):
            return None
        
        return self.matrix[firstId][secondId]

    def addDistancePair(self, firstId, secondId, distance):
        if (firstId not in self.matrix):
            self.matrix[firstId] = {}
        self.matrix[firstId][secondId] = float(distance)

    def containsId(self, id):
        return id in self.matrix
    
class DisulfideDrpDistanceMatrix(DisulfideDrpDistances, DistanceMatrix):
    def addDistance(self, firstCode, secondCode, distance, firstOrder, secondOrder):
        self.addDistancePair(firstCode, secondCode, distance)
        self.addDistancePair(secondCode, firstCode, distance)

class DrpDistanceMatrix(DistanceMatrix):
    def readFile(self, fileName):
        reader = PtFileReader(fileName)
        lines = reader.getLines()
        self.distances = {}

        for (lineCount, line) in enumerate(lines):

            line = line.rstrip("\n\r")
            [firstPdbCode, secondPdbCode, distanceType, distance] = line.split('\t')
            self.addDistancePair(firstPdbCode, secondPdbCode, distance)
            self.addDistancePair(secondPdbCode, firstPdbCode, distance)

    def createDrpList(self):
        drpList = []
        for drpCode in self.matrix:
            nextDrp = Drp(drpCode)
            drpList.append(nextDrp)
        return drpList

class FastaWriter:
    def __init__(self, fileName):
        self.fh = open(fileName, 'w')

    def writeFastaLine(self, id, sequence):
        self.fh.write(">%s\n" % id)
        self.fh.write("%s\n\n" % sequence)

class ClusterSingletons:
    def __init__(self):
        self.singletonDict = {}

    def readShorterSingletonFile(self, fileName, distanceCutoff):
        self.readSingletonFile(fileName, distanceCutoff, "shorter")

    def readLongerSingletonFile(self, fileName, distanceCutoff):
        self.readSingletonFile(fileName, distanceCutoff, "longer")

    def isSingleton(self, drpCode):
        return drpCode in self.singletonDict

    def readSingletonFile(self, fileName, distanceCutoff, singletonType):
        reader = PtFileReader(fileName)
        lines = reader.getLines()
        for line in lines:
            [singletonDrp, refDrp, distance] = line.split('\t')
            if (float(distance) >= float(distanceCutoff)):
                print 'cluster singleton set singleton %s ref %s distance %s type %s' % (singletonDrp, refDrp, distance, singletonType)
                self.singletonDict[singletonDrp] = SingletonTuple(singletonDrp, refDrp, distance, singletonType)
    
    def getSingletonTuple(self, drpCode):
        return self.singletonDict[drpCode]
        
                          

class ClusterError(Exception):
    def __init__(self, msg):
        self.msg = msg


class DrpClusterSet:

    """Class for handling a group of clusters; mostly easy storage"""

    def __init__(self, clusterIndexList=None):
        """@param clusterIndexList: optional list of cluster indices that will be read in from the 
        cluster member file (others will be discarded)
        """
        self.clusterDict = {}
        finalClusterIndexList = None
        if (clusterIndexList is not None):
            finalClusterIndexList = []
            for clusterIndex in clusterIndexList:
                finalClusterIndexList.append(int(clusterIndex) - 1)
        self.clusterIndexList = finalClusterIndexList

    def getClusterIndexList(self):
        return self.clusterDict.keys()

    def addCluster(self, cluster):
        #todo -- check to make sure don't already have cluster with this index
        self.clusterDict[int(cluster.getClusterIndex())] = cluster

    def getDrp(self, drpCode):
        """Search all my clusters for the drp represented by <drpCode>"""
        for cluster in self.getClusterList():
            if (cluster.hasDrp(drpCode)):
                return cluster.getDrp(drpCode)
        return None

    def getClusterIndexForDrpCode(self, drpCode):
        for cluster in self.getClusterList():
            if (drpCode in cluster.getDrpCodeList()):
                return cluster.getClusterIndex()
        return None
        
    def renumberBySize(self):
        """Reassign indices for all clusters according to their size, starting with 0 for the largest """
        for cluster in self.getClusterList():
            if (len(cluster.getMemberList()) == 0):
                del self.clusterDict[cluster.getClusterIndex()]
        
        tempClusterList = []
        for (i, cluster) in enumerate(sorted(self.getClusterList(), key=lambda x: x.getSize(), reverse=True)):
            del self.clusterDict[cluster.getClusterIndex()]
            cluster.updateClusterIndex(i)
            tempClusterList.append(cluster)

        for cluster in tempClusterList:
            self.addCluster(cluster)

    def readSingletonFiles(self, longerFile, shorterFile, cutoff):
        """Assign singleton information to DRPs in my clusters."""
        clusterSingletons = ClusterSingletons()
        clusterSingletons.readLongerSingletonFile(longerFile, cutoff)
        clusterSingletons.readShorterSingletonFile(shorterFile, cutoff)

        for cluster in self.getClusterList():
            for drp in cluster.getMemberList():
                if (clusterSingletons.isSingleton(drp.drpCode)):
                    drp.setSingletonTuple(clusterSingletons.getSingletonTuple(drp.drpCode))

    def getClusterList(self):
        return self.clusterDict.values()

    def getClusterCount(self):
        return len(self.clusterDict)

    def getCluster(self, clusterIndex):
        return self.clusterDict[int(clusterIndex)]

    def writeSingletonFile(self, fileName, singletonType):
        outputFh = open(fileName, 'w')
        for cluster in self.clusterDict.values():
            
            for drp in cluster.getMemberList():
                #print 'next %s drp test: %s' % (singletonType, drp.drpCode)
                if (drp.isSingleton() and drp.getSingletonTuple().singletonType == singletonType):  #1/20/16: fixed bug that wrote all singletons regardless of type
                    outputLine = [drp.drpCode, drp.getSingletonTuple().refDrp, drp.getSingletonTuple().distance]
                    outputFh.write("%s\n" % '\t'.join(str(x) for x in outputLine))
            
    def writeClusterMemberFile(self, fileName):
        outputFh = open(fileName, 'w')
        if (sorted(self.clusterDict.keys())[0] != 0):
            raise ClusterError("Writing clusters; first cluster index should be 0")

        for clusterIndex, cluster in sorted(self.clusterDict.iteritems(), key=itemgetter(0)):
            for drp in cluster.getMemberList():
                oneBasedClusterIndex = cluster.getClusterIndex() + 1
                outputFh.write("%s\t%s\n" % (oneBasedClusterIndex, drp.drpCode))

    def hasClusterIndex(self, clusterIndex):
        return clusterIndex in self.clusterDict

    def readClusterMemberFile(self, fileName, skipList=None):
        """Read a previously generated cluster member file. If self.clusterIndexList was populated upon
        initialization, only include those clusters.

        @param skipList -- List of DRP codes that will be skipped (sort of a hack to deal with problematic DRPs)
        """

        print "reading cluster members"
        reader = PtFileReader(fileName)
        lines = reader.getLines()
        for line in lines:
            [clusterIndex, drpCode] = line.split('\t')
            if (clusterIndex == "0"):
                raise ClusterError("Read clusters in from file %s; have unexpected cluster index of 0 (lowest index should be 1)")
            clusterIndex = int(clusterIndex) - 1
            if (skipList is not None and drpCode in skipList):
                continue
            if (self.includeCluster(clusterIndex)):
                clusterIndex = int(clusterIndex)
                if (clusterIndex not in self.clusterDict):
                    self.clusterDict[clusterIndex] = DrpCluster(clusterIndex)
                
                self.clusterDict[clusterIndex].addMember(Drp(drpCode))

    def includeCluster(self, clusterIndex):
        if (self.clusterIndexList is not None):
            return clusterIndex in self.clusterIndexList
        return True

    def getFirstDrpFromAllClusters(self):
        """Get list of representatives for clusters. 

        Each rep is the first DRP in the cluster, sorted alphabetically by PDB Code
        """
        sortedClusterRepList = []
        for drpCluster in self.clusterDict.values():
            sortedClusterRepList.append(drpCluster.getFirstDrp())
        return sortedClusterRepList

class Cluster:
    def __init__(self, clusterId):
        self.clusterId = int(clusterId)
        self.memberList = []
        self.isActive = True
        self.centroid = None

    def getClusterIndex(self):
        return self.clusterId

    def updateClusterIndex(self, i):
        self.clusterId = i

    def setClusterIndex(self, clusterIndex):
        self.clusterId = int(clusterIndex)

    def addMember(self, clusterMember):
        self.memberList.append(clusterMember)

    def inactivate(self):
        self.isActive = False

    def getMemberList(self):
        return self.memberList

    def calculateCentroidDistances(self, distanceMatrix):
        for member in self.memberList:
            distance = 0
            if (member.getClusterMemberId() != self.centroid.getClusterMemberId()):
                distance = distanceMatrix.getDistance(member.getClusterMemberId(), self.centroid.getClusterMemberId())
            member.setCentroidDistance(distance)

    def setCentroid(self, centroid):
        self.centroid = centroid

    def getCentroid(self):
        return self.centroid

    def calculateCentroid(self, distanceMatrix, distanceComparer):
        bestDistance = distanceComparer.getInitialDistance()
        self.centroid = None
        if (len(self.memberList) == 1):
            self.centroid = self.memberList[0]
        else:
            for firstMember in self.memberList:
                averageDistance = self.getAverageDistanceForMember(firstMember, distanceMatrix, distanceComparer)
                print "cluster %s average disance for drp %s is %s" % (self.getClusterIndex(), firstMember.getClusterMemberId(), averageDistance)
                if (distanceComparer.betterDistance(averageDistance, bestDistance)):
                    self.centroid = firstMember
                    bestDistance = averageDistance
        self.bestDistance = bestDistance
        self.centroid.setCentroid()

    def getCentroidDistance(self):
        return self.bestDistance

    def getAverageDistanceForMember(self, member, distanceMatrix, distanceComparer):
        totalDistance = 0
        count = 0
        averageDistance = distanceComparer.getInitialDistance()
        for secondMember in self.memberList:
            if (secondMember.getClusterMemberId() != member.getClusterMemberId()):
                nextDistance = distanceMatrix.getDistance(member.getClusterMemberId(), secondMember.getClusterMemberId())
                if (member.getClusterMemberId() == "2crdA"):
                    print "2crda with %s is %s" % (secondMember.getClusterMemberId(), nextDistance)
                if (nextDistance is not None):
                    totalDistance += nextDistance
                    count += 1
        if (count > 0):
            averageDistance = float(totalDistance) / float(count)
            if (member.getClusterMemberId() == "2crdA"):
                print "total distance %s count %s" % (totalDistance, count)
        return averageDistance

    def getSortedFirstMember(self):
        sortedList = sorted(self.memberList, key=lambda x: x.getClusterMemberId())
        return sortedList[0]

    def getSize(self):
        return len(self.memberList)

    
class DrpCluster(Cluster):

    """Store information about a group of DRPs """ 
   
    def getSingletonDrpList(self):
        """Get the list of all singleton DRPs in this cluster"""
        singletonDrpList = []
        for drp in self.memberList:
            if (drp.isSingleton()):
                singletonDrpList.append(drp)

        return singletonDrpList
        
    def getDrpCodeList(self):
        """Return the list of strings representing DRP codes for all DRPs in the cluster"""

        drpCodeList = []
        for drp in self.memberList:
            drpCodeList.append(drp.drpCode)
        return drpCodeList

    def getFirstDrp(self):  
        """Return the first DRP in the cluster (sort order is alphabetical by PDB code"""
        return sorted(self.getDrpList())[0]

    def hasDrp(self, drpCode):
        return drpCode in self.getDrpCodeList()

    def getDrp(self, drpCode):
        for drp in self.memberList:
            if (drp.drpCode == drpCode):
                return drp

    def salignDrps(self, pdbDirectory, outputDir, writeFit, skipSingletons=True, manualSingletonAlignment=True):
        """Run the SALIGN Program which performs an iterative structure / sequence alignment over all DRPs in this cluster.
        Must have MODELLER installed for this to function
        """
        from modeller import *
        import modeller
        from modeller import scripts

        import salign_tools
        env = environ()
        env.io.atom_files_directory = [pdbDirectory]
        aln = alignment(env)
        logger.info("Running SALIGN on Cluster %s" % int(self.getClusterIndex() + 1))
        for drp in self.getMemberList():
            if (skipSingletons and drp.isSingleton()):
                continue
            aln.append_model(drp.getModel(pdbDirectory), atom_files=drp.drpCode, align_codes=drp.drpCode)


        saveStdout = sys.stdout
        sys.stdout = open(os.devnull, "w")
        
        salign_tools.myIterativeStructuralAlign(aln, True) #keeping writeFit as True for now to avoid odd SALIGN error
        sys.stdout = saveStdout

        alignmentFastaFile = os.path.join(pdbDirectory, "%s_alignment.fasta" % self.getClusterIndex())
        
        aln.write(alignmentFastaFile, alignment_format='FASTA')
        self.aln = aln

        if (manualSingletonAlignment):
            #hack -- need to align each singleton to its reference DRP. Reference drps don't have their model
            #coordinates changed, only the _fit file that was written during SALIGN (which is in CWD). Create
            #new drp object for the fit fil and align the singleton to it. Even hackier as the atom files for 
            #each DRP are in two different directories
            env.io.atom_files_directory = [pdbDirectory, os.getcwd()]
            for singletonDrp in self.getSingletonDrpList():
                singletonAln = alignment(env)

                if (not os.path.exists(singletonDrp.getSingletonTuple().refDrp + "_fit.pdb")):
                    print "singleton %s did not get ref drp %s"  % (singletonDrp.drpCode, singletonDrp.getSingletonTuple().refDrp + "_fit")

                fitRefDrp = Drp(singletonDrp.getSingletonTuple().refDrp + "_fit")

                fitRefDrp.getModel(".")
                singletonDrp.getModel(pdbDirectory)

                singletonAln.append_model(fitRefDrp.getModel(pdbDirectory), atom_files=fitRefDrp.drpCode, align_codes=fitRefDrp.drpCode)
                singletonAln.append_model(singletonDrp.getModel(pdbDirectory), atom_files=singletonDrp.drpCode, align_codes=singletonDrp.drpCode)
                saveStdout = sys.stdout
                sys.stdout = open(os.devnull, "w")
                modeller.salign.iterative_structural_align(singletonAln)

                atmsel = selection(fitRefDrp.getModel(pdbDirectory)).only_atom_types('CA')
                    

                r = atmsel.superpose(singletonDrp.getModel(pdbDirectory), singletonAln)
                sys.stdout = saveStdout                                     
                singletonDrp.getModel(pdbDirectory).write("%s_fit.pdb" % singletonDrp.drpCode)

class SingletonProcessor:
    
    """
    Class that processes singletons -- takes a set of 'reference' clusters all containing a large number of 
    DRPs and a set of target clusters each containing a small number of DRPs, and for each target DRP, assigns
    it to the cluster containing the reference DRP with the highest similarity to that target DRP, above a cutoff.
    This process consolidates clusters a bit in cases where the target DRP is more similar to other DRPs in its 
    target cluster but still could go into the larger reference DRP. The distance can use the length of the longer DRP
    as the reference (better quality) or the length of the shorter DRP as the reference (not as good quality).

    @param topClusterOrder -- DRPs from reference clusters (note these are not the reference clusters themselves)
    @param bottomClusterSet -- set of target clusters
    """

    def __init__(self, distanceMatrix, topClusterOrder, bottomClusterSet, singletonCutoff, drpToLength=None):
        self.ldc = LargerDistanceComparer(singletonCutoff)
        self.topClusterOrder = topClusterOrder
        self.bottomClusterSet = bottomClusterSet
        self.distanceMatrix = distanceMatrix
        self.drpToLengthMap = drpToLength

    def processSingletons(self):
        """Top level method for processing"""
        for cluster in self.bottomClusterSet.getClusterList():
            for drp in cluster.getMemberList():
                self.addSingletonInfo(drp)

    def addSingletonInfo(self, targetDrp):
        """Get the closest reference DRP for this target DRP and add that information to the target DRP"""
        refDrp = self.getBestRefDrp(targetDrp)
        
        if (refDrp is not None):
            bestDistance = self.distanceMatrix.getDistance(targetDrp.drpCode, refDrp.drpCode)

            singletonTuple = SingletonTuple(targetDrp.drpCode, refDrp.drpCode, bestDistance, self.getSingletonTypeString())
            targetDrp.setSingletonTuple(singletonTuple)

    def getBestRefDrp(self, targetDrp):
        """Get closest reference DRP to the target DRP that is above self.singletonCutoff"""
        refDrp = None
        bestDistance = self.ldc.getInitialDistance()
        for topDrp in self.topClusterOrder:
            distance = self.distanceMatrix.getDistance(targetDrp.drpCode, topDrp.drpCode)
            if (self.passesAllCutoffs(distance, bestDistance, targetDrp, topDrp)):
                bestDistance = distance
                refDrp = topDrp
        return refDrp

    def passesAllCutoffs(self, distance, bestDistance, targetDrp, topDrp):
        """Ensure that the thresholds for <targetDrp> to be considered a singleton for <topDrp> are met"""
        if (self.ldc.betterDistance(distance, bestDistance) and 
            self.ldc.strictPassesCutoff(distance) and 
            self.passesSingletonType(targetDrp, topDrp)):
            
            return True
        return False

class LongerFractionSingletonProcessor(SingletonProcessor):
    """
    Subclass for processing whether a target DRP passes the similarity cutoff to a reference DRP, using the length
    of the longer DRP as the reference in the similarity calculation
    """  
    def getSingletonTypeString(self):
        return "longer"

    def passesSingletonType(self, targetDrp, refDrp):
        return True

class ShorterFractionSingletonProcessor(SingletonProcessor):
    """
    Subclass for processing whether a target DRP passes the similarity cutoff to a reference DRP, using the length
    of the shorter DRP as the reference in the similarity calculation and also ensuring that the target DRP is smaller
    than the reference DRP to ensure the calculation makes sense
    """
    def getSingletonTypeString(self):
        return "shorter"
    
    def passesSingletonType(self, targetDrp, refDrp):
        targetLength = self.drpToLengthMap[targetDrp.getDrpCode()]
        refLength = self.drpToLengthMap[refDrp.getDrpCode()]        
        if (targetLength > refLength):
            return False
        return True

class DistanceComparer:
    def __init__(self, cutoff):
        self.cutoff = float(cutoff)

class SmallerDistanceComparer(DistanceComparer):
    def betterDistance(self, firstDistance, secondDistance):
        return firstDistance < secondDistance

    def betterOrEqualDistance(self, firstDistance, secondDistance):
        return firstDistance <= secondDistance

    def getInitialDistance(self):
        return 1000000000

    #a little asymmetric since <= for cutoff but < for betterDistance
    def passesCutoff(self, distance):
        return distance <= self.cutoff

    def strictPassesCutoff(self, distance):
        return distance < self.cutoff

    def getNullCrossProductValue(self):
        return 1000000000

class LargerDistanceComparer(DistanceComparer):
    def betterDistance(self, firstDistance, secondDistance):
        return firstDistance > secondDistance

    def betterOrEqualDistance(self, firstDistance, secondDistance):
        return firstDistance >= secondDistance

    def getInitialDistance(self):
        return -1000000000

    def passesCutoff(self, distance):
        return distance >= self.cutoff

    def strictPassesCutoff(self, distance):
        return distance > self.cutoff

    def getNullCrossProductValue(self):
        return 0

        
class HierarchicalClusters:
    def __init__(self, distanceMatrix, clusterList, distanceComparer):
        self.clusterDistanceMatrix = {}
        self.distanceMatrix = distanceMatrix
        self.distanceComparer = distanceComparer
        self.clusterList = clusterList
        self.initializeClusterDistances()
        self.bestDistance = None
        
    def initializeClusterDistances(self):
        for i in range (len(self.clusterList)):
            self.clusterDistanceMatrix[i] = {}
        for i in range (len(self.clusterList)):
            self.clusterDistanceMatrix[i] = {}
            
            for j in range (i+1, len(self.clusterList)):
                self.addInitialDistancePair(i, j)
                self.addInitialDistancePair(j, i)

    def addInitialDistancePair(self, firstClusterIndex, secondClusterIndex):
        firstMember = self.clusterList[firstClusterIndex].getMemberList()[0]
        secondMember = self.clusterList[secondClusterIndex].getMemberList()[0]
        distance = self.distanceMatrix.getDistance(firstMember.getClusterMemberId(), secondMember.getClusterMemberId())

        self.clusterDistanceMatrix[firstClusterIndex][secondClusterIndex] = distance

    def clusterHierarchical(self):
        while (self.keepClustering()):
            
            [bestFirstClusterIndex, bestSecondClusterIndex, bestDistance] = self.getBestDistance()

            if (not self.distanceComparer.passesCutoff(bestDistance)):
                break

            self.mergeClusters(bestFirstClusterIndex, bestSecondClusterIndex)
            
            self.clusterList[bestSecondClusterIndex].inactivate()

            self.updateDistances(bestFirstClusterIndex)

    def updateDistances(self, clusterIndex):
        secondClusterIndices = self.clusterDistanceMatrix[clusterIndex]
        for secondClusterIndex in secondClusterIndices:
            self.updateAverageLinkageDistances(clusterIndex, secondClusterIndex)

    def updateAverageLinkageDistances(self, firstClusterIndex, secondClusterIndex):
        firstCluster = self.clusterList[firstClusterIndex]
        secondCluster = self.clusterList[secondClusterIndex]
        if (secondCluster.isActive):
            allDistances = self.getDistanceCrossProduct(firstCluster, secondCluster)
            average = self.distanceComparer.getNullCrossProductValue()
            if (len(allDistances) > 0):
                average = float(sum(allDistances)) / float(len(allDistances))
            self.clusterDistanceMatrix[firstClusterIndex][secondClusterIndex] = average 
            self.clusterDistanceMatrix[secondClusterIndex][firstClusterIndex] = average

    def getDistanceCrossProduct(self, firstCluster, secondCluster):
        
        allDistances = []
        for firstMember in firstCluster.getMemberList():
            for secondMember in secondCluster.getMemberList():
                distance = self.distanceMatrix.getDistance(firstMember.getClusterMemberId(), secondMember.getClusterMemberId())
                if (distance is not None):
                    allDistances.append(distance)
        return allDistances

    def mergeClusters(self, targetClusterIndex, fromClusterIndex):
        for member in self.clusterList[fromClusterIndex].getMemberList():
            self.clusterList[targetClusterIndex].addMember(member)
        
    def keepClustering(self):
        activeCount = 0
        for cluster in self.clusterList:
            if (cluster.isActive):
                activeCount += 1
        if (activeCount % 10 == 0):
            distanceOutput = None
            if (self.bestDistance):
                distanceOutput = round(self.bestDistance, 4)
            else:
                distanceOutput = "N/A"
                
            logger.info("Number of active clusters: %s next distance to merge: %s" % (activeCount, distanceOutput))
            
        return activeCount > 1

    def getBestDistance(self):
        bestDistance = self.distanceComparer.getInitialDistance()
        
        bestFirstIndex = None
        bestSecondIndex = None
        
        for (firstClusterIndex, secondClusterIndices) in sorted(self.clusterDistanceMatrix.iteritems(), key=itemgetter(0)):
            for secondClusterIndex in sorted(secondClusterIndices):

                if (self.betterDistance(firstClusterIndex, secondClusterIndex, bestDistance)):
                    bestDistance = self.clusterDistanceMatrix[firstClusterIndex][secondClusterIndex]
                    bestFirstIndex = firstClusterIndex
                    bestSecondIndex = secondClusterIndex
        self.bestDistance = bestDistance
        return [bestFirstIndex, bestSecondIndex, bestDistance] 

    def betterDistance(self, firstClusterIndex, secondClusterIndex, bestDistance):
        if (self.clusterList[firstClusterIndex].isActive and
            self.clusterList[secondClusterIndex].isActive and
            self.distanceComparer.betterDistance(self.clusterDistanceMatrix[firstClusterIndex][secondClusterIndex], bestDistance)):
            return True
        return False

    def getActiveClusters(self):
        activeClusterList = []
        for cluster in self.clusterList:
            if (cluster.isActive):
                activeClusterList.append(cluster)
        return activeClusterList
    
    
