from operator import itemgetter
import time
import re
import subprocess
import math
import ptpTools
import collections
import configobj
import sys
import os
import itertools
import ptpDrpClusters
import ptpPymol
from Bio import PDB
import logging
log = logging.getLogger("cluster_annotation")

class ClusterAnnotator(ptpTools.Runner):

    def readUniversalAnnotation(self, clusterIndexList=None, drpSkipList=None):

        self.clusterSet = ptpDrpClusters.DrpClusterSet(clusterIndexList)
        print 'read member file'
        self.clusterSet.readClusterMemberFile(self.config["cluster_member_file"], drpSkipList)
        print 'cluster centroids'
        self.clusterSet.readClusterCentroids(self.config["cluster_stat_file"])

        self.clusterSet.readFastaFile(self.config["drp_fasta_file"])

        self.clusterSet.readSingletonFiles(self.config["longer_singleton_file"], self.config["shorter_singleton_file"], self.config["singleton_merge_cutoff"])
        print 'done universal'
class ClusterTextAnnotator(ClusterAnnotator):

    def readClusterFile(self):
        
        self.readUniversalAnnotation()

        self.clusterSet.readUniprotFiles(self.config["pdb_to_uniprot_file"], self.config["uniprot_fasta_file"])
        print 'read uniprot'
        
        self.clusterSet.readSpeciesFile(self.config["species_file"])

        self.clusterSet.readDisulfideBondFile(self.config["drp_list_file"])
        print 'done read disulfide'
        self.clusterSet.readDrpComplexFile(self.config["drp_complex_file"])
        print 'done read drp complex'
        self.clusterSet.readDistanceFile(self.config["distance_file"])
        print 'done read distance'
        self.clusterSet.readClusterFoldFile(self.config["cluster_fold_file"])
        print 'done read fold'
        self.clusterSet.readSequenceLengthFile(self.config["drp_list_file"])
        print 'done read length'
    def readIdenticalDrps(self):
        filterFile = self.config["identical_drp_cluster_file"]
        filteredClusterSet = ptpDrpClusters.DrpClusterSet()
        filteredClusterSet.readClusterMemberFile(filterFile)
        for cluster in self.clusterSet.getClusterList():
            for drp in cluster.getMemberList():
                identicalCluster = filteredClusterSet.getCluster(filteredClusterSet.getClusterIndexForDrpCode(drp.drpCode))
                identicalDrpList = []
                for identicalDrpCode in identicalCluster.getDrpCodeList():
                    if (identicalDrpCode != drp.drpCode):
                        identicalDrpList.append(identicalDrpCode)

                drp.setIdenticalDrpCodeList(identicalDrpList)

    def writeAnnotationFile(self):
        self.readIdenticalDrps()

        outputFile = self.getFullDatedOutputFile("cluster_annotation.txt")
        outputFh = open(outputFile, 'w')
        columnHeaderLine = ["DRP Code", "Cluster Index", "DRP Name", "Cluster Name", "Singleton?", "PDB Length", 
                            "Distance to Centroid", "Uniprot Accession", "Species", "Full Sequence Length", 
                            "PDB Sequence", "Disulfide Bond Count", "Identical DRPs"]
        outputString = '\t'.join(str(x) for x in columnHeaderLine)
        outputFh.write("%s\n" % outputString)
        for cluster in self.clusterSet.getClusterList():
            for drp in cluster.getMemberList():
                
                outputLine = [drp.drpCode, str(cluster.getClusterIndex() + 1), drp.getUniprotName(), cluster.getFoldName(), 
                              drp.getSingletonAnnotationString(), drp.getSequenceLength(),  drp.getDistanceToCentroid(), 
                              drp.getUniprotAccession(), drp.getSpecies(), drp.getUniprotLength(), drp.getSequence(), 
                              drp.getDisulfideBondCount(), drp.getIdenticalDrpCodeListString()]
                outputString = '\t'.join(str(x) for x in outputLine)
                outputFh.write("%s\n" % outputString)
                    

    def writeStatsForPaper(self):
        self.writeKnottinReclustering()

        self.writeDrpPerSpecies()

        self.writeDisulfideBondPatterns()

        self.writeLengthTable()

    def writeLengthTable(self): #added 1/20/16
        outputFile = self.getFullDatedOutputFile("cluster_lengths.txt")
        outputFh = open(outputFile, 'w')
        for cluster in self.clusterSet.getClusterList():
            lengthList = []
            for drp in cluster.getMemberList():
                lengthList.append(drp.getSequenceLength())
            avg = round(float(sum(lengthList)) / float(len(lengthList)), 1)
            outputList = [cluster.getClusterIndex() + 1, avg, len(lengthList)]
            outputFh.write('%s\n' % '\t'.join(str(x) for x in outputList))
        outputFh.close()

    def writeDrpPerSpecies(self):
        from collections import defaultdict
        outputFile = self.getFullDatedOutputFile("drpPerSpecies.txt")
        outputFh = open(outputFile, 'w')
        columnList = ["Cluster Index"]
        for level in ["Kingdom", "Phylum",  "Class", "Order", "Family", "Genus", "Species"]:
            columnList.append("Unique %s" % level)
            columnList.append("DRPs per %s" % level)
        outputFh.write("%s\n" % '\t'.join(columnList))

        for cluster in self.clusterSet.getClusterList():
            levelToCount = defaultdict(dict)
        
            for drp in cluster.getMemberList():
                index = cluster.getClusterIndex()
                for level in ["kingdom", "phylum",  "tClass", "order", "family", "genus", "species"]:
                    tt = drp.getTaxonTuple()
                    if (tt):
                        value = getattr(tt, level)
                        if (value.strip() == ''):
                            continue
                    else:
                        continue
                    if (value not in levelToCount[level]):
                        levelToCount[level][value] = 0
                    levelToCount[level][value] += 1
            #levelToCount[species][homo_sapien] = 7
            outputList = [cluster.getClusterIndex() + 1]
            for level in ["kingdom", "phylum",  "tClass", "order", "family", "genus", "species"]:
                totalCount = sum(levelToCount[level].values())
                levelCount = len(levelToCount[level])
                print 'cluster %s level: %s' % (cluster.getClusterIndex() + 1, level)
                print '\n'.join(levelToCount[level].keys())
                if (levelCount == 0):
                    outputList += [levelCount, 'N/A']
                else:
                    outputList += [levelCount, round(float(totalCount) / float(levelCount), 2)]                    

            outputFh.write("%s\n" %  '\t'.join(str(x) for x in outputList))                


    def writeDisulfideBondPatterns(self):
        outputFile = self.getFullDatedOutputFile("disulfideBondPatterns.txt")
        outputFh = open(outputFile, 'w')
        clusterIndexToMotif = {}
        for clusterIndex in self.config['disulfide_bond_pattern_clusters']:
            cluster = self.clusterSet.getCluster(int(clusterIndex) - 1)
            for drp in cluster.getMemberList():
                clusterIndex = cluster.getClusterIndex()
                drpCode = drp.drpCode
                sequence = drp.getSequence()

                if (clusterIndex not in clusterIndexToMotif):
                    clusterIndexToMotif[clusterIndex] = {}
                if (sequence.count('C') == 3):
                    continue
                    #print "skipping line %s" % line
                elif (sequence.count('C') == 2):
                    motifRegex = re.compile('.*C(.*)C.*')
                    foundMatch = motifRegex.search(sequence)
                    if (foundMatch):
                        innerString = foundMatch.group(1)
                        motif = len(innerString)
                        #print "got motif %s string %s from seq %s" % (motif, innerString, sequence)
                        if (motif not in clusterIndexToMotif[clusterIndex]):
                            clusterIndexToMotif[clusterIndex][motif] = []
                        clusterIndexToMotif[clusterIndex][motif].append(innerString)
        for clusterIndex in sorted(clusterIndexToMotif):
            for motif in sorted(clusterIndexToMotif[clusterIndex]):
                outputList = [clusterIndex + 1, cluster.getFoldName(), "CX%sC" % motif, len(clusterIndexToMotif[clusterIndex][motif])]
                outputFh.write("%s\n" % '\t'.join(str(x) for x in outputList))


    def writeKnottinReclustering(self):

        nativeOverlapClusterFile = self.config["native_overlap_cluster_member_file"]
        knottinClusterFile = self.config["knottin_recluster_cluster_member_file"]

        nativeOverlapSet = ptpDrpClusters.DrpClusterSet()
        nativeOverlapSet.readClusterMemberFile(nativeOverlapClusterFile)

        knottinSet = ptpDrpClusters.DrpClusterSet()
        knottinSet.readClusterMemberFile(knottinClusterFile)

        initialKnottinClusterIndexList = self.config["knottin_cluster_index_list"]
        knottinClusterIndexList = [int(x) - 1 for x in initialKnottinClusterIndexList]
        outputFile = self.getFullDatedOutputFile("knottin_reclustering.txt")
        outputFh = open(outputFile, 'w')
        knottinStats = {}
        for clusterIndex in knottinSet.getClusterIndexList():
            knottinStats[clusterIndex] = {}
        columnList = [""]
        for oldKci in knottinClusterIndexList:
            oldCluster = nativeOverlapSet.getCluster(oldKci)
            for clusterIndex in knottinSet.getClusterIndexList():
                knottinStats[clusterIndex][oldKci] = 0
            for drp in oldCluster.getDrpCodeList():
                foundDrp = 0
                for knottinCluster in knottinSet.getClusterList():
                    if (knottinCluster.hasDrp(drp)):
                        foundDrp = 1
                        knottinStats[knottinCluster.getClusterIndex()][oldCluster.getClusterIndex()] += 1
            columnList.append(oldKci + 1)
        outputFh.write("%s\n" % '\t'.join(str(x) for x in columnList))

        for newCluster in knottinStats:
            outputList = [newCluster + 1]
            for (oldCluster, count) in sorted(knottinStats[newCluster].iteritems(), key=itemgetter(0)):
                outputList.append("%s | %s" % (knottinStats[newCluster][oldCluster], nativeOverlapSet.getCluster(oldCluster).getSize()))
            outputFh.write("%s\n" % '\t'.join(str(x) for x in outputList))


    
