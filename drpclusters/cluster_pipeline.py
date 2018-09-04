from operator import itemgetter
import time
import re
import subprocess
import logging
import math
import collections
import argparse
import sys
import os
import itertools
import cluster_lib
from operator import itemgetter


class ClusterPipelineRunner(cluster_lib.Runner):
    def __init__(self, config):
        super(ClusterPipelineRunner, self).__init__(config)
        self.nativeOverlapDm = self.readDistanceFile(self.config.native_overlap_distance_file)
        self.shorterFractionDm = self.readDistanceFile(self.config.shorter_fraction_distance_file)

    def execute(self):

        print "Begin step filter sequence"
        self.filterSequences()

        print "\nBegin step cluster native overlap"
        self.clusterNativeOverlap()

        print "\nBegin step make knottin cluster subset"
        self.makeKnottinClusterSubset()

        print "\nBegin step merge knottin clusters"
        self.mergeKnottinClusters()

        print "\nBegin step longer singletons"
        self.processLongerFractionSingletons()

        print "\nBegin step shorter singletons"
        self.processShorterFractionSingletons()

    """
    Class that implements the primary steps of the clustering pipeline.
    In order these are:
    1. Filter DRPs by 100% sequence and structure identity
    2. Cluster the filtered DRPs hierarchically by their native overlap similarity
    3. Identify clusters containing DRPs annotated primarily as knottins 
    4. Recluster knottin clusters by the amount their disulfide bonds overlap to consolidate them
    5. Merge the knottin clusters with the non-knottin clusters
    6. Consolidate clusters further by processing 'longer fraction' singletons
    7. Consolidate even furthre by processing 'shorter fraction' singletons
    
    These steps are called by an executable (i.e. ../bin/runFullClusterPipeline.py). At any step, the results of
    the previous clustering step can be read in from a file (i.e., 'readPreviousFilterClusters'). All clusters 
    are handled as DrpClusterSets and can be compared to previous runs using self.compareToPreviousClusters()
    """

    #DELETE
    def readShorterFractionDistanceFile(self):
        if (self.shorterFractionDm is None):
            distanceFile = self.config.shorter_fraction_distance_file
            self.shorterFractionDm = self.readDistanceFile(distanceFile)
        return self.shorterFractionDm

    #DELETE
    def readNativeOverlapDistanceFile(self):
        if (not self.nativeOverlapDm):
            distanceFile = self.config.native_overlap_distance_file
            self.nativeOverlapDm = self.readDistanceFile(distanceFile)
        return self.nativeOverlapDm

    
    def readDisulfideDistances(self):
        disulfideDistanceMatrix = cluster_lib.DisulfideDrpDistanceMatrix()
        disulfideDistanceMatrix.readFile(self.config.disulfide_distance_file)
        return disulfideDistanceMatrix

    def readDistanceFile(self, distanceFile):
        distanceMatrix = cluster_lib.DrpDistanceMatrix()
        distanceMatrix.readFile(distanceFile)
        return distanceMatrix

    #DELETE
    def createUnfilteredDrpList(self):
        """Read the results of the PDB query that obtained all DRPs. The current parameters are all solved structures
        between 1 and 50 residues and 1-4 disulfide bonds
        """
        dqi = cluster_lib.DrpQueryInfo()
        dqi.readFile(self.config.drp_query_file)
        drpList = dqi.createDrpList()
        return drpList

    def filterSequences(self):
        """Filter DRPs by 100% sequence and structure identity"""
        filteredSeqsDistances = self.readDistanceFile(self.config.filter_seqs_distance_file)

        reader = cluster_lib.PtFileReader(self.config.drp_query_file)
        drpList = [cluster_lib.Drp(x) for x in reader.getLines()]
        
        ldc = cluster_lib.LargerDistanceComparer(self.config.filtering_cutoff)
        
        self.filteredClusterSet = self.clusterDrps(drpList, filteredSeqsDistances, ldc)
        
        self.filteredDrps = self.getFilteredClusterRepresentatives(self.filteredClusterSet)

        self.filteredClusterSet.writeClusterMemberFile(self.getFullOutputFile("filterSequences_cluster_members.txt"))

    def clusterDrps(self, drpList, distanceMatrix, comparer, clusterOrder=None):
        """Run the main hierarchical clustering step.
        
        @param clusterOrder -- order in which to initialize DRPs, returned by self.readClusterOrder()
        """
        clusterList = []
        print 'type of drps is %s' % drpList[0]
        sortedDrpList = self.createDrpListForClustering(drpList, distanceMatrix, clusterOrder)

        for (i, drp) in enumerate(sortedDrpList):
            cluster = cluster_lib.DrpCluster(i)
            cluster.addMember(drp)
            clusterList.append(cluster)

        hc = cluster_lib.HierarchicalClusters(distanceMatrix, clusterList, comparer)
        
        hc.clusterHierarchical()

        clusterSet = cluster_lib.DrpClusterSet()
        for cluster in hc.getActiveClusters():
            clusterSet.addCluster(cluster)
            
        return clusterSet

    def createDrpListForClustering(self, drpList, distanceMatrix, clusterOrder):
        """Return a sorted list of DRPs. Either sort alphabetically or according to <clusterOrder>"""
        initialDrpList = []
        for drp in drpList:
            if (distanceMatrix.containsId(drp.drpCode)):
                initialDrpList.append(drp)
        sortedDrpList = sorted(initialDrpList, key=lambda x: x.drpCode)
        
        if (clusterOrder is not None):
            sortedDrpList = sorted(initialDrpList, key=lambda x: clusterOrder.index(x.drpCode))
        return sortedDrpList

    def getFilteredClusterRepresentatives(self, filteredClusterSet):
        """Return a list of filtered DRPs from <filteredClusterSet>"""
        repDrps = []
        for cluster in filteredClusterSet.getClusterList():
            repDrps.append(cluster.getSortedFirstMember())
        return repDrps

    def clusterNativeOverlap(self):
        """Cluster the set of filtered DRPs according to their native overlap (i.e., the fraction of equivalent
        alpha-carbon atoms within 3.5 angstroms of each other across a pair of DRPs
        """
        ldc = cluster_lib.LargerDistanceComparer(self.config.native_overlap_cutoff)

        self.nativeOverlapClusterSet = self.clusterDrps(self.filteredDrps, self.nativeOverlapDm, ldc)

        self.nativeOverlapClusterSet.renumberBySize()

        self.nativeOverlapClusterSet.writeClusterMemberFile(self.getFullOutputFile("nativeOverlap_cluster_members.txt"))
        
    def identifyKnottinClusters(self):
        """Return the list of cluster indexes for clusters primarily made up of Knottins as annotated by SCOP"""
        knottinClusterIndexList = []
        
        for noc in self.nativeOverlapClusterSet.getClusterList():
            knottinCount = 0
            for drp in noc.getMemberList():
                if (drp.getScopFold() == "g.3"):
                    knottinCount += 1
            if (self.passesKnottinCutoffs(knottinCount, len(noc.getMemberList()))):
                knottinClusterIndexList.append(noc.getClusterIndex())

        return knottinClusterIndexList

    def passesKnottinCutoffs(self, knottinCount, totalCount):
        """Determine if there are enough knottins in this cluster for it to be considered a knottin cluster"""
        knottinFraction = float(knottinCount) / float(totalCount)
        cutoffFraction = float(self.config.knottin_cutoff_fraction)
        cutoffAbsolute = int(self.config.knottin_cutoff_absolute)
        
        if (knottinFraction >= cutoffFraction and knottinCount >= cutoffAbsolute):
            return True
        return False

    def makeKnottinDrpList(self, knottinClusterIndexList):
        """Return a list of DRP codes for all Knottin DRPs in <knottinClusterIndexList>"""
        knottinList = []
        for cluster in self.nativeOverlapClusterSet.getClusterList():
            if (cluster.getClusterIndex() in knottinClusterIndexList):
                knottinList += cluster.getMemberList()

        return knottinList
        
    def populateScopData(self, clusterSet):
        """Read SCOP file and assign SCOP annotation to all DRPs in the native overlap cluster set"""
        self.scopInfo = cluster_lib.ScopInfo(self.config.scop_file)
        for cluster in clusterSet.getClusterList():
            for drp in cluster.getMemberList():
                drp.setScopFold(self.scopInfo.getScopFoldId(drp.drpCode))
                
    def makeKnottinClusterSubset(self):
        """Recluster those clusters primarily containing knottin DRPs according to the distances between their 
        equivalent disulfide bonds
        """
        self.populateScopData(self.nativeOverlapClusterSet)
        
        self.knottinClusterIndexList = self.identifyKnottinClusters()
        
        knottinList = self.makeKnottinDrpList(self.knottinClusterIndexList)
        disulfideDistanceMatrix = cluster_lib.DisulfideDrpDistanceMatrix()
        disulfideDistanceMatrix.readFile(self.config.disulfide_distance_file)
        disulfideDistanceMatrix = self.readDisulfideDistances()

        sdc = cluster_lib.SmallerDistanceComparer(self.config.knottin_disulfide_cutoff)

        self.knottinClusterSet = self.clusterDrps(knottinList, disulfideDistanceMatrix, sdc)
        self.knottinClusterSet.renumberBySize()
        self.knottinClusterSet.writeClusterMemberFile(self.getFullOutputFile("clusterKnottins_cluster_members.txt"))

    def mergeKnottinClusters(self):
        """Merge clusters in self.nativeOverlapClusterSet not containing knottins with clusters in self.knottinClusterSet
        and renumber them by size
        """
        self.mergedClusterSet = cluster_lib.DrpClusterSet()
        nonKnottinClusters = [x for x in self.nativeOverlapClusterSet.getClusterList()
                              if x.getClusterIndex() not in self.knottinClusterIndexList]
        
        for (i, cluster) in enumerate(nonKnottinClusters + self.knottinClusterSet.getClusterList()):
            cluster.setClusterIndex(i)
            self.mergedClusterSet.addCluster(cluster)
            
        self.mergedClusterSet.renumberBySize()

        self.mergedClusterSet.writeClusterMemberFile(self.getFullOutputFile("mergeKnottins_cluster_members.txt"))

    def processTopClusterOrder(self, topClusterSet, paramFileName):
        """Return the order in which DRPs will be considered as the reference DRPs for singletons"""
        drpList = []
        for cluster in topClusterSet.getClusterList():
            drpList += cluster.getMemberList()
        return drpList
    
    def processLongerFractionSingletons(self):
        """Take all DRPs in the less populated clusters and add them to the greater populated cluster that has the closest 
        reference DRP passing a cutoff, taking the length of the longer DRP as the reference in the similarity fraction
        """
        [topClusterSet, bottomClusterSet] = self.divideClusters(self.mergedClusterSet)

        topClusterDrpList = list(itertools.chain(*[x.getMemberList() for x in topClusterSet.getClusterList()]))

        singletonProcessor = cluster_lib.LongerFractionSingletonProcessor(self.nativeOverlapDm, topClusterDrpList, bottomClusterSet, self.config.singleton_cutoff)
        singletonProcessor.processSingletons()

        self.longerSingletonClusterSet = self.addSingletonsToClusters(self.mergedClusterSet)

        self.longerSingletonClusterSet.writeClusterMemberFile(self.getFullOutputFile("processLongerSingletons_cluster_members.txt"))
        
        self.longerSingletonClusterSet.writeSingletonFile(self.getFullOutputFile("processLongerSingletons_singleton_pairs.txt"), 'longer')

    def addSingletonsToClusters(self, inputClusterSet):
        """Remove the DRPs newly identified as singletons from their current clusters and add each to the cluster
        containing its reference DRP.
        """
        [topClusterSet, bottomClusterSet] = self.divideClusters(inputClusterSet)
        newClusterSet = cluster_lib.DrpClusterSet()

        #process bottom clusters
        for cluster in bottomClusterSet.getClusterList():
            newCluster = cluster_lib.DrpCluster(cluster.getClusterIndex())
            for drp in cluster.getMemberList():
                if (drp.isSingleton()):
                    #drp gets added to top cluster
                    refDrpCode = drp.getSingletonTuple().refDrp
                    topClusterIndex = topClusterSet.getClusterIndexForDrpCode(refDrpCode)
                    topClusterSet.getCluster(topClusterIndex).addMember(drp)
                else:
                    #drp still in bottom cluster
                    newCluster.addMember(drp)
            if (len(newCluster.getMemberList()) > 0):
                newClusterSet.addCluster(newCluster)

        #process top clusters
        for cluster in topClusterSet.getClusterList():
            newClusterSet.addCluster(cluster)

        newClusterSet.renumberBySize()
        return newClusterSet

    def readDrpLengthFile(self, fn):
        drpToLength = {}
        pfr = cluster_lib.PtFileReader(fn)
        for line in pfr.getLines():
            [drpName, length] = line.split('\t')
            drpToLength[drpName] = int(length)
            
        return drpToLength
                        
    def processShorterFractionSingletons(self):
        """Take all DRPs in the less populated clusters and add them to the greater populated cluster that has the closest 
        reference DRP passing a cutoff, taking the length of the shorter DRP as the reference in the similarity fraction
        """
        [topClusterSet, bottomClusterSet] = self.divideClusters(self.longerSingletonClusterSet)
        drpToLengthMap = self.readDrpLengthFile(self.config.drp_length_file)
        topClusterDrpList = list(itertools.chain(*[x.getMemberList() for x in topClusterSet.getClusterList()]))
        singletonProcessor = cluster_lib.ShorterFractionSingletonProcessor(self.shorterFractionDm, topClusterDrpList, bottomClusterSet, self.config.singleton_cutoff, drpToLengthMap)
        singletonProcessor.processSingletons()

        self.shorterSingletonClusterSet = self.addSingletonsToClusters(self.longerSingletonClusterSet)

        self.shorterSingletonClusterSet.writeClusterMemberFile(self.getFullOutputFile("processShorterSingletons_cluster_members.txt"))

        for nextCluster in sorted(self.shorterSingletonClusterSet.getClusterList(), key=lambda x: len(x.getMemberList()), reverse=True):
            print "cluster %s length %s" % (nextCluster.getClusterIndex(), len(nextCluster.getMemberList()))

        self.shorterSingletonClusterSet.writeSingletonFile(self.getFullOutputFile("processShorterSingletons_singleton_pairs.txt"), 'shorter')

    def divideClusters(self, inputClusterSet):
        """Divide the clusters into a 'top' and 'bottom' set specified by the top_cluster_count parameter.
        If there is a tie for the last top cluster, add all clusters involved in the tie. 
        """
        topClusterIndex = int(self.config.top_cluster_count) - 1
        topClusterSet = cluster_lib.DrpClusterSet()
        bottomClusterSet = cluster_lib.DrpClusterSet()
        sizeOfTopCluster = len(inputClusterSet.getCluster(topClusterIndex).getMemberList())
        for nextCluster in sorted(inputClusterSet.getClusterList(), key=lambda x: len(x.getMemberList()), reverse=True):
            if (sizeOfTopCluster <= len(nextCluster.getMemberList())):
                topClusterSet.addCluster(nextCluster)
            else:
                bottomClusterSet.addCluster(nextCluster)
        return [topClusterSet, bottomClusterSet]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-r", "--run_directory", help="Directory to which output is written. Will be created if does not exists", required=True)
    parser.add_argument("-q", "--drp_query_file", help="Text file with input set of DRPs. One DRP per line, specified as a DRP code (5 characters; first 4 are PDB ID and 5th is chain, eg 1zdcA)", required=True)
    parser.add_argument("-f", "--filter_seqs_distance_file", help="Distance matrix for similarity product filtering (see documentation accompanying publication for details on creating this file)", required=True)
    parser.add_argument("-n", "--native_overlap_distance_file", help="Distance matrix for native overlap clustering (see documentation accompanying publication for details on creating this file)", required=True)
    parser.add_argument("-d", "--disulfide_distance_file", help="Distance matrix for knottin reclustering (see documentation accompanying publication for details on creating this file)", required=True)
    parser.add_argument("-s", "--shorter_fraction_distance_file", help="Distance matrix for shorter singleton postprocessing step (see documentation accompanying publication for details on creating this file)", required=True)
    parser.add_argument("-c", "--filtering_cutoff", help="Cutoff for the initial filtering step (recommended to set to 99 to filter identical DRPs", required=True, default=99)
    parser.add_argument("-k", "--knottin_disulfide_cutoff", help="Cutoff for Knottin reclustering step in Angstroms", required=True, default=2.0)
    parser.add_argument("-t", "--native_overlap_cutoff", help="Cutoff for the initial native overlap clustering step", required=True, default=0.70)
    parser.add_argument("-v", "--knottin_cutoff_fraction", help="Fraction of DRPs in a cluster annotated as Knottins by SCOP for it to be considered a 'knottin' cluster and thus be reclustered in the disulfide overlap step (cluster must also pass knottin_cutoff_absolute filter)", required=True, default=.01)
    parser.add_argument("-b", "--knottin_cutoff_absolute",  help="Number of DRPs in a cluster annotated as Knottins by SCOP for it to be considered a 'knottin' cluster and thus be reclustered in the disulfide overlap step (cluster must also pass knottin_cutoff_fraction filter)", required=True, default=4)
    parser.add_argument("-l", "--top_cluster_count", required=True, default=25, help="In the singleton steps for the pipeline, DRPs in the top N clusters sorted by size are set as references for DRPs not in those clusters. Use this to set N")
    parser.add_argument("-g", "--singleton_cutoff", required=True, default=0.7, help="In the singleton steps for the pipeline, DRPs are added from the smaller clusters to larger clusters if they are within this cutoff. Recommend to be identical to native_overlap_cutoff")
    parser.add_argument("-e", "--drp_length_file", required=True, help="tab-delimited file of DRP codes to sequence lengths")
    parser.add_argument("-p", "--scop_file", required=True, help="SCOP fold file, format identical to http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.des.scop.txt_1.75")

    if (len(sys.argv) < 2):
        print "Please run with '-h' for full usage"
        sys.exit()
    config = parser.parse_args(sys.argv[1:])

    runner = ClusterPipelineRunner(config)
    runner.execute()