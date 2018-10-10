from operator import itemgetter
import time
import re
import subprocess
import math
import shutil
import collections
import argparse
import sys
import os
import itertools
import cluster_lib
import logging

class ClusterVisAnnotationRunner(cluster_lib.Runner):
    def initSubclass(self):
      
        self.shutil = cluster_lib.PtpShutils()

        self.drpOrder = {}

    def sleepUntilDone(self, fileName, predicate):
        """Sleep until predicate involving fileName is true; useful for avoiding race conditions in file manipulation"""
        sleepTime = 0
        while (predicate(fileName)):
            print "sleep 1 second"
            time.sleep(1)
            sleepTime += 1
            if (sleepTime > 10):
                raise Exception("Timeout on file %s" % fileName)
    
    def readClusterFile(self):
        self.clusterSet = cluster_lib.DrpClusterSet(self.config.cluster_index_list)

        self.clusterSet.readClusterMemberFile(self.config.cluster_member_file)

        self.clusterSet.readSingletonFiles(self.config.longer_singleton_file, self.config.shorter_singleton_file, self.config.singleton_merge_cutoff)

    def runFinalAlignments(self):
        self.withinClusterAl2CoResults = {}
        for cluster in self.clusterSet.getClusterList():
            clusterOutputDir = self.getFullOutputFile(str(cluster.getClusterIndex() + 1))
            if (not os.path.exists(clusterOutputDir)):
                os.mkdir(clusterOutputDir)
            
            cluster.salignDrps(self.config.drp_pdb_directory, None, True)
            
            for drpCode in cluster.getDrpCodeList():
                sourcePdb = os.path.join(os.getcwd(), "%s_fit.pdb" % drpCode)
                if (os.path.exists(sourcePdb)):
                    finalFile = os.path.join(clusterOutputDir, "%s_fit.pdb" % drpCode)
                    shutil.move(sourcePdb, finalFile)

    def makeLogPrefix(self):
        return "clusterVisAnnotation"

def getParser():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="Align all DRPs in each cluster for visualization purposes")
    
    parser.add_argument("-r", dest="run_directory", metavar="<directory>", required=True,
                        help="Directory to which output is written. Will be created if does not exists")
    
    parser.add_argument("-i", dest="cluster_index_list", metavar="<space-separated list>", nargs='+', required=True,
                        help="List of cluster indices to align")
    
    parser.add_argument("-c", dest="cluster_member_file", metavar="<file>", required=True,
                        help="Final 'cluster_members.txt' file written by cluster_pipeline.py (usually shorter singletons file")
    
    parser.add_argument("-l", dest="longer_singleton_file", metavar="<file>", required=True,
                        help="processLongerSingletons_singleton_pairs.txt file written by cluster_pipeline.py")
    
    parser.add_argument("-f", dest="shorter_singleton_file", metavar="<file>", required=True,
                        help="processLongerSingletons_singleton_pairs.txt file written by cluster_pipeline.py")
    
    parser.add_argument("-m", dest="singleton_merge_cutoff", metavar="<float>", required=True,
                        help="Singleton merge cutoff; should be identical to the one set in cluster_pipeline.py")
    
    parser.add_argument("-p", dest="drp_pdb_directory", metavar="<directory>", required=True,
                        help="Location of PDB files. Created in setup_pdb.py; each file should be named after its DRP code\n\n")

    return parser

if __name__ == '__main__':

    parser = getParser()
    config = parser.parse_args(sys.argv[1:])

    runner = ClusterVisAnnotationRunner(config)

    runner.readClusterFile()
    
    runner.runFinalAlignments()
