import sys
import os
import unittest
import tempfile
import subprocess
import cluster_lib
import cluster_pipeline
import shutil
from nose.plugins.attrib import attr

class DrpStub(cluster_lib.Drp):
    def __init__(self, memberId):
        self.memberId = memberId
        
    def getClusterMemberId(self):
        return self.memberId

@attr(type='unit')
class TestHierarchicalClustering(unittest.TestCase):
    def test_clustering(self):
        dm = cluster_lib.DistanceMatrix()
        for (id1, id2, d) in [('m1', 'm2', .20), ('m1', 'm3', .21), ('m1', 'm4', .22), ('m1', 'm5', .23),
                              ('m2', 'm3', .89), ('m2', 'm4', .66), ('m2', 'm5', .67),
                              ('m3', 'm4', .68), ('m3', 'm5', .69),
                              ('m4', 'm5', .9)]:
            
            dm.addDistancePair(id1, id2, d)
            dm.addDistancePair(id2, id1, d)
        cl = []
        for (i, m) in enumerate(['m1', 'm2', 'm3', 'm4', 'm5']):
            c = cluster_lib.DrpCluster(i)
            c.addMember(DrpStub(m))
            cl.append(c)
        ldc7 = cluster_lib.LargerDistanceComparer(.7)
        hc = cluster_lib.HierarchicalClusters(dm, cl, ldc7)
        hc.clusterHierarchical()

        #cutoff of .7 clusters m2 and m3, m4 and m5, and m1 as a singleton
        expectedClusterMembers = {0: ['m1'], 1: ['m2', 'm3'], 3: ['m4', 'm5']}
        for c in hc.getActiveClusters():
            self.assertEquals(expectedClusterMembers[c.getClusterIndex()],
                              [x.getClusterMemberId() for x in c.getMemberList()])
        self.assertEquals(hc.clusterDistanceMatrix[1][3], .675)

        ldc5 = cluster_lib.LargerDistanceComparer(.5)
        hc = cluster_lib.HierarchicalClusters(dm, cl, ldc5)
        hc.clusterHierarchical()

        #cutoff of .5 is more lenient and merges clusters 1 and 3
        expectedClusterMembers = {0: ['m1'], 1: ['m2', 'm3', 'm4', 'm5']}
        for c in hc.getActiveClusters():
            self.assertEquals(expectedClusterMembers[c.getClusterIndex()],
                              [x.getClusterMemberId() for x in c.getMemberList()])

def runFullClusterPipeline(exampleDataDir, experimentDir):
    inputFileList = [os.path.join(exampleDataDir, 'exampleDistances', x) for x in ['drp_lengths.txt', 'similarity_product.txt', 'longer_fraction.txt', 'disulfides.txt', 'shorter_fraction.txt']]
    inputFileList += [os.path.join(exampleDataDir, x) for x in ['drp_list.txt', 'scop_family_assignment.txt']]

    subprocess.check_output(['cp'] + inputFileList + [experimentDir])
    p = cluster_pipeline.getParser()

    argList = ['-r', experimentDir,
               '-q', os.path.join(experimentDir, 'drp_list.txt'),
               '-f', os.path.join(experimentDir, 'similarity_product.txt'),
               '-n', os.path.join(experimentDir, 'longer_fraction.txt'),
               '-d', os.path.join(experimentDir, 'disulfides.txt'),
               '-s', os.path.join(experimentDir, 'shorter_fraction.txt'),
               '-e', os.path.join(experimentDir, 'drp_lengths.txt'),
               '-p', os.path.join(experimentDir, 'scop_family_assignment.txt'),
               '-c', '99', '-t', '.9', '-k', '2.0', '-v', '0.01',  '-b', '4', '-l', '7', '-g', '.7']
    config = p.parse_args(argList)
    runner = cluster_pipeline.ClusterPipelineRunner(config)
    runner.execute()        

@attr(type='unit')            
class ClusterPipelineTest(unittest.TestCase):
    def setUp(self):
        self.experimentDir = tempfile.mkdtemp()
        self.testModuleDir = os.path.dirname(os.path.abspath(sys.modules[self.__class__.__module__].__file__))
        self.exampleDataDir = os.path.join(self.testModuleDir, '..', 'exampleData')
        print self.experimentDir
    
    def test_cluster_pipeline(self):
        runFullClusterPipeline(self.exampleDataDir, self.experimentDir)
        pfr = cluster_lib.PtFileReader(os.path.join(self.experimentDir, 'processShorterSingletons_cluster_members.txt'))
        self.assertEquals(len(pfr.getLines()), 95)
        self.assertEquals(pfr.getLines()[-1], "8	2ny9X")

    def tearDown(self):
        shutil.rmtree(self.experimentDir)
        
@attr(type='integration')            
class ClusterIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.experimentDir = tempfile.mkdtemp()
        self.testModuleDir = os.path.dirname(os.path.abspath(sys.modules[self.__class__.__module__].__file__))
        self.miniCondaCmd = os.getenv('CONDA', 'python')
        #ideally the test dir would have its own data but since we already have an example, this avoids duplicating files
        self.exampleDataDir = os.path.join(self.testModuleDir, '..', 'exampleData')
        self.scriptDir = os.path.join(self.testModuleDir, '..')
        print self.experimentDir

    def tearDown(self):
        shutil.rmtree(self.experimentDir)
        
    def copyPdbInput(self):
        print 'copy pdb'
        subprocess.check_output(['cp', os.path.join(self.exampleDataDir, 'dividedPdb.tar.gz'), self.experimentDir])
        print 'unzip pdb'
        subprocess.check_output(['tar', '-C', self.experimentDir, '-xzf', os.path.join(self.experimentDir, 'dividedPdb.tar.gz')])

    def copyDrpListFile(self):
        subprocess.check_output(['cp', os.path.join(self.exampleDataDir, 'drp_list.txt'), self.experimentDir])

    def runFullSetupPdb(self):
        self.copyPdbInput()
        self.copyDrpListFile()
        os.mkdir(os.path.join(self.experimentDir, 'drpDir'))
        subprocess.check_output([self.miniCondaCmd, os.path.join(self.scriptDir, 'setup_pdb.py'),
                                 '-r', self.experimentDir,
                                 '-q', os.path.join(self.experimentDir, 'drp_list.txt'),
                                 '-p', os.path.join(self.experimentDir, 'dividedPdb'),
                                 '-l', os.path.join(self.experimentDir, 'drp_lengths.txt'),
                                 '-o', os.path.join(self.experimentDir, 'drpDir')])

    def test_setup_pdb(self):
        print self.experimentDir
        self.runFullSetupPdb()
        self.assertTrue(os.path.exists(os.path.join(self.experimentDir, 'drpDir', '1a0mB.pdb')))
        
    def test_align_drps(self):
        self.runFullSetupPdb()
        subprocess.check_output([self.miniCondaCmd, os.path.join(self.scriptDir, 'pairwise_align.py'),
                                 '-r', self.experimentDir,                                 
                                 '-q', os.path.join(self.experimentDir, 'drp_list.txt'),
                                 '-p', os.path.join(self.experimentDir, 'drpDir'),
                                 '-m', 'full_drp',
                                 '-x', '3',
                                 '-o', os.path.join(self.experimentDir, 'pairwise.txt')])
        pfr = cluster_lib.PtFileReader(os.path.join(self.experimentDir, 'pairwise.txt'))
        self.assertEquals(len(pfr.getLines()), 39)
        self.assertEquals(pfr.getLines()[0], '1m2sA	2ktxA	native_overlap	36')

        subprocess.check_output([self.miniCondaCmd, os.path.join(self.scriptDir, 'pairwise_align.py'),
                                 '-r', self.experimentDir,                                 
                                 '-q', os.path.join(self.experimentDir, 'drp_list.txt'),
                                 '-p', os.path.join(self.experimentDir, 'drpDir'),
                                 '-m', 'disulfides',
                                 '-x', '3',
                                 '-o', os.path.join(self.experimentDir, 'disulfides.txt')])
        pfr = cluster_lib.PtFileReader(os.path.join(self.experimentDir, 'disulfides.txt'))
        self.assertEquals(len(pfr.getLines()), 3)
        self.assertEquals(pfr.getLines()[0], "1m2sA	2ktxA	cysteine_rms	0.994954764843	['8', '29', '14', '34', '18', '36']	['8', '28', '14', '33', '18', '35']")

    def test_vis_annotation(self):
        self.runFullSetupPdb()
        runFullClusterPipeline(self.exampleDataDir, self.experimentDir)

        subprocess.check_output([self.miniCondaCmd, os.path.join(self.scriptDir, 'cluster_vis_annotation.py'),
                                 '-r', os.path.join(self.experimentDir, 'visAnnotation'),
                                 '-c', os.path.join(self.experimentDir, 'processShorterSingletons_cluster_members.txt'),
                                 '-l', os.path.join(self.experimentDir, 'processLongerSingletons_singleton_pairs.txt'),
                                 '-f', os.path.join(self.experimentDir, 'processShorterSingletons_singleton_pairs.txt'),
                                 '-p', os.path.join(self.experimentDir, 'processShorterSingletons_singleton_pairs.txt'),                                 
                                 '-p', os.path.join(self.experimentDir, 'drpDir'),
                                 '-i', '1', '2', '3', '4', '5', '6', '7', '8', '-m', '.7'])
        fileList = os.listdir(os.path.join(self.experimentDir, 'visAnnotation', '1'))
        self.assertEquals(len([x for x in fileList if x.endswith('_fit.pdb')]), 29)

        
if __name__ == '__main__':
    unittest.main()
        
