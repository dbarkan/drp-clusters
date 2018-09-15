# drp-clusters
This project presents an automated protocol for clustering small disulfide-rich peptides (DRPs) using amino acid sequence and molecular structure features. Applying this protocol on DRPs in the Protein Databank yields insight into their evolution, and also results in a small, diverse set of representative DRPs that provide basis for further bioengineering efforts such as phage display. A detailed description of this approach, as well as its results and experimental application is available in [DT Barkan *et al*, *BMC Bioinformatics*, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27881076). This code is available as supplemental files in that project; it is reproduced here for increased accessibility.

todo: image of clustered DRPs??

## Prerequisites
- drp-clusters is run with Python 2.7.15 in a 64-bit Linux environment
- [MODELLER libraries](https://salilab.org/modeller/) are required for large components of the software. 

## Installing
- The drpclusters package can be installed using pip:

```pip install drpclusters```

- MODELLER is available for download [here](https://salilab.org/modeller/download_installation.html). If you use the [Ana/Miniconda scientific Python distribution](https://www.continuum.io/downloads), installing MODELLER as follows:

```conda config --add channels salilab
conda install modeller
```

Other ways to install MODELLER are listed in the above link.

MODELLER is free for academic use provided you [register for a license](https://salilab.org/modeller/registration.html). You will be prompted after installation to edit a file to add your Modeller license key.

If you run MODELLER using the python executable provided with your Conda installation, you may have to tell it where to find drpclusters, as it may ignore your PYTHONPATH. Do this either by specifying a target when installing with pip:

```pip install --target=<conda_root>/lib/python2.7/site-packages/```

todo -- make sure that works

If you install drpclusters elsewhere, it can be imported by your Conda installation as follows:

1. Create the following file: `<conda_root>/lib/python2.7/site-packages/drpclusters.pth`
2. Add one line to this file listing the full path to the drpclusters package: <full_path>/drpclusters/


## Running the drpclusters pipeline
### Overview
Below is a step-by-step guide to running the full pipeline start to finish. The drp-clusters package includes an example set of 100 DRP PDB files for demonstration and testing purposes. Each step below includes the commands to apply the pipeline to this example set of DRPs. 

Most of these steps require using libraries from MODELLER. An exception is step xyz, which runs the actual clustering pipeline. If you don't have access to MODELLER but want to demo the pipeline, example distance matrix files are also provided. 

In this example, all steps are run in the `drp-clusters/example/` directory, which comes with the input files described in step 1 (using the two-letter PDB directory structure convention).

Many commands are executed with `condapython`; this means they should be run with the python executable that comes with your conda distribution and includes MODELLER libraries. (Here, `condapython` is my alias to that executable).

### 1. Compile pipeline input

#### DRP file
Create a file with a list of DRPs, one per line. Each DRP should be represented by its PDB ID and
the chain in the PDB entry associated with the DRP. For example, for the DRP omega-grammotoxin SIA,
which is Chain A in PDB 1KOZ, add the line '1kozA' (no quotes).

#### PDB files
You will need local access to the PDB coordinate files that are listed in the DRP file. This can be structured in one of two ways:

todo: probably better formatting
1. Mirrored copy of the PDB - many institutions have a local mirrored copy of the PDB, using the [middle two character directory format](https://www.rcsb.org/pages/download/ftp). If all entries in the DRP list are accounted for here, you're good to go.

2. Alternatively, you can store the PDB entries for all DRPs in your input file in a single directory. They must be named with their standard PDB identifier (xyz case-insensitive?) (i.e. 1koz.pdb).

#### Example
No step is run here, but the PDB files need to be unpacked:

```
tar -xzf dividedPdb.tar.gz
```


todo - optimal column width?
todo - check 1koz case sensitive. Also syntax highlighting

### 2. Finalize PDB input
Run the setup_pdb.py script to extract the coordinates of the DRP chains in each PDB entry and write them out as a separate PDB file.

#### Example

```
condapython drp-clusters/drpclusters/setup_pdb.py  -q drpList.txt -p dividedPdbDir/
```
todo -- output dir too

### 3. Align DRP PDB files
The protocol creates pairwise distances matrices using two methods, Native Overlap and Equivalent Disulfides. These matrices must be prepared prior to running the full pipeline. These are ideally prepared on a distributed compute cluster as the computation time scales exponentially, but if there is a tractable number of DRPs, it's possible to use a single CPU  (for reference, 100 DRPs takes xyz on a xyz system). Scripts for both methods are described.

#### Example: Single processor
*This method iterates through all pairs of DRPs in `drpList.txt` and appends results to the `pairwise.txt` file*
```
condapython drp-clusters/drpclusters/pairwise_align.py  -q drpList.txt -p drpPdb -o pairwise.txt -m full_drp
condapython drp-clusters/drpclusters/pairwise_align.py  -q drpList.txt -p drpPdb -o disulfides.txt -m disulfides
grep longer_fraction pairwise.txt > longerFraction.txt
grep longer_sequence_product pairwise.txt > similarityProduct.txt
grep shorter_fraction pairwise.txt > shorterFraction.txt
```

todo -- decide whether to grep these into 'distances' dir

#### Example: Distributed system
*Alternatively, this method explicitly runs one command for each pair of DRPs. Each pair gets its own output file. To keep things clean, all output is stored in subdirectories. Since each pair of DRPs is processed with a single command, it can be run quickly on a cluster, in accordance with your cluster environment. These example commands only show all-vs-all commands for three DRPs; enumerating across all 100 pairs is left as an exercise to the user*

```
mkdir nativeOverlapWork
mkdir disulfideWork
condapython drp-clusters/drpclusters/align_native_overlap.py -f 1b45A -s 1dfnA -p drpPdb -o nativeOverlapWork/1b45A_1dfnA_no.txt
condapython drp-clusters/drpclusters/align_native_overlap.py -f 1b45A -s 1hjeA -p drpPdb -o nativeOverlapWork/1b45A_1dfnA_no.txt
condapython drp-clusters/drpclusters/align_native_overlap.py -f 1dfnA -s 1hjeA -p drpPdb -o nativeOverlapWork/1b45A_1dfnA_no.txt
condapython drp-clusters/drpclusters/align_native_overlap.py -f 1b45A -s 1dfnA -p drpPdb -o nativeOverlapWork/1b45A_1dfnA_no.txt
condapython drp-clusters/drpclusters/align_native_overlap.py -f 1b45A -s 1hjeA -p drpPdb -o nativeOverlapWork/1b45A_1dfnA_no.txt
condapython drp-clusters/drpclusters/align_native_overlap.py -f 1dfnA -s 1hjeA -p drpPdb -o nativeOverlapWork/1b45A_1dfnA_no.txt
```

*After all individual pairwise distance files have been generated, merge them into a final set of files with the following:*
```
grep longer_fraction nativeOverlapWork/*_no.txt > longerFraction.txt
grep longer_sequence_product nativeOverlapWork/*_no.txt > similarityProduct.txt
grep shorter_fraction nativeOverlapWork/*_no.txt > shorterFraction.txt
```
todo -- decide whether to grep these into 'distances' dir 

### 4. Run Cluster Pipeline

### 5. Create Visualization Sessions

### 6. Cluster Text annotation (coming soon)










```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
