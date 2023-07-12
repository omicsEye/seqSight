import os
import sys
import re
import logging
import shutil
import argparse
import warnings

import os, sys
seqSightdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(seqSightdir)
sys.path.insert(0, seqSightdir)

from seqSight.tools import seqSight_map, seqSight_id
class testseqSightMapOptions:
    MAX_REF_FILE_SIZE = 4.3e9
    verbose = False
    outDir = ""
    indexDir = ""
    numThreads = 8
    outAlignFile = "Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData/outAlign.sam"
    inReadFile = "/Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData/ex1.fastq"
    inReadFilePair1 = ""
    inReadFilePair2 = ""
    targetRefFiles = []
    filterRefFiles = []
    targetIndexPrefixes = ["Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData/genomes"]
    filterIndexPrefixes = []
    targetAlignFiles = []
    filterAlignFiles = []
    targetAlignParameters = None
    filterAlignParameters = None
    btHome = None
    exp_tag = ""


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help="files contains the sequences", type=str, required=True)
    parser.add_argument('--output', '-o', help="output directory to write teh results", type=str, required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    warnings.filterwarnings("ignore")
    warnings.simplefilter('ignore')

    # args = parse_arguments()
    # print(args)  # printing Namespace

    """
    step1 - MAP function;
    Input -  ex1.fastq;
    Output - ex1.sam;
    Tool - bowtie2;
    """
    outAlignFile_test = seqSight_map.processseqSightMap(testseqSightMapOptions)
    """
     step2 - ID function;
     Input -  sam file from the MAP function;
     Output - report.tsv;
     Tool - reassign+em;
     """
    seqSight_id.seqSight_reassign(True, 0.01, "testset", "sam", outAlignFile_test, "/Users/xinyang/Documents/Github/seqSight/seqSight/Test/TestData", 10,
                         not(False), 0, 0, False, False, emEpsilon=1e-7)
    # seqSight_id.seqSight_reassign(out_matrix, scoreCutoff, expTag, ali_format, ali_file, output, maxIter,
                                  # upalign, piPrior, thetaPrior, noCutOff, verbose, emEpsilon=0.01)

    return print('done!')




# main()
if __name__ == "__main__":
    main()
