import os, shutil, math, sys, argparse;

# TODO - Need to transfer to seqSight defalut directory
# pathoscopedir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.insert(0,pathoscopedir)


from seqSight.tools.bowtie2Wrap import bowtie2Wrap
from seqSight.tools.utils import seqParse


def parse_arguments(args):
    """
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description="Map the input files(single or paired-end files) to the Reference Database \n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--MAX_REF_FILE_SIZE",
        help="MAX_REF_FILE_SIZE\n",
        default="4.3e9")
    parser.add_argument(
        "-v","--verbose",
        default="False"),
    parser.add_argument(
        "--outDir",
        default=".",
        help="")
    parser.add_argument(
        "--indexDir",
        default=".",
        help="")
    parser.add_argument(
        "--numThreads",
        default="8",
        help="Number of the threads.\n", )
    parser.add_argument(
        "--outAlignFile",
        default="outalign.sam")
    parser.add_argument(
        "-U", "--inReadFile",
        help="The single input file from your single file\n",
        default="",
        required=False)
    parser.add_argument(
        "-1", "--inReadFilePair1",
        default="",
        help="The R1 input file from your paired-end file\n",
        required=False)
    parser.add_argument(
        "-2", "--inReadFilePair2",
        default="",
        help="The R2 input file from your paired-end file\n",
        required=False)
    parser.add_argument(
        "--targetRefFiles",
        default="",
        help="The target reference file did not ends in bowtie2 format\n",
        required=False)
    parser.add_argument(
        "--filterRefFiles",
        default="",
        help="The filter reference file did not ends in bowtie2 format\n",
        required=False)
    parser.add_argument(
        "--targetIndexPrefixes",
        help="The target reference file ends in bowtie2 format and put their prefixes\n",
        default="",
        required=False)
    parser.add_argument(
        "--filterIndexPrefixes",
        default="",
        help="The filter reference file ends in bowtie2 format and put their prefixes\n",
        required=False)
    parser.add_argument(
        "--targetAlignFiles",
        help="",
        required=False)
    parser.add_argument(
        "--filterAlignFiles",
        help="",
        required=False)
    parser.add_argument(
        "--targetAlignParameters",
        default="None")
    parser.add_argument(
        "--filterAlignParameters",
        default="None")
    parser.add_argument(
        "--btHome",
        default="None")
    parser.add_argument(
        "--exp_tag",
        default="",
        required=False)

    return parser.parse_args()


def main():
    # Make a copy of core seqSightOptions
    def copyseqSightMapOptions(seqSightMapOptions):
        procseqSightMapOptions = seqSightMapOptions()
        procseqSightMapOptions.verbose = seqSightMapOptions.verbose
        procseqSightMapOptions.outDir = seqSightMapOptions.outDir
        procseqSightMapOptions.indexDir = seqSightMapOptions.indexDir
        procseqSightMapOptions.numThreads = seqSightMapOptions.numThreads
        procseqSightMapOptions.outAlignFile = seqSightMapOptions.outAlignFile
        procseqSightMapOptions.inReadFile = seqSightMapOptions.inReadFile
        procseqSightMapOptions.inReadFilePair1 = seqSightMapOptions.inReadFilePair1
        procseqSightMapOptions.inReadFilePair2 = seqSightMapOptions.inReadFilePair2
        procseqSightMapOptions.targetRefFiles = seqSightMapOptions.targetRefFiles
        procseqSightMapOptions.filterRefFiles = seqSightMapOptions.filterRefFiles
        procseqSightMapOptions.targetIndexPrefixes = seqSightMapOptions.targetIndexPrefixes
        procseqSightMapOptions.filterIndexPrefixes = seqSightMapOptions.filterIndexPrefixes
        procseqSightMapOptions.targetAlignFiles = seqSightMapOptions.targetAlignFiles
        procseqSightMapOptions.filterAlignFiles = seqSightMapOptions.filterAlignFiles
        procseqSightMapOptions.targetAlignParameters = seqSightMapOptions.targetAlignParameters
        procseqSightMapOptions.filterAlignParameters = seqSightMapOptions.filterAlignParameters
        procseqSightMapOptions.btHome = seqSightMapOptions.btHome
        procseqSightMapOptions.exp_tag = seqSightMapOptions.exp_tag
        return procseqSightMapOptions

    # If the given file size is greater than maxSize, then it splits
    # and returns a list of split file paths where each file is less than maxSize
    def splitCheck(filePath, maxSize):
        files = []
        fileSize = os.stat(filePath).st_size
        nSplit = 1
        if (fileSize > maxSize):
            nSplit = int(math.ceil(1.0 * fileSize / float(maxSize)))
        if nSplit == 1:
            files.append(filePath)
            return files
        (base, ext) = os.path.splitext(filePath)
        # check if we have already done this splitting
        for i in range(nSplit):
            fiPath = base + '_' + str(i) + ext
            splitReq = False
            if not os.path.exists(fiPath):
                splitReq = True
                break
        fps = []
        for i in range(nSplit):
            fiPath = base + '_' + str(i) + ext
            files.append(fiPath)
            if splitReq:
                fps.append(open(fiPath, 'w'))
        if splitReq:
            with open(filePath, 'r') as fp:
                j = 0
                if ext == '.fq':
                    for r in seqParse.parse(fp, 'fastq'):
                        fps[j % nSplit].write('>%s %s\n%s\n%s\n' % (r.id, r.description, r.seq, r.qual))
                        j += 1
                else:
                    for r in seqParse.parse(fp, 'fasta'):
                        fps[j % nSplit].write('>%s %s\n%s\n' % (r.id, r.description, r.seq))
                        j += 1
            for i in range(nSplit):
                fps[i].close()
        return files

    def filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile):
        return bowtie2Wrap.filter_alignment(targetAlignFile, filterAlignFiles, outAlignFile)

    # Appends all the appendFiles to outfile
    def append_file(outfile, appendFiles):
        with open(outfile, 'w') as out1:
            for file1 in appendFiles:
                if (file1 is not None):
                    with open(file1, 'r') as in2:
                        for ln in in2:
                            out1.write(ln)

    # Appends all the sam appendFiles to outfile
    def append_sam_file(outfile, appendFiles):
        with open(outfile, 'w') as out1:
            # First, writing the header by merging headers from all files
            for file1 in appendFiles:
                if (file1 is not None):
                    with open(file1, 'r') as in2:
                        for ln in in2:
                            if ln[0] == '@':
                                out1.write(ln)
            # Writing the body by merging body from all files
            for file1 in appendFiles:
                if (file1 is not None):
                    with open(file1, 'r') as in2:
                        for ln in in2:
                            if ln[0] != '@':
                                out1.write(ln)

    args = parse_arguments(sys.argv)

    # ===========================================================
    class seqSightMapOptions:
        MAX_REF_FILE_SIZE = float(args.MAX_REF_FILE_SIZE)
        verbose = args.verbose
        outDir = args.outDir
        indexDir = args.indexDir
        numThreads = args.numThreads
        outAlignFile = args.outAlignFile
        inReadFile = args.inReadFile
        inReadFilePair1 = args.inReadFilePair1
        inReadFilePair2 = args.inReadFilePair2
        targetRefFiles = args.targetRefFiles
        filterRefFiles = args.filterRefFiles
        targetIndexPrefixes = [args.targetIndexPrefixes]
        filterIndexPrefixes = [args.filterIndexPrefixes]
        targetAlignFiles = [args.targetAlignFiles]
        filterAlignFiles = [args.filterAlignFiles]
        targetAlignParameters = args.targetAlignParameters
        filterAlignParameters = args.filterAlignParameters
        btHome = args.btHome
        exp_tag = args.exp_tag

    # Main entry function to seqSightMap that does all the processing
    # def processseqSightMap(seqSightMapOptions):
    procseqSightMapOptions = copyseqSightMapOptions(seqSightMapOptions)

    print("targetAlignFiles", seqSightMapOptions.targetAlignFiles)

    # Splitting reference files if bigger than MAX_REF_FILE_SIZE
    ptargetRefFiles = []
    if len(seqSightMapOptions.targetRefFiles) > 1:
        for filePath in seqSightMapOptions.targetRefFiles.split(','):
            if seqSightMapOptions.verbose:
                print("Checking whether the file: " + filePath + " needs to be split")
            files = splitCheck(filePath, seqSightMapOptions.MAX_REF_FILE_SIZE)
            for f in files:
                ptargetRefFiles.append(f)
    elif len(seqSightMapOptions.targetRefFiles) == 1:
        for filePath in seqSightMapOptions.targetRefFiles:
            if seqSightMapOptions.verbose:
                print("Checking whether the file: " + filePath + " needs to be split")
            files = splitCheck(filePath, seqSightMapOptions.MAX_REF_FILE_SIZE)
            for f in files:
                ptargetRefFiles.append(f)

    procseqSightMapOptions.targetRefFiles = ptargetRefFiles
    print("###############", procseqSightMapOptions.targetRefFiles)

    pfilterRefFiles = []
    if len(seqSightMapOptions.filterRefFiles) > 1:
        for filePath in seqSightMapOptions.filterRefFiles.split(','):
            print(filePath)
            if seqSightMapOptions.verbose:
                print("Checking whether the file: " + filePath + " needs to be split")
            files = splitCheck(filePath, seqSightMapOptions.MAX_REF_FILE_SIZE)
            for f in files:
                pfilterRefFiles.append(f)
    elif len(seqSightMapOptions.filterRefFiles) == 1:
        for filePath in seqSightMapOptions.filterRefFiles:
            print(filePath)
            if seqSightMapOptions.verbose:
                print("Checking whether the file: " + filePath + " needs to be split")
            files = splitCheck(filePath, seqSightMapOptions.MAX_REF_FILE_SIZE)
            for f in files:
                pfilterRefFiles.append(f)
        procseqSightMapOptions.filterRefFiles = pfilterRefFiles

    # Creating Index if it does not exist
    bowtie2Options = bowtie2Wrap.Bowtie2Options()
    bowtie2Options.verbose = procseqSightMapOptions.verbose
    bowtie2Options.btHome = procseqSightMapOptions.btHome
    bowtie2Options.indexDir = procseqSightMapOptions.indexDir
    print("ptargetRefFiles", ptargetRefFiles)
    """
    Clean the procseqSightMapOptions.targetIndexPrefixes list to be empty, remove the default value there 
    If the input has procseqSightMapOptions.targetIndexPrefixes value then need to add the case
    """
    if procseqSightMapOptions.targetIndexPrefixes[0] == "":
        procseqSightMapOptions.targetIndexPrefixes = []
    for filePath in ptargetRefFiles:
        bowtie2Options.refFile = filePath
        (_, tail) = os.path.split(filePath)
        (base, _) = os.path.splitext(tail)
        bowtie2Options.btIndexPrefix = base
        if seqSightMapOptions.verbose:
            print("Creating bowtie2 index for: " + filePath)
        bowtie2Wrap.create_bowtie2_index(bowtie2Options)
        print("base", base)
        procseqSightMapOptions.targetIndexPrefixes.append(base)
        # print("procseqSightMapOptions.targetIndexPrefixes", procseqSightMapOptions.targetIndexPrefixes)
    for filePath in pfilterRefFiles:
        bowtie2Options.refFile = filePath
        (_, tail) = os.path.split(filePath)
        (base, _) = os.path.splitext(tail)
        bowtie2Options.btIndexPrefix = base
        if seqSightMapOptions.verbose:
            print("Creating bowtie2 index for: " + filePath)
        bowtie2Wrap.create_bowtie2_index(bowtie2Options)
        procseqSightMapOptions.filterIndexPrefixes.append(base)

    # Creating the Alignment file
    bowtie2Options = bowtie2Wrap.Bowtie2Options()
    bowtie2Options.verbose = procseqSightMapOptions.verbose
    bowtie2Options.btHome = procseqSightMapOptions.btHome
    bowtie2Options.numThreads = procseqSightMapOptions.numThreads
    bowtie2Options.outDir = procseqSightMapOptions.outDir
    bowtie2Options.indexDir = procseqSightMapOptions.indexDir
    bowtie2Options.readFile = procseqSightMapOptions.inReadFile
    bowtie2Options.readFilePair1 = procseqSightMapOptions.inReadFilePair1
    bowtie2Options.readFilePair2 = procseqSightMapOptions.inReadFilePair2

    print("procseqSightMapOptions.targetIndexPrefixes", procseqSightMapOptions.targetIndexPrefixes)

    if (len(procseqSightMapOptions.inReadFilePair1) > 0 and
            len(procseqSightMapOptions.inReadFilePair2) > 0 and
            len(procseqSightMapOptions.inReadFile) > 0):
        bowtie2Options.bothReadFlag = True  # newly added
    elif (len(procseqSightMapOptions.inReadFilePair1) > 0 and
          len(procseqSightMapOptions.inReadFilePair2) > 0):
        bowtie2Options.pairedReadFlag = True
    elif len(procseqSightMapOptions.inReadFile) > 0:
        bowtie2Options.singleReadFlag = True  # newly added
    if procseqSightMapOptions.targetAlignParameters is not None:
        bowtie2Options.additionalOptions = procseqSightMapOptions.targetAlignParameters

    '''
    Two situations for the length of the targetIndexPrefiexes
    1. procseqSightMapOptions.targetIndexPrefixes ['demo1_Testdata']
    2. procseqSightMapOptions.targetIndexPrefixes ['demo1_Testdata,demo2_Testdata....']
    '''

    print("11111111111procseqSightMapOptions.targetAlignFiles", procseqSightMapOptions.targetAlignFiles)
    print(len(procseqSightMapOptions.targetAlignFiles))
    if len(procseqSightMapOptions.targetIndexPrefixes[0].split(',')) > 1:
        for indexPrefix in procseqSightMapOptions.targetIndexPrefixes[0].split(','):
            bowtie2Options.btIndexPrefix = procseqSightMapOptions.indexDir + os.sep + indexPrefix
            bowtie2Options.outAlignFile = procseqSightMapOptions.exp_tag + indexPrefix + ".sam"
            print("bowtie2Options.outAlignFile", bowtie2Options.outAlignFile)
            if seqSightMapOptions.verbose:
                print("Creating bowtie2 alignment: " + bowtie2Options.outAlignFile)
            bowtie2Wrap.run_bowtie2(bowtie2Options)
            procseqSightMapOptions.targetAlignFiles.append(procseqSightMapOptions.outDir + os.sep +
                                                           bowtie2Options.outAlignFile)
    if len(procseqSightMapOptions.targetIndexPrefixes[0].split(',')) == 1:
        for indexPrefix in procseqSightMapOptions.targetIndexPrefixes:
            bowtie2Options.btIndexPrefix = procseqSightMapOptions.indexDir + os.sep + indexPrefix
            bowtie2Options.outAlignFile = procseqSightMapOptions.exp_tag + indexPrefix + ".sam"
            print("bowtie2Options.outAlignFile", bowtie2Options.outAlignFile)
            if seqSightMapOptions.verbose:
                print("Creating bowtie2 alignment: " + bowtie2Options.outAlignFile)
            bowtie2Wrap.run_bowtie2(bowtie2Options)
            procseqSightMapOptions.targetAlignFiles.append(procseqSightMapOptions.outDir + os.sep +
                                                           bowtie2Options.outAlignFile)

    print("222222222procseqSightMapOptions.targetAlignFiles", procseqSightMapOptions.targetAlignFiles)

    # Appending the Alignment files and Filtering
    if len(procseqSightMapOptions.targetAlignFiles) > 1:
        appendAlignFile = procseqSightMapOptions.outDir + os.sep + procseqSightMapOptions.exp_tag + "appendAlign.sam"
        if seqSightMapOptions.verbose:
            print("Appending alignment files to: " + appendAlignFile)
        print("procseqSightMapOptions.targetAlignFiles", procseqSightMapOptions.targetAlignFiles)
        append_sam_file(appendAlignFile, procseqSightMapOptions.targetAlignFiles)
    #     TODO Add a case for there is only one target
    else:
        appendAlignFile = procseqSightMapOptions.targetAlignFiles[1]

    if len(procseqSightMapOptions.filterIndexPrefixes) > 0:
        bowtie2Options.readFile = procseqSightMapOptions.outDir + os.sep + procseqSightMapOptions.exp_tag + "appendAlign.fq"
        bowtie2Options.readFilePair1 = ""
        bowtie2Options.readFilePair2 = ""
        bowtie2Options.bothReadFlag = False
        bowtie2Options.pairedReadFlag = False
        bowtie2Options.singleReadFlag = True
        if procseqSightMapOptions.filterAlignParameters is not None:
            bowtie2Options.additionalOptions = procseqSightMapOptions.filterAlignParameters
        bowtie2Wrap.extractRead(appendAlignFile, bowtie2Options.readFile)
        if len(procseqSightMapOptions.filterIndexPrefixes) > 1:
            for indexPrefix in procseqSightMapOptions.filterIndexPrefixes:
                print("indexPrefix", indexPrefix)
                bowtie2Options.btIndexPrefix = procseqSightMapOptions.indexDir + os.sep + indexPrefix
                bowtie2Options.outAlignFile = procseqSightMapOptions.exp_tag + indexPrefix + ".sam"
                if seqSightMapOptions.verbose:
                    print("Creating bowtie2 alignment: " + bowtie2Options.outAlignFile)
                bowtie2Wrap.run_bowtie2(bowtie2Options)
                procseqSightMapOptions.filterAlignFiles.append(procseqSightMapOptions.outDir + os.sep +
                                                              bowtie2Options.outAlignFile)
        #TODO elif  len(procseqSightMapOptions.filterIndexPrefixes) == 1:

    # Filtering the Alignment file
    outAlignFile = procseqSightMapOptions.outDir + os.sep + procseqSightMapOptions.outAlignFile
    if seqSightMapOptions.verbose:
        print("Filtering and creating the alignment: " + outAlignFile)
    if len(procseqSightMapOptions.filterAlignFiles) > 0:
        filter_alignment(appendAlignFile, procseqSightMapOptions.filterAlignFiles, outAlignFile)
    elif ((len(procseqSightMapOptions.targetAlignFiles) > 1) or \
          (len(procseqSightMapOptions.targetIndexPrefixes) > 0)):
        os.rename(appendAlignFile, outAlignFile)
    else:  # Input appendAlignFile provided by user, hence make a copy for outAlignFile
        shutil.copy(appendAlignFile, outAlignFile)


if __name__ == '__main__':
    main()
