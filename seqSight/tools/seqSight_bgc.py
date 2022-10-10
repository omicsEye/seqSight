import os
import sys
import re
import argparse
import logging
import shutil
from .. import utilities
from .. import config
from ..tools.seqSight_id import seqSight_reassign

logger = logging.getLogger(__name__)
def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(
        description="SeqSight to profile BGCs ",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i", "--fastq",
        required = True,
        help="reads file (fastq)",
    )
    parser.add_argument(
        "-o", "--output",
        required = True,
        help="Path for outputs",
    )
    parser.add_argument(
        "--output-basename",
        dest = 'output_basename',
        help="the basename for the output files\n[DEFAULT: " +
             "input file basename]",
        default=config.file_basename,
        metavar="<sample_name>")
    parser.add_argument(
        "-r","--resume",
        help="bypass commands if the output files exist\n",
        action="store_true",
        default=config.resume)
    parser.add_argument(
        "--threads",
        help="number of threads/processes\n[DEFAULT: " + str(config.threads) + "]",
        metavar="<" + str(config.threads) + ">",
        type=int,
        default=config.threads)
    parser.add_argument(
        "--one-contig",
        dest ="one_contig",
        help="If there is only one contig file for all samples, then this option should be provided",
        action="store_true")

    args = parser.parse_args()
    return args


def main():
    args = get_args()
    config.threads = args.threads
    config.resume = args.resume
    config.one_contig = args.one_contig

    # Set the basename of the output files if specified as an option
    if args.output_basename:
        config.file_basename = args.output_basename
    else:
        # Determine the basename of the input file to use as output file basename
        input_file_basename=os.path.basename(args.fastq)
        # Remove gzip extension if present
        if re.search('.gz$',input_file_basename):
            input_file_basename='.'.join(input_file_basename.split('.')[:-1])
        # Remove input file extension if present
        if '.' in input_file_basename:
            input_file_basename='.'.join(input_file_basename.split('.')[:-1])

        config.file_basename=input_file_basename
    config.output_folder = args.output
    config.temp_dir= config.output_folder+'/'+config.file_basename+'_temp'

    #Steps:

    #make a directory or outputs
    utilities.make_directory(config.output_folder)
    #utilities.make_directory(config.output_folder+'/no_hits')
    #utilities.make_directory(config.output_folder+'/hits')
    utilities.make_directory(config.temp_dir)
    # Steps

    #make a directory or outputs
    utilities.make_directory(config.temp_dir)
    temp_dir = config.temp_dir
    Input(R1.fastq+R2/R.fastq) –> map –> id –> qc -> .tsv/qc graph
    # map reads to reference database


    # reassignment model: quantifies abundance
    config.temp_dir = temp_dir+'/map_output/'
    utilities.make_directory(config.temp_dir)
    #seqSight_map -outDir TestData -indexDir /TestData -outAlignFile TestOutalign.sam -1 TestData/10R1.fastq -2 /Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/TestData/10R2.fastq -targetIndexPrefixes demo1_Testdata,demo2_Testdata
    seqSight_reassign(args.) #xinyang provide teh parameters

    # visualize QC results


    
    # make directory for bowtie2 output
    config.temp_dir = temp_dir+'/bowtie2_output/'
    utilities.make_directory(config.temp_dir)

    # make index database using bowtie2-build
    index_name = utilities.index(new_contig_file)

    # reads alignment using Bowtie2
    gene_alignment_file = utilities.alignment(args.fastq, index_name)

    # make directory for featureCounts abundance output
    config.temp_dir = temp_dir+'/featureCounts_output/'
    utilities.make_directory(config.temp_dir)

    # gene abundance using featureCounts
    config.temp_dir = temp_dir
    abundance_file = utilities.abundance(genes_file_gff, gene_alignment_file)
    shutil.move(abundance_file, config.output_folder+'/'+os.path.basename(os.path.normpath(abundance_file)))
    # if each sample has its own contig
    if not args.one_contig:
        # Run diamond
        uniref_alignment_file = utilities.diamond_alignment(genes_file_faa, args.uniref )

        # Infer abundance for sufficient hits to  uniref90 and no_hits
        hits, no_hits, uniref_gene_map = utilities.Infer_aligmnets(uniref_alignment_file, config.temp_dir)

        # select sequence for insufficient hits
        no_hits_genes_faa = utilities.select_sequnces(genes_file_faa, no_hits, output_name = '_no_hits.faa')

        # select sequence for sufficient hits
        hits_genes_faa = utilities.select_sequnces(genes_file_faa, hits, output_name = '_hits.faa')
    else:
        if not os.path.isfile(config.output_folder+'/no_hits/no_hits.txt') or \
                not os.path.isfile(config.output_folder+'/hits/hits.txt') or \
                not os.path.isfile(config.output_folder+'/hits/uniref_gene_map.txt'):
            # Run diamond
            uniref_alignment_file = utilities.diamond_alignment(genes_file_faa, args.uniref )

            # Infer abundance for sufficient hits to  uniref90 and no_hits
            hits, no_hits, uniref_gene_map = utilities.Infer_aligmnets(uniref_alignment_file, config.temp_dir)

            # select sequence for insufficient hits
            no_hits_genes_faa = utilities.select_sequnces(genes_file_faa, no_hits, output_name = '_no_hits.faa')

            # select sequence for sufficient hits
            hits_genes_faa = utilities.select_sequnces(genes_file_faa, hits, output_name = '_hits.faa')

            # move the three main output under main output folder from temp files
            # if there is more than one contig and the prodigal outputs ahvn't been produces (first sample run)

            shutil.move(no_hits_genes_faa, config.output_folder+'/no_hits/' + 'no_hits.faa')
            shutil.move(hits_genes_faa, config.output_folder+'/hits/' + 'hits.faa')
            shutil.move(no_hits, config.output_folder+'/no_hits/' + 'no_hits.txt')
            shutil.move(hits, config.output_folder+'/hits/' + 'hits.txt')
            shutil.move(uniref_gene_map, config.output_folder+'/hits/' + 'uniref_gene_map.txt')
            print ("Main output files for ppanini_press are written in: \n%s\n%s\n%s\n%s\n%s\n%s")% (config.output_folder+'/'+os.path.basename(os.path.normpath(abundance_file)),
                                                                                                     config.output_folder+'/no_hits/' + 'no_hits.faa',
                                                                                                     config.output_folder+'/no_hits/' + 'no_hits.txt',
                                                                                                     config.output_folder+'/hits/' + 'hits.faa',
                                                                                                     config.output_folder+'/hits/' + 'hits.txt',
                                                                                                     config.output_folder+'/hits/' + 'uniref_gene_map.txt')

    if not args.one_contig:
        # move the three main output under main output folder from temp files

        shutil.move(no_hits_genes_faa, config.output_folder+'/no_hits/'+ config.file_basename+ '_no_hits.faa')
        shutil.move(hits_genes_faa, config.output_folder+'/hits/'+ config.file_basename+ '_hits.faa')
        shutil.move(no_hits, config.output_folder+'/no_hits/'+ config.file_basename+ '_no_hits.txt')
        shutil.move(hits, config.output_folder+'/hits/'+ config.file_basename+ '_hits.txt')
        shutil.move(uniref_gene_map, config.output_folder+'/hits/'+ config.file_basename+ '_uniref_gene_map.txt')
        print ("Main output files for ppanini_press are written in: \n%s\n%s\n%s\n%s\n%s\n%s")% (config.output_folder+'/'+os.path.basename(os.path.normpath(abundance_file)),
                                                                                                 config.output_folder+'/no_hits/'+ config.file_basename+ '_no_hits.faa',
                                                                                                 config.output_folder+'/no_hits/'+ config.file_basename+ '_no_hits.txt',
                                                                                                 config.output_folder+'/hits/'+ config.file_basename+ '_hits.faa',
                                                                                                 config.output_folder+'/hits/'+ config.file_basename+ '_hits.txt',
                                                                                                 config.output_folder+'/hits/'+ config.file_basename+ '_uniref_gene_map.txt')


if __name__ == '__main__':
    main()