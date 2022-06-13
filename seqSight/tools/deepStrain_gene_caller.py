import os
import shutil

import utilities
import config


def main():
    # make directories
    utilities.make_directory(config.output_folder + '/no_hits')
    utilities.make_directory(config.output_folder + '/hits')

    temp_dir = config.temp_dir

    # get inputs reads in fastq file
    # new_contig_file = utilities.append_filename2contignames(args.contig)
    # rename the inputs to bowtie2_output folder
    os.system(
        "cp /Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/TestData/demo.fna  "
        "/Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/bowtie2_output/demo.fna")

    config.temp_dir = config.output_folder + '/prodigal_output/'
    utilities.make_directory(config.temp_dir)
    new_contig_file = "/Test/bowtie2_output/demo.fna"
    genes_file_gff, genes_file_fna, genes_file_faa = utilities.genecall(new_contig_file)

    # make bowtie2 output
    config.temp_dir = config.output_folder + '/bowtie2_output/'
    utilities.make_directory(config.temp_dir)
    index_name = utilities.index(new_contig_file)

    # reads alignment using Bowtie2
    gene_alignment_file = utilities.alignment(
        "/Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/TestData/demo1_Testdata.fasta", index_name)
    config.temp_dir = config.output_folder + '/featureCounts_output/'
    utilities.make_directory(config.temp_dir)
    # gene abundance using featureCounts
    print("start")
    abundance_file = utilities.abundance(genes_file_gff, gene_alignment_file)

    # os.system("featureCounts -T 8 -g ID -t CDS  -a
    # Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/prodigal_output/.gff -o
    # /Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/featureCounts_output/counts.txt
    # /Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/bowtie2_output/.sam")

    print("end")
    shutil.move(abundance_file, config.output_folder + '/' + os.path.basename(os.path.normpath(abundance_file)))

    # Run diamond alignment
    # TODO - FAILED - NO SPACE LEFT
    uniref_alignemnt_file = utilities.diamond_alignment(genes_file_faa, uniref_db)
    # Infer abundance for sufficient hits to  uniref90 and no_hits
    hits, no_hits, uniref_gene_map = utilities.Infer_aligmnets(uniref_alignemnt_file, config.temp_dir)

    # select sequence for insufficient hits
    no_hits_genes_faa = utilities.select_sequnces(genes_file_faa, no_hits, output_name='_no_hits.faa')

    # select sequence for sufficient hits
    hits_genes_faa = utilities.select_sequnces(genes_file_faa, hits, output_name='_hits.faa')


if __name__ == '__main__':
    uniref_db = "/Users/xinyangzhang/Documents/GitHub/seqSight/seqSight/Test/TestData/uniref90.fasta"
    main()
