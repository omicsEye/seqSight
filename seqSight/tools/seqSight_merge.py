import pandas as pd
import argparse
import csv
import sys
import os
from itertools import takewhile


def merge(aaastrIn, ostm, gtdb):
    listmpaVersion = set()
    profiles_list = []
    merged_tables = None

    for f in aaastrIn:
        headers = [x.strip() for x in takewhile(lambda x: x.startswith('#'), open(f))]
        listmpaVersion.add(headers[0])
        names = headers[-1].split('#')[1].strip().split('\t')

        iIn = pd.read_csv(f, sep='\t', skiprows=len(headers), names=names, usecols=[0, 2] if not gtdb else range(2),
                          index_col=0)
        profiles_list.append(pd.Series(data=iIn['relative_abundance'], index=iIn.index,
                                       name=os.path.splitext(os.path.basename(f))[0].replace('_profile', '')))

    merged_tables = pd.concat([merged_tables, pd.concat(profiles_list, axis=1).fillna(0)], axis=1).fillna(0)
    ostm.write(list(listmpaVersion)[0] + '\n')
    merged_tables.to_csv(ostm, sep='\t')

    return


argp = argparse.ArgumentParser(prog="merge_seqSight_tables.py",
                               description="Performs a table join on one or more seqSight output files.")
argp.add_argument("aistms", metavar="input.txt", nargs="*", help="One or more tab-delimited text tables to join")
argp.add_argument("-l", help="Name of file containing the paths to the files to combine")
argp.add_argument('-o', metavar="output.txt", help="Name of output file in which joined tables are saved")
argp.add_argument('--overwrite', default=False, action='store_true', help="Overwrite output file if exists")
argp.add_argument('--gtdb_profiles', action='store_true', default=False,
                  help=("To specify when running the script with GTDB-based profiles"))

argp.usage = (argp.format_usage() + "\nPlease make sure to supply file paths to the files to combine.\n\n" +
              "If combining 3 files (Table1.txt, Table2.txt, and Table3.txt) the call should be:\n" +
              "   ./merge_seqSight_tables.py Table1.txt Table2.txt Table3.txt > output.txt\n\n" +
              "A wildcard to indicate all .txt files that start with Table can be used as follows:\n" +
              "    ./merge_seqSight_tables.py Table*.txt > output.txt")


def main():
    args = argp.parse_args()

    if args.l:
        args.aistms = [x.strip().split()[0] for x in open(args.l)]

    if not args.aistms:
        print('merge_seqSight_tables: no inputs to merge!')
        return

    if args.o and os.path.exists(args.o) and not args.overwrite:
        print('merge_seqSight_tables: output file "{}" exists, specify the --overwrite param to ovrewrite it!'.format(
            args.o))
        return

    merge(args.aistms, open(args.o, 'w') if args.o else sys.stdout, args.gtdb_profiles)


if __name__ == '__main__':
    main()

    # # Get all file names from the result folder
    # file_path = "/Users/xinyang/Library/CloudStorage/Box-Box/Cancer_Mircrobiome/seqSight_BGC_All_Results/BGC_Frankel/"
    # path_list = os.listdir(file_path)
    # path_name = []
    # for i in path_list:
    #   path_name.append(i.split("-sam")[0])
    # print(sorted(path_name))
    #
    # # Set the first one as base
    # samp2 = pd.read_csv(sorted(path_name)[0], delimiter='\t', header=0, skiprows=1)
    # s2 = samp2[["Genome", "Final Guess"]]
    # # lst = ["Genome_ID", "SRR5930493"]
    # # s2.columns = lst
    #
    # # Apply for the rest others
    # for i in sorted(path_name)[1:]:
    #     if i != ".DS_Store":
    #         print(i)
    #         samp = pd.read_csv("seqSight_BGC_All_Results/BGC_Frankel/"+i+"-sam-report.tsv", delimiter='\t', header=0,skiprows=1)
    #         temp = samp[["Genome","Final Guess"]]
    #         renamelst =  ["Genome_ID",i]
    #         temp.columns = renamelst
    #         s2 = pd.merge(s2,temp,how='outer',on=['Genome_ID'])
    # # pd.set_option('display.max_rows',None)
    # print(s2)
    # s2.set_index('Genome_ID',inplace=True)
    # # s2['count'] = s2.apply(lambda x: x.count(), axis=1)
    # # s2.loc["count1"] = s2.count()
    #
    # # put the path you wish to save the output
    # s2.to_csv("seqSight_BGC_All_Results/BGC_Frankel_Merged.tsv")
