import os
import sys
import re
import logging
import shutil
import argparse
import warnings


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help="files contains the sequences", type=str, required=True)
    parser.add_argument('--output', '-o', help="output directory to write teh results", type=str, required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    warnings.filterwarnings("ignore")
    warnings.simplefilter('ignore')

    args = parse_arguments()
    print(args)  # printing Namespace
    return print('done!')


# main()
if __name__ == "__main__":
    main()
