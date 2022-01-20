"""This script is used to convert a vcf file to a impute2 input file.

Author:  Shyam Gopalakrishnan
Date:    1st March 2019
Version: 1.0
"""


import argparse
import gzip
import os


progress = 10000


def parse_args():
    """Parse command line args.

    This function parses the input arguments and
    returns an object with them.

    Returns
    -------
    args : namespace
        Contains parsed command line arguments.

    """
    parser = argparse.ArgumentParser(description="VCF to impute2 converter")
    parser.add_argument("-i", "--inputFile", help="Input in vcf format",
                        required=True)
    parser.add_argument("-p", "--plColumn", type=int, help="Column of PLs", 
                        required=True)
    parser.add_argument("-o", "--outFile", help="Out filename", required=True)

    args = parser.parse_args()
    return(args)


def process_one_vcf_line(vcf_line, num_samples, pl_column):
    """Process one non-header line from vcf.

    Process one line of a vcf file (non-header version). Write information to
    legend file and update the strings for the haps file.

    Parameters
    ----------
    vcf_line : string
        Line from vcf file
    num_samples : int
        Number of samples
    pl_column : int
        Column in which the PLs are stored in a comma-seperated list

    Returns
    -------
    string
        output line for this vcf line - no newline present

    Raises
    ------
    ValueError
        when something is fishy about the vcf line.

    """
    vcf_toks = vcf_line.strip().split()
    header = vcf_toks[0] + "\t" + vcf_toks[2] + "\t"
    header += vcf_toks[1]+ "\t" + vcf_toks[3] + "\t" 
    header += vcf_toks[4] + "\t"
    if len(vcf_toks) != (num_samples + 9):
        raise ValueError("Something is rotten in the vcf line\n" + vcf_line)
    all_pls = []
    for token in vcf_toks[9:]:
        if token.count("./.") >= 1: ## missing PL.
            pls = [0, 0, 0]
        else:
            pl_token = token.split(":")[4]
            pls = [float(x) for x in pl_token.split(",")]
        pls = [10**(-x/10.0) for x in pls]
        sum_pls = sum(pls)
        pls = [x/sum_pls for x in pls]
        all_pls.extend(pls)
    return(header + "\t".join([str(round(x,5)) for x in all_pls]))

def main():
    """Do all the work here."""
    args = parse_args()
    if args.inputFile[-3:] == ".gz":
        infile = gzip.open(args.inputFile)
    else:
        infile = open(args.inputFile)
    # Skip the entire header section
    num_samples = 0
    for line in infile:
        toks = line.strip().split()
        if toks[0] == "#CHROM":
            num_samples = len(toks) - 9
            break
    
    # Open the output files
    if args.outFile[-3:] == ".gz":
        output_impute2 = gzip.open(args.outFile, "wb")
    else:
        output_impute2 = open(args.outFile, "w")

    args.plColumn -= 1
    # Now process one vcf line at a time
    cnt = 0
    outline = ""
    for line in infile:
        outline += process_one_vcf_line(line, num_samples, args.plColumn)
        outline += "\n"
        cnt = cnt + 1
        if cnt % progress == 0:
            #print cnt, "lines done."
            output_impute2.write(outline)
            outline = ""
    # check for leftovers
    if outline != "":
        output_impute2.write(outline)
        outline = ""
    output_impute2.close()
    infile.close()
    #print args.inputFile, "done processing."

main()
