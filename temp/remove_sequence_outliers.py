#!/usr/bin/env python3

# by Remi Allio, Damien de Vienne and Theo Tricou

import sys, subprocess, time, shlex, os, shutil, argparse, re, glob, subprocess
from shutil import copyfile
from argparse import RawTextHelpFormatter
from datetime import datetime

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Phylter is a tool that allows detecting, removing and visualizing outliers in phylogeneomics dataset by iteratively removing taxa in genes and optimizing a score of concordance between individual matrices.''', formatter_class=SmartFormatter)
    parser.add_argument('-j', '--job-name', help='''Job name also used as output directory''', default="PhylteR", dest='jobname')
    parser.add_argument('-s', '--sequences', help='''Directory containing FASTA gene sequences (aligned or not).''', dest='sequences')
    parser.add_argument('-o', '--outliers', help='''A 'phylter.out' file containing a list of outliers.''', dest='outliers')
    args = parser.parse_args()

    module_dir = os.path.dirname(__file__)
    module_dir = os.path.abspath(module_dir)


print("Starting sequences pruning step" + '\n')

phylter_out=open(args.outliers)
dico={}
phylter=0
for line in phylter_out:
	if not line.startswith("#"):
		line=line.rstrip()
		aln=line.split("\t")[0]
		sp=line.split("\t")[1]
		if aln in dico:
			dico[aln]=dico.get(aln)+',@'+sp+"@"
		else:
			dico[aln]="@"+sp+"@"
			phylter+=1

pathToOut = os.path.join(args.jobname, "sequences_PhylteR/")
if os.path.exists(pathToOut): shutil.rmtree(pathToOut)

os.makedirs(pathToOut)

if phylter >= 1:
    print(str(phylter) + " sequences files to filter" + '\n')
    list_ali=sorted(glob.glob(args.sequences+"/*"))
    same = []
    for ali in list_ali:
        match = next((x for x in list(dico.keys()) if x in ali), False)
        if match != False:
            ID=".".join(ali.split('/')[-1].split(".")[:-1])
            fout=open(os.path.join(pathToOut,ID+"_phylter.fasta"),"w")
            for name,seq in read_fasta(open(ali)):
                if not "@"+name.replace(">","")+"@" in dico.get(match):
                    fout.write(name+"\n"+seq+"\n")
            fout.close()
            if os.stat(os.path.join(pathToOut,ID+"_phylter.fasta")).st_size == 0:
                os.remove(os.path.join(pathToOut,ID+"_phylter.fasta"))
                print(" ".join(["The sequences file from gene", match, "is empty after outliers are removed"]))
        else:
            same.append(ali)
    if len(same) >= 1:
        print("Copying other sequences files"+'\n')
        for file in same:
            shutil.copy(file, pathToOut)
else:
    print("No sequences outliers detected by PhylteR to remove"+'\n')

tt=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print("Sequences have been removed - "+tt+'\n')
