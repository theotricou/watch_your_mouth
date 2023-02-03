#!/usr/bin/env python3

# by Remi Allio, Damien de Vienne and Theo Tricou

import subprocess
from subprocess import Popen
import sys, subprocess, time, shlex, os, shutil
from shutil import copyfile
import glob
import collections
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import random
import string
import tempfile
from threading import Thread
import re
from argparse import RawTextHelpFormatter
import argparse, os, shlex, shutil, sys
import os.path
from argparse import RawTextHelpFormatter
from datetime import datetime
import time
import operator

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
    parser.add_argument('-t', '--trees', help='''Directory containing gene tree files (only)''', default="", dest='trees')
    parser.add_argument('-s', '--sequences', help='''Directory containing FASTA gene sequences (aligned or not). Optional but mandatory for cleaning step. FASTA files must be .fasta, .fst, .faa or .fna files.''', default="NULL", dest='sequences')
    parser.add_argument('-bv', '--bvalue', help='''Nodes with a support below this value will be collapsed prior to the outlier detection''', default=0, type=float, dest='bvalue')
    parser.add_argument('-d', '--distance', help = 'Type of distance used to compute the pairwise matrices for each tree. Can be \"patristic\" (sum of branch lengths separating tips, the default) or nodal (number of nodes separating tips).The \"nodal\" option should only be used if all species are present in all genes.', default="patristic", dest='distance')
    parser.add_argument('-k', help='''Strength of outlier detection. The higher this value the less outliers detected. (Default=3)''', default=3, type=int, dest='k')
    parser.add_argument('-k2', help='''Same as k for complete gene outlier detection. To preserve complete genes from being discarded, k2 can be increased . By default, k2 = k. (see above)''', default=3, type=int, dest='k2')
    parser.add_argument('--norm', help='''Should the matrices be normalized prior to the complete analysis and how. If \"median\" (default), matrices are divided by their median, if \"mean\" they are divided by their mean, if \"none\", no normalization if performed. Normalizing ensures that fast-evolving (and slow-evolving) genes are not treated as outliers. Normalization by median is a better choice as it is less sensitive to outlier values.''', default="median", dest='norm')
    parser.add_argument('--norm-cutoff', help='''Value of the median (if Norm="median") or the mean (if Norm="mean") of phylogenetic distances in the matrix below which a gene is discarded from the analysis. This prevents dividing by 0, and allows getting rid of genes that contain mostly branches of length 0 and are therefore uninformative anyway. Discarded genes, if any, are listed in the output (out$DiscardedGenes). Default to 0.000001''', default=0.000001, type=float, dest='normcutoff')
    parser.add_argument('--test-island', help='''If TRUE (the default), only the highest value in an ’island’ of outliers is consid- ered an outlier. This prevents non-outliers hitchhiked by outliers to be considered outliers themselves.''', default="TRUE", dest='island')
    parser.add_argument('--verbose', help='''If TRUE (the default), messages are written during the filtering process to get information on what is happening.''', default="TRUE", dest="verbose")
    parser.add_argument('--stop-criteria', help='''The optimization stops when the gain (quality of compromise) between round n and round n+1 is smaller than this value. Default to 0.00001.''', default=0.00001, type=float, dest="stop")
    parser.add_argument('--initial-only', help='''If TRUE, only the Initial state of the data is computed.''', default="FALSE", dest="InitialOnly")
    parser.add_argument('--normalizeby', help='''Should the gene x species matrix be normalized prior to outlier detection. By \"row\" (=species; default) or by \"col\" (=gene)''', default="row", dest='normalizeby')
    parser.add_argument('--parallel', help='''Logical. Should the computations be parallelized when possible? Default to TRUE. Note that the number of threads cannot be set by the user when parallel = TRUE. It uses all available cores on the machine.''', default="TRUE", dest="parallel")
    parser.add_argument('-v', '--version', help=''''Version 0.9.6''', default=False, dest='versionCheck', action='store_true')
    parser.add_argument("-g", '--graphical', help=''''Logical. Should a full report of the phylter analysis be output in a pdf file''', default="FALSE", dest='graph')
    parser.add_argument("-p",'--prune', help=''''Logical. Should gene trees be pruned of outliers detected by PhylteR''', default="NULL", dest='tprune')
    # to be completed
    parser.add_argument('--example', help='''Print getting started examples''', default=False, dest='example', action='store_true')
    # to be completed
    parser.add_argument('--citation', help='''How to cite PhylteR''', default=False, dest='citation', action='store_true')
    args = parser.parse_args()

    module_dir = os.path.dirname(__file__)
    module_dir = os.path.abspath(module_dir)

    if args.example == True:
        print("\n # One basic example:\nPhylteR --trees [.tre directory \
        containing all newick formatted tree files] \\\n\t ")
        exit()
    if args.citation == True:
        txt="\n\nIf you use PhylteR, please cite:\n-TO BE UPDATED.\n\n"
        txt=txt.replace(u"\u2013","-")
        print(txt)
        exit()

    if args.versionCheck == True:
        command = "Rscript -e 'library(pylter); packageVersion(phylter)'"
        output = subprocess.check_output(command, shell = True)
        print(paste("PhylteR current version on this singularity container is: ", re.split('‘|’', output.decode('utf-8').strip())[1], sep =""))
        exit()

    start_time=datetime.now()
    initial_path= os.getcwd()+"/"
    pathToWork = os.getcwd()+"/"+args.jobname+"/"

    Logfile = os.path.join(initial_path,args.jobname+"_PhylteR.log")

    args.genes="NULL"

    logfile=open(Logfile,"a")

    if os.path.isdir(args.jobname):
        print("\nERROR: Job name is already being used."+"\n")
        logfile.write("\nERROR: Job name is already being used."+"\n"+"\n")
        exit()

    if args.jobname == "":
        print("\nERROR: Job name is required (-j option)."+"\n")
        logfile.write("\nERROR: Job name is required (-j option)."+"\n"+"\n")
        exit()

    if not os.path.exists(pathToWork):
        os.makedirs(pathToWork)
        # os.chdir(pathToWork) #to be updated

    if args.trees == "":
        print("\nERROR: Tree files directory is required (--trees option)."+"\n")
        logfile.write("\nERROR: Tree files directory is required (--trees option).\n"+"\n")
        exit()

    if args.trees != "":
        if not os.path.isdir(os.path.join(initial_path,args.trees)):
            print("ERROR: "+str(os.path.join(initial_path,args.trees)+ " no such directory.\nAborting.\n\n"))
            logfile.write("ERROR: "+str(os.path.join(initial_path,args.trees)+ " no such directory.\nAborting.\n\n"))
            exit()

    if args.sequences != "NULL":
        if not os.path.isdir(os.path.join(initial_path,args.sequences)):
            print("ERROR: "+str(os.path.join(initial_path,args.sequences)+ " no such directory.\nAborting.\n\n"))
            logfile.write("ERROR: "+str(os.path.join(initial_path,args.sequences)+ " no such directory.\nAborting.\n\n"))
            exit()

    if args.normalizeby != "row" and args.normalizeby != "col":
        print("\nERROR: detection mode is unclear (--detection-mode option). Need to be \"row\" or \"col\""+"\n")
        logfile.write("\nERROR: detection mode is unclear (--detection-mode option). Need to be \"row\" or \"col\""+"\n"+"\n")
        exit()


    print('\nStarting PhylteR ')
    logfile.write('Starting PhylteR '+ "\n"+"\n")
    tt=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Start time : "+tt)
    logfile.write("Start time : "+tt+"\n"+"\n")
    print('\nResult files will be saved here: ')
    print(str(pathToWork))
    logfile.write('\nResult files will be saved here: '+"\n"+str(pathToWork)+"\n"+"\n")




    ##################################################
    ######        export tree files            #######
    ##################################################

    if os.path.isdir(os.path.join(initial_path,args.trees)):
        print("\nThe following directory was provided to collect the trees: "+str(os.path.join(initial_path,args.trees))+"\n")
        logfile.write("\nThe following was provided to collect the trees: "+str(os.path.join(initial_path,args.trees))+"\n"+"\n")
        trees_dir = os.path.join(initial_path,args.trees)
        all_tree_output = os.path.join(args.jobname, "PhylteR_all_tree_named")
        lt = glob.glob(trees_dir+"/*")
        lst = [os.path.splitext(x)[0] for x in lt]
        # test if the extension from the name of the tree file can be removed
        if len(lt) != len(set(lst)):
            command = 'for i in ' + trees_dir + '/*; do echo `echo $i | rev | cut -d'/' -f 1 | rev``cat $i`; done > ' + all_tree_output
        else:
            command = 'for i in ' + trees_dir + "/*; do echo `echo $i | rev | cut -d'.' -f 2 | cut -d'/' -f 1 | rev``cat $i`; done > " + all_tree_output
        print("Creating a MultiPhylo file containing names"+"\n")
        logfile.write("Creating a MultiPhylo file containing names"+"\n"+"\n")
        concat_tree = Popen(command, stdout=logfile, stderr=logfile, shell=True)
        concat_tree.wait()

        if not os.path.isfile(all_tree_output):
            print("ERROR: No trees were found in "+str(os.path.join(initial_path,args.trees)+ ".\nThe files should ideally be .nwk, .tre or .treefile files.\nAborting.\n\n"))
            logfile.write("ERROR: No trees were found in "+str(os.path.join(initial_path,args.trees)+ ".\nThe files should ideally be .nwk, .tre or .treefile files.\nAborting.\n\n"))
            exit()
        trees=all_tree_output
    else:
        print("ERROR: "+str(os.path.join(initial_path,args.trees)+ " no such directory.\nAborting.\n\n"))
        logfile.write("ERROR: "+str(os.path.join(initial_path,args.trees)+ " no such directory.\nAborting.\n\n"))
        exit()

    #####################################
    ######       Run PhylteR       ######
    #####################################

    command = 'Rscript '+ module_dir + '/PhylteR_run.R ' + str(trees) + ' ' + str(args.bvalue) + ' ' + str(args.distance) + ' ' + str(args.k) + ' ' + str(args.k2) + ' ' + str(args.norm) + ' ' + str(args.normcutoff) + ' ' + args.genes + ' ' + str(args.island) + ' ' + str(args.verbose) + ' ' + str(args.stop) + ' ' + str(args.InitialOnly) + ' ' + str(args.normalizeby) + ' ' + str(args.parallel) + ' ' + str(args.graph) + ' ' + str(args.jobname)


    print("PhylteR command used = "+command+"\n")
    logfile.write("PhylteR command used = "+command +"\n"+"\n")
    phylter = Popen(command, stdout= subprocess.PIPE, stderr = logfile, shell=True)
    for line in phylter.stdout:
        print(line.decode('utf-8'))
        logfile.write(line.decode('utf-8'))
    phylter.wait()

    if os.path.isfile(os.path.join(pathToWork,"phylter.out")):
        print("PhylteR is done - "+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'\n')
        logfile.write("PhylteR is done - "+datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'\n'+"\n")
    else:
        print("An error occured. Please see "+Logfile)
        logfile.write("An error occured. Please see "+Logfile+"\n"+"\n")
        exit()

    #####################################
    ######       Prune Trees       ######
    #####################################

    if args.tprune != "NULL":
        print("\nStarting tree pruning step" + '\n')
        logfile.write("\nStarting tree pruning step" + '\n'+"\n")

        command = 'prune_outliers.R ' + str(args.trees) + ' ' + str(args.jobname)
        print("Tree pruning command used = "+command+"\n")
        logfile.write("Tree pruning command used = "+command +"\n"+"\n")

        prune_trees = Popen(command, stdout= subprocess.PIPE, stderr = logfile, shell=True)
        for line in prune_trees.stdout:
            print(line.decode('utf-8'))
            logfile.write(line.decode('utf-8'))
        prune_trees.wait()

        tt=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("Trees have been pruned - "+tt+'\n')
        logfile.write("Trees have been pruned - "+tt+'\n'+"\n")

    #####################################
    ######     Prune sequences     ######
    #####################################

    if args.sequences != "NULL":
        print("Starting sequences pruning step" + '\n')
        logfile.write("Starting sequences pruning step"+'\n'+"\n")

        phylter_out=open(os.path.join(args.jobname, "phylter.out"))
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

        pathToOut = os.path.join(args.jobname, "seqs_PhylteR/")
        if os.path.exists(pathToOut): shutil.rmtree(pathToOut)

        os.makedirs(pathToOut)

        if phylter >= 1:
            print(str(phylter) + " sequences files to filter" + '\n')
            logfile.write(str(phylter) + " sequences files to filter" + '\n'+"\n")

            aln_dir="alis"
            list_ali=sorted(glob.glob(aln_dir+"/*"))
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
                        print(" ".join(["The sequences file from gene", match, "is empty after outliers are removed"])+'\n')
                        logfile.write(" ".join(["The sequences file from gene", match, "is empty after outliers are removed"])+'\n'+"\n")

                else:
                    same.append(ali)

            if len(same) >= 1:
                print("Copying other sequences files"+'\n')
                logfile.write("Copying other sequences files"+'\n'+"\n")
                for file in same:
                    shutil.copy(file, pathToOut)
                # command = 'cp {' + ",".join(same) + '} ' + pathToOut
                # phylter = Popen(command, stdout=logfile, stderr=logfile, shell=True)
        else:
            print("No sequences outliers detected by PhylteR to remove"+'\n')
            logfile.write("No sequences outliers detected by PhylteR to remove"+'\n'+"\n")

        tt=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print("Sequences have been removed - "+tt+'\n')
        logfile.write("Sequences have been removed - "+tt+'\n'+"\n")

    tt=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Everything is done - "+tt+'\n')
    logfile.write("Everything is done - "+tt+'\n'+"\n")
    logfile.close()
