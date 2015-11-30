#!/usr/bin/python
from pyfaidx import Fasta
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess

def return_bioseq_obj(sequence):
    """ Returns biopython sequence object
    """
    return Seq(sequence,'generic_dna')

#### LOAD THE REFERENCE FILE #### 

mzebra_ref = Fasta(sys.argv[1])
inversion_file = sys.argv[2]
blast_database = sys.argv[3]

with open(inversion_file,'r') as INV_FILE:
    for line in INV_FILE:
        line = line.rstrip('\n')
        scaffold = line.split('\t')[0]
        start = line.split('\t')[1]
        stop = line.split('\t')[2]
        print (scaffold)
        inversion = mzebra_ref[scaffold][int(start):int(stop)].seq
        inversion_seq = return_bioseq_obj(inversion)
        with open('temp.fasta','w') as TEMP:
            TEMP.write(">%s_%s_%s\n"%(scaffold,start,stop))
            TEMP.write(inversion)
            out_file = scaffold+"_"+start+"_"+stop 
            subprocess.check_output("blastn -db %s -query temp.fasta -out %s -outfmt 5 -num_threads 2 -evalue 1e-20"%(sys.argv[3],out_file),shell=True)
