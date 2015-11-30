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
        if(start==0):
            inversion_bp_1 = mzebra_ref[scaffold][int(start):int(start)+1000].seq
        elif(stop == len(inversion)):
            inversion_bp_2 = mzebra_ref[scaffold][int(stop)-1000:int(stop)].seq
        else:
            inversion_bp_1 = mzebra_ref[scaffold][int(start):int(start)+1000].seq
            inversion_bp_2 = mzebra_ref[scaffold][int(stop)-1000:int(stop)].seq
        with open('temp.fasta','w') as TEMP, open('scaffold.fasta','w') as SCAFFOLD:
            #TEMP.write(">%s_%s_%s\n"%(scaffold,start,stop))
            #TEMP.write(inversion)
            TEMP.write(">bp1\n")
            TEMP.write(inversion_bp_1+"\n")
            TEMP.write(">bp2\n")
            TEMP.write(inversion_bp_2)
            out_file = scaffold+"_"+start+"_"+stop+".xml"
            SCAFFOLD.write(">%s\n"%(scaffold))
            SCAFFOLD.write("%s"%(mzebra_ref[scaffold][0:].seq))
            subprocess.check_output("makeblastdb -in scaffold.fasta -input_type fasta -dbtype nucl -title m_zebra -out temp_mz_v0_db",shell=True)
            subprocess.check_output("blastn -db %s -query temp.fasta -out %s -outfmt 5 -num_threads 2 -evalue 1e-20"%(sys.argv[3],out_file),shell=True)
            #subprocess.check_output("rm temp_mz*",shell=True)
            subprocess.check_output("rm scaffold.fasta",shell=True)
