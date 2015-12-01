#!/usr/bin/env python
import matplotlib.pyplot as plt 
import sys
import numpy as np
import glob
from matplotlib.widgets import Slider, Button, RadioButtons
from collections import defaultdict,Counter
import os
import re
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans,k_means
from sklearn.ensemble import ExtraTreesClassifier

def create_bed(log_file,output_file):
    """ Creates a bed file from the logfile.txt file generated earlier
    """

    with open(log_file,'r') as LOG,open(output_file,'w') as OUT:
        for line in LOG:
            line = line.rstrip('\n')
            scaffold = line.split('\t')[0].split(':')[0]
            start = line.split('\t')[0].split(':')[1].split('-')[0]
            start = start.replace(',','')
            stop = line.split('\t')[0].split(':')[1].split('-')[1]
            stop = stop.replace(',','')
            OUT.write(scaffold+"\t"+start+"\t"+stop+"\n")

def cluster_strains(folder):
    
    #abs_alt_inv_set = defaultdict(dict)
    #abs_ref_inv_set = defaultdict(dict)
    #perc_alt_inv_set = defaultdict(dict)
    #perc_ref_inv_set = defaultdict(dict)

    [perc_alt_inv_set,abs_alt_inv_set] = parse_geno_file(folder,True)

    #print prec_alt_inv_set['MC'].keys()

    num_features = len(perc_alt_inv_set['MC'].keys())

    print (num_features)

    feature_matrix = np.zeros((13,num_features))
    strain_labels = []
    inversion_features = []

    j = 0 
    for strain in perc_alt_inv_set.keys():
        strain_labels.append(strain)

    for inversion in perc_alt_inv_set['MC'].keys():
        i = 0 
        for strain in strain_labels:
            feature_matrix[i][j] = perc_alt_inv_set[strain][inversion]
            inversion_features.append(inversion)
            i = i + 1 
        j = j + 1

    ## Clustering Stuff
    estimator = KMeans(init='k-means++', n_clusters=2, n_init=10)
    estimator.fit(feature_matrix)
    print (estimator.labels_)
    print (strain_labels)
    #centroids,labels,inertia = k_means(feature_matrix,2)
    
    ## Random Forests for extracting important features
    forest = ExtraTreesClassifier(n_estimators=250,random_state=0)
    forest.fit(feature_matrix, estimator.labels_)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],axis=0)
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    #print("Feature ranking:")
    
    #for f in range(feature_matrix.shape[1]):
        #print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
    
    with open('log_file.txt','a') as LOG:
        for i in range(0,10):
            LOG.write(inversion_features[indices[i]]+"\n")

def output_inversions(folder):
    """ Output to a file all the inversions above a certain threshold 
    folder : The input folder
    """
    
    [perc_alt_alleles,abs_alt_alleles] = parse_geno_file(folder,True) ## Call the parser, the returned object is a dictionary of dictionaries

    cluster_strains(folder)

    strains = [' ','MC','CL','CM','CN','TI','DC','MS','CV','PN','AC','LF','MP','MZ']
    strain_types = ['Sand','Sand','Sand','Sand','Sand','Sand','Sand','Sand','Rock','Rock','Rock','Rock','Rock']

    with open('out_file1.txt','w') as OUT1, open('out_file2.txt','w') as OUT2:
        OUT1.write("\t".join(strains)+"\n")
        OUT2.write("\t".join(strains)+"\n")
        with open('log_file.txt','r') as LOG:
            for line in LOG:
                line = line.strip('\n')
                log_inversion = line.split('\t')[0]
                output1 = []
                output2 = []
                output1.append(line)
                output2.append(line)
                for strain in strains:
                    if strain != " ":
                        output1.append(str(perc_alt_alleles[strain][log_inversion]))
                        output2.append(str(abs_alt_alleles[strain][log_inversion]))
                OUT1.write("\t".join(output1)+"\n")
                OUT2.write("\t".join(output2)+"\n")

    #subprocess.check_output('cat out_file1.txt > out_file.txt',shell=True)
    #subprocess.check_output('cat out_file2.txt >> out_

def parse_geno_file(folder,return_flag):
    """ Parses the outputs from svviz to get the allele counts
        folder : input folder to search for the output files
        return_inversions : a boolean to specify the return object, either a dictionary of lists in case of the plot function or a dict of
                            dicts when called by the output_inversions function
    """

    perc_alt = defaultdict(list)
    perc_ref = defaultdict(list)
    abs_alt = defaultdict(list)
    abs_ref = defaultdict(list)

    perc_alt_inv = defaultdict(dict)
    perc_ref_inv = defaultdict(dict)
    abs_alt_inv = defaultdict(dict)
    abs_ref_inv = defaultdict(dict)

    for geno_file in glob.glob(folder+'*_test_summary.tsv'):
        strain = geno_file.split('/')[-1].split('_')[0]
        #print strain
        prev_coordinate = "0"
        count = 0
        alt_allele = {}
        amb_allele = {}
        ref_allele = {}
        flag = 0 

        TEMP_HANDLE = open(geno_file,'r')
        for line in TEMP_HANDLE:
            line = line.rstrip('\n')

            if(line[0]!='v'): ## Skip the header
                coordinate = line.split('\t')[0].split('::')[-1]
                if(coordinate != prev_coordinate):
                    #prev_coordinate = coordinate
                    count = count + 1
                    if(count == 1):
                        if(line.split('\t')[-3]!='alt'): ## No reads supporting the alternate allele
                            flag = 1 
                            alt_allele[coordinate] = 0
                            amb_allele[coordinate] = int(line.split('\t')[-1])
                            #print line
                        else:
                            alt_allele[coordinate] = int(line.split('\t')[-1])
                    if(count == 2):
                                amb_allele[coordinate] = int(line.split('\t')[-1])
                    if(count == 3):
                        if(line.split('\t')[-3]!='ref'): ## No reads supporting the reference allele (all are ambiguous)
                            ref_allele[coordinate] = 0
                        else:
                            ref_allele[coordinate] = int(line.split('\t')[-1])
                        prev_coordinate = coordinate
                        count = 0
                    if(flag == 1): ## The case where there are no alternate allele reads, counter is incremented to account for changed numbering
                        count = count + 1 
                        flag = 0 

        
        for key in alt_allele:
            if(alt_allele[key]+ref_allele[key]!= 0): ## Check to see if the denominator is not zero
                abs_alt[strain].append(float(alt_allele[key]))
                abs_ref[strain].append(float(ref_allele[key]))
                perc_alt[strain].append(float(alt_allele[key])/(alt_allele[key]+ref_allele[key]))
                perc_ref[strain].append(float(ref_allele[key])/(alt_allele[key]+ref_allele[key]))


                abs_alt_inv[strain][key] = float(alt_allele[key])
                abs_ref_inv[strain][key] = float(ref_allele[key])
                perc_alt_inv[strain][key] = float(alt_allele[key])/(alt_allele[key]+ref_allele[key])
                perc_ref_inv[strain][key] = float(ref_allele[key])/(alt_allele[key]+ref_allele[key])
     
                

    ## Keep only the common inversions, i.e. those between MC and the rest 
    all_inversions = []
    common_inversions = []
    abs_alt_set = defaultdict(list)
    perc_alt_set = defaultdict(list)

    abs_alt_inv_set = defaultdict(dict)
    perc_alt_inv_set = defaultdict(dict)
    abs_ref_inv_set = defaultdict(dict)
    perc_ref_inv_set = defaultdict(dict)

    Sand = ['CN','CL','CM','CN','TI','CV','MC','MS']
    Rock = ['MZ','AC','PN','LF','MP']


    sand_inversions = []
    rock_inversions = []

    for strain in abs_alt_inv.keys():
        for inversion in abs_alt_inv[strain].keys():
            if(strain in Rock):
                rock_inversions.append(inversion)
            else:
                sand_inversions.append(inversion)
            all_inversions.append(inversion)
    
    
    common_inversions_sand = Counter(sand_inversions)
    common_inversions_rock = Counter(rock_inversions)
    #count_sand = 0
    common_inversions = Counter(all_inversions)
    return_inversions = []
    
    
    #print common_inversions
    for inversion in common_inversions.keys():
        if(common_inversions[inversion]==13):
            return_inversions.append(inversion)
            for strain in abs_alt_inv.keys():
                abs_alt_set[strain].append(abs_alt_inv[strain][inversion])
                perc_alt_set[strain].append(perc_alt_inv[strain][inversion])

                abs_alt_inv_set[strain][inversion] = abs_alt_inv[strain][inversion]
                perc_alt_inv_set[strain][inversion] = perc_alt_inv[strain][inversion]
                abs_ref_inv_set[strain][inversion] = abs_ref_inv[strain][inversion]
                perc_ref_inv_set[strain][inversion] = perc_ref_inv[strain][inversion]

    with open('log_file.txt','w') as LOG_FILE:
        for inversion in abs_alt_inv_set['MC']:
            alternate_allele_sum_rock = 0
            reference_allele_sum_rock = 0
            alternate_allele_sum_sand = 0
            reference_allele_sum_sand = 0  
            for strain in Rock:
                alternate_allele_sum_rock = alternate_allele_sum_rock + abs_alt_inv_set[strain][inversion]
                reference_allele_sum_rock = reference_allele_sum_rock + abs_ref_inv_set[strain][inversion]

            for strain in Sand:
                alternate_allele_sum_sand = alternate_allele_sum_sand + abs_alt_inv_set[strain][inversion]
                reference_allele_sum_sand = reference_allele_sum_sand + abs_ref_inv_set[strain][inversion]

            abs_alt_set['Rock'].append(alternate_allele_sum_rock)
            allele_freq_rock = float((alternate_allele_sum_rock)/(alternate_allele_sum_rock + reference_allele_sum_rock))
            allele_freq_sand = float((alternate_allele_sum_sand)/(alternate_allele_sum_sand + reference_allele_sum_sand))
            perc_alt_set['Rock'].append(allele_freq_rock)
        
            abs_alt_set['Sand'].append(alternate_allele_sum_sand)
            perc_alt_set['Sand'].append(allele_freq_sand)
    
        
            if(allele_freq_rock > float(sys.argv[2]) or allele_freq_sand > float(sys.argv[2])):
                LOG_FILE.write(inversion+"\t"+str(allele_freq_rock)+"\t"+str(allele_freq_sand))
                LOG_FILE.write("\n")
                #print (inversion)
    

    #print "Sand : "+str(count_sand)

    if return_flag == True:
        #print len([abs_alt_inv_set,abs_ref_inv_set,perc_alt_inv_set,perc_ref_inv_set])
        return [perc_alt_inv_set,abs_alt_inv_set]
    else:
        return [abs_alt_set,perc_alt_set]

def plot_mapqs(folder):
    """Plots the Distribution of MAPQ Scores for the Inversion Regions in the respective species
       folder : input folder to look for the files
    """
    
    mapqs_counts = defaultdict(list)
    mapqs_labels = defaultdict(list)

    strain_matcher = re.compile('(.*) \(.*\)')  ## regex to get strain name from radio button click
 
    ## Parse the information from the mapq count files and store the information in a dictionary of lists 
    for mapq_file in glob.glob(folder+"*.count.txt"):
        strain = mapq_file.split("/")[-1].split(".")[0]
        TEMP_HANDLE = open(mapq_file,'r')
        for line in TEMP_HANDLE:
            line = line.strip('\n')
            mapqs_counts[strain].append(int(line.split(' ')[0]))
            mapqs_labels[strain].append(line.split(' ')[1])
    
    fig, ax = plt.subplots()
    indexes = np.arange(len(mapqs_counts[S]))
    rects = plt.bar(indexes,mapqs_counts[S],0.5)
    ax.set_xticks(indexes+1*0.2)
    ax.set_xticklabels(mapqs_labels[S],rotation=30)

    rax1 = plt.axes([0.92, 0.1, 0.08, 0.8])
    radio1 = RadioButtons(rax1, ('MC (Sand)','CL (Sand)','CM (Sand)','CN (Sand)','TI (Sand)','DC (Sand)','MS (Sand)','CV (Sand)','PN (Rock)','AC (Rock)','LF (Rock)','MP (Rock)','MZ (Rock)'), active=0)
    
    def update1(strain):
        match = re.match(strain_matcher,strain)
        strain = match.group(1)
        #indexes = np.arange(len(mapqs_counts[S]))
        for rect,h in zip(rects,mapqs_counts[strain]):
            rect.set_height(h)
        ax.set_xticks(indexes+1*0.2)
        ax.set_xticklabels(mapqs_labels[strain],rotation=30)
        ax.relim()
        ax.autoscale()
        fig.canvas.draw_idle()
        #global S 
        #print S
        #S = strain

    radio1.on_clicked(update1)
    plt.show()


def plot_insertsize():
    """
    Plot the scatter plots for the mean insert size for reference vs alt orientation
    """

def plot_alleles(folder):
    """ Plot reference vs alternate allele frequencies 
        folder : input folder to search for the svviz files
    """

    abs_alt,perc_alt = parse_geno_file(folder,False)

    strain_matcher = re.compile('(.*) \(.*\)')  ## regex to get strain name from radio button click
    strain_matcher2 = re.compile('(.*)-(.*)')
    
    fig, ax = plt.subplots()
    l, = plt.plot(perc_alt['MC'],perc_alt[S],'bs')
    ax.set_ylabel('Alternate Allele Percentage (MC)')
    ax.set_xlabel('Alternate Allele Percentage (MC)')
    
    rax1 = plt.axes([0.92, 0.1, 0.08, 0.8])
    radio1 = RadioButtons(rax1, ('MC (Sand)','CL (Sand)','CM (Sand)','CN (Sand)','TI (Sand)','DC (Sand)','MS (Sand)','CV (Sand)','PN (Rock)','AC (Rock)','LF (Rock)','MP (Rock)','MZ (Rock)','Sand-MC','Rock-MC','Sand-MZ','Rock-MZ','Rock-Sand'), active=0)

    rax2 = plt.axes([0.92, 0.9, 0.08, 0.1])
    radio2 = RadioButtons(rax2, ('%', 'Abs'), active=0)

    #for strain in perc_alt.keys():
        #print strain
        
    #ylabel = plt.axes([1,1,0.08,0.08])

    def update1(strain):
        if strain not in ["Rock-MZ","Sand-MC","Rock-MC","Sand-MZ","Rock-Sand"]:
            match = re.match(strain_matcher,strain)
            strain = match.group(1)
            global base
            base = 'MC'
        else:
            match = re.match(strain_matcher2,strain)
            strain = match.group(1)
            #global base 
            base = match.group(2)

        #print D
        if(D == "%"):
            l.set_xdata(perc_alt[base])
            l.set_ydata(perc_alt[strain])
            ax.set_ylabel('Alternate Allele Percentage (%s)'%(strain))
            ax.set_xlabel('Alternate Allele Percentage (%s)'%(base))
        else:
            l.set_xdata(abs_alt[base])
            l.set_ydata(abs_alt[strain])
            ax.set_ylabel('Alternate Allele Absolute (%s)'%(strain))
            ax.set_xlabel('Alternate Allele Absolute (%s)'%(base))

        ax.relim()
        ax.autoscale()
        fig.canvas.draw_idle()
        global S 
        #print S
        S = strain


    def update2(data):
        if(data == '%'):
            l.set_xdata(perc_alt[base])
            l.set_ydata(perc_alt[S])
            ax.relim()
            ax.autoscale()
            ax.set_ylabel('Alternate Allele Percentage (%s)'%(S))
            ax.set_xlabel('Alternate Allele Percentage (%s)'%(base))
            fig.canvas.draw_idle()
        if(data == 'Abs'):
            l.set_xdata(abs_alt[base])
            l.set_ydata(abs_alt[S])
            ax.relim()
            ax.autoscale()
            ax.set_ylabel('Alternate Allele Absolute (%s)'%(S))
            ax.set_xlabel('Alternate Allele Absolute (%s)'%(base))
            fig.canvas.draw_idle()
            
        global D
        D = data

    radio1.on_clicked(update1)
    radio2.on_clicked(update2)
    plt.show() 
   
""" Dont use this for now
def plot_alleles(folder):
     Plot reference vs alternate allele frequencies 
        folder : input folder to search for the svviz files
    

    abs_alt,perc_alt = collections.deafultdict(dict)

    abs_alt['MC'],perc_alt['MC'] = parse_geno_file(sys.argv[1],False)
    abs_alt['MZ'],perc_alt['MZ'] = parse_geno_file(sys.argv[2],False)
    
    strain_matcher = re.compile('(.*) \(.*\)')  ## regex to get strain name from radio button click
    strain_matcher2 = re.compile('(.*)-(.*)')
    
    fig, ax = plt.subplots()
    l, = plt.plot(perc_alt['MC']['MC'],perc_alt['MC'][S],'bs')
    ax.set_ylabel('Alternate Allele Percentage (MC)')
    ax.set_xlabel('Alternate Allele Percentage (MC)')
    
    rax1 = plt.axes([0.92, 0.1, 0.08, 0.8])
    radio1 = RadioButtons(rax1, ('MC (Sand)','CL (Sand)','CM (Sand)','CN (Sand)','TI (Sand)','DC (Sand)','MS (Sand)','CV (Sand)','PN (Rock)','AC (Rock)','LF (Rock)','MP (Rock)','MZ (Rock)','Sand-MC','Rock-MC','Sand-MZ','Rock-MZ','Rock-Sand'), active=0)

    rax2 = plt.axes([0.92, 0.9, 0.08, 0.1])
    radio2 = RadioButtons(rax2, ('%', 'Abs'), active=0)

    #for strain in perc_alt.keys():
        #print strain
        
    #ylabel = plt.axes([1,1,0.08,0.08])

    def update1(strain):
        if strain not in ["Rock-MZ","Sand-MC","Rock-MC","Sand-MZ","Rock-Sand"]:
            match = re.match(strain_matcher,strain)
            strain = match.group(1)
            global base
            base = match.group(2)
        else:
            match = re.match(strain_matcher2,strain)
            strain = match.group(1)
            #global base 
            base = match.group(2)

        #print D
        if(D == "%"):
            l.set_xdata(perc_alt[base][base])
            l.set_ydata(perc_alt[base][strain])
            ax.set_ylabel('Alternate Allele Percentage (%s)'%(strain))
            ax.set_xlabel('Alternate Allele Percentage (%s)'%(base))
        else:
            l.set_xdata(abs_alt[base][base])
            l.set_ydata(abs_alt[base][strain])
            ax.set_ylabel('Alternate Allele Absolute (%s)'%(strain))
            ax.set_xlabel('Alternate Allele Absolute (%s)'%(base))

        ax.relim()
        ax.autoscale()
        fig.canvas.draw_idle()
        global S 
        #print S
        S = strain


    def update2(data):
        if(data == '%'):
            l.set_xdata(perc_alt[base][base])
            l.set_ydata(perc_alt[base][S])
            ax.relim()
            ax.autoscale()
            ax.set_ylabel('Alternate Allele Percentage (%s)'%(S))
            ax.set_xlabel('Alternate Allele Percentage (%s)'%(base))
            fig.canvas.draw_idle()
        if(data == 'Abs'):
            l.set_xdata(abs_alt[base][base])
            l.set_ydata(abs_alt[base][S])
            ax.relim()
            ax.autoscale()
            ax.set_ylabel('Alternate Allele Absolute (%s)'%(S))
            ax.set_xlabel('Alternate Allele Absolute (%s)'%(base))
            fig.canvas.draw_idle()
            
        global D
        D = data

    radio1.on_clicked(update1)
    radio2.on_clicked(update2)
    plt.show()    
"""

############## MAIN ##################
S = 'MC' ## Global Plot Variables for button initializations
D = "%"
base = 'MC'
if(len(sys.argv) < 3):
    print ("Run the program as : ./create_plot.py <input_folder> <threshold>")
    sys.exit()

folder = sys.argv[1]
threshold = float(sys.argv[2])

plot_alleles(folder)

S = 'MC'
#plot_mapqs(folde
output_inversions(folder)
create_bed('log_file.txt','log_file.bed')
