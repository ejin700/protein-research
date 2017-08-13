# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 13:01:54 2016

@author: ewj
"""
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO, SeqIO
import seq_utils
import pandas as pd
import os, itertools
import numpy as np
import scipy.stats as stat
import matplotlib.pyplot as plt


def graph_data(data,out_file):
    data.sort()
    mean = np.mean(data)
    std = np.std(data)    
    
    #axes = plt.gca()
    #axes.set_xlim([0,7])
    
    pdf = stat.norm.pdf(data,mean,std)
    plt.plot(data,pdf)
    plt.text(0,0, 'Mean = %s \nStandard Deviation = %s' %(mean,std))
    plt.savefig(out_file)
    #plt.close()


def sort_weird(path,out_file):
    bad_files = os.path.join(path,'Weird')
    if not os.path.exists(bad_files):
        os.makedirs(bad_files)
     
    weird = []
    bad = []
    files = [file for file in os.listdir(path) if 'Weird' not in file]
    for file in files:  
        name = file[8:-6] #for aligned_NAME.fasta
        in_file = os.path.join(path,file)
        
        #check for experimental sites that don't have standard consesus motif
        aln = AlignIO.read(in_file,'fasta')
        for seq in aln:
            s = str(seq.seq)
            pat,pos = seq_utils.find_seq_pattern(s,pattern='[n][-]*?[a-oq-z][-]*?[a-bd-ru-z]',ret_start=False)
            for p in pos:
                if len(p) != 0:
                    weird.append((name,seq.id.split('|')[1],p[0],s[p[0]:p[1]]))
                    bad.append(file)
                    
            pat1,pos1 = seq_utils.find_seq_pattern(s,pattern='[n][-]*?[p][-]*?[stc]',ret_start=False)
            for p1 in pos1:
                if len(p1) != 0:
                    weird.append((name,seq.id.split('|')[1],p1[0],s[p1[0]:p1[1]]))
                    bad.append(file)
    
    #return csv with info on weird proven glycosylation sites
    df = pd.DataFrame(weird,columns=['gene','id','position','sequence'])
    df.to_csv(out_file)
    
    b = set(bad)
    for f in b:
        os.rename(os.path.join(path,f),os.path.join(bad_files,f))
            

def conserved_site(in_path,out_file):      
    data = []    
    files = [file for file in os.listdir(in_path) if 'Bad' not in file]
    for file in files:
        name = file[8:-6] #for aligned_NAME.fasta
        in_file = os.path.join(in_path,file)
        
        aln = AlignIO.read(in_file,'fasta')
        seq_count = 0
        exp = []
        
        #find all positions with experimentally proven sites
        for seq in aln:
            seq_count += 1
            s = str(seq.seq)
            pat,pos = seq_utils.find_seq_pattern(s,pattern='[n][-]*?[a-oq-z][-]*?[stc]',ret_start=False)
            for p in pos:
                if len(p) != 0:
                    exp.append(p)
        
        ref = sorted(set(exp))
        #for each of those positions, count how many sequences have glycosylation sites there
        for site in ref:
            count = 0
            for seq in aln:
                s = str(seq.seq)
                pat1,pos1 = seq_utils.find_seq_pattern(s,pattern='[nN][-]*?[a-oA-Oq-zQ-Z][-]*?[stcSTC]',ret_start=False)
                for p in pos1:
                    #allow for plus or minus 2 position shift                    
                    if abs(p[0] - site[0]) <= 2:
                        count +=1
            data.append((name,site[0],count/seq_count))
    
    df = pd.DataFrame(data,columns=['gene','position','percent'])
    df = df.drop_duplicates()
    df.to_csv(out_file)
    

def comp_entire_seq(in_path,evidence_file,out_file):    
    msk = pd.read_csv(evidence_file)
    msk = msk[msk.gene_name != '_none_']
    evidence = set([x.upper() for x in msk.id])
    
    data = []
    files = [file for file in os.listdir(in_path) if 'Bad' not in file]
    for file in files:
        name = file[8:-6] #for aligned_NAME.fasta
        file_in = os.path.join(in_path,file)
        entries = []
        good = []
        
        for seq in SeqIO.parse(file_in,'fasta'):
            entries.append(seq.id.split('|')[1])
         
        #store entries with experimental data 
        for x in entries:
            if x in evidence:
                good.append(x)
                
        #do clustalo alignment on fasta file, output percent difference matrix
        in_file, out_file, matrix = file_in, 'out.fasta','matrix'
        clustalomega_cline = ClustalOmegaCommandline(infile=in_file,outfile=out_file,distmat_out=matrix,force=True,distmat_full=True)
        clustalomega_cline()
        
        with open('matrix','r') as m:
            num = int(m.readline())
            l = m.read().split('\n')
            
            #iterate through entries with experimental data
            for entry in good:
                for row in l:
                    if entry in row:
                        #find mean percent difference relative to selected entry
                        dist = [float(z) for z in row.split(' ') if '0.' in z]
                        avg = 1 - sum(dist)/num
                        data.append((name,entry,avg))
    
    df = pd.DataFrame(data, columns=['gene','entry','percent'])
    df = df.drop_duplicates()
    #print(df)
    df.to_csv('total_conserved.csv') #wont work with out_file?!?

                
def site_to_seq(tot_file,sites_file,out_file): #code varies based on intended analysis
    totals = pd.read_csv(tot_file)
    unique = set(totals['gene'])
    ref = {}
    #finds mean total sequence conservation by averaging %conservations of sequences with exp data
    for gene in unique:
        msk = totals[totals['gene'] == gene]
        avg = np.mean(msk['percent'])
        ref[gene] = avg
    
    data = []
    
    sites = pd.read_csv(sites_file)
    for gene in unique:
        msk1 = sites[sites['gene'] == gene]
        tot_per = ref.get(gene)
        
        #commented out code is to find the target list (%potential position conserved/%total conserved)
        #for x,y,z in itertools.zip_longest(msk1['position'],msk1['percent'],msk1['type']):
         #   data.append((gene,x,y/tot_per,z))
        
        #iterates through to find %proven position conserved/%total conserved
        for x,y in itertools.zip_longest(msk1['position'],msk1['percent']):
            data.append((gene,x,y/tot_per))
            
    
    #df = pd.DataFrame(data,columns=['gene','position','percent','type'])
    df = pd.DataFrame(data,columns=['gene','position','percent'])
    df.to_csv(out_file)


def flanking(n,direction,in_path):
    
    data = []    
    files = [file for file in os.listdir(in_path) if 'Bad' not in file]
    for file in files:
        name = file[8:-6] #for aligned_NAME.fasta
        in_file = os.path.join(in_path,file)
        
        aln = AlignIO.read(in_file,'fasta')        
        left = []
        right = []
        seq_count = 0
        #find position of proven sites
        for seq in aln:
            seq_count += 1
            s = str(seq.seq)
            pat,pos = seq_utils.find_seq_pattern(s,pattern='[n][-]*?[a-oq-z][-]*?[stc]',ret_start=False)
        
            for p in pos:
                if len(p) != 0:
                    #find and store bases to the left and right of proven sites
                    z = seq_utils.extract_info_before_after(s,p,n)
                    left.append((p,z[0][0]))
                    right.append((p,z[0][2]))

        if direction == 'left':
            for tup in left:
                count = 0
                for seq in aln:
                    s = str(seq.seq)
                    #check to see if flanking base is same as that for exp proven site
                    if s[tup[0][0]-n:tup[0][0]] == tup[1]:
                        count += 1
                data.append((name,tup[0][0],count/seq_count))
            
        if direction == 'right':
            for tup in right:
                count = 0
                for seq in aln:
                    s = str(seq.seq)
                    if s[tup[0][1]:tup[0][1]+n] == tup[1]:
                        count += 1
                data.append((name,tup[0][0],count/seq_count))
            
    data = sorted(data)    
    df = pd.DataFrame(data,columns=['gene','position','percent'])
    #take mean in case there are two exp sites with different flanking bases
    clean = pd.DataFrame(df.groupby(['gene','position'])['percent'].mean())
    clean.to_csv('%s_%d.csv' %(direction,n))


def flanking_to_sites(site_file,flank_file,out_file):      
    data = [] 
    #read in info for site conservation and flanking base conservation
    df_flank = pd.read_csv(flank_file)
    df_site = pd.read_csv(site_file)
    unique = list(set(df_flank['gene']))
    
    for gene in unique:
        flank = df_flank[df_flank['gene'] == gene]
        site = df_site[df_site['gene'] == gene]
        
        #iterate through and divide %flanking conserved/%site conserved
        for f,s in itertools.zip_longest(flank['percent'],site['percent']):
            data.append(f/s)

    #graph_data(data,out_file)


def potential_sites(in_path,out_file): #can probably be combined with conserved_site 
    data = []    
    files = [file for file in os.listdir(in_path) if 'Bad' not in file]
    for file in files:
        name = file[8:-6] #for aligned_NAME.fasta
        in_file = os.path.join(in_path,file)
        
        aln = AlignIO.read(in_file,'fasta')
        seq_count = 0
        a = []
        
        #find all potential glycosylation sites
        for seq in aln:
            seq_count += 1
            s = str(seq.seq)
            pat,pos = seq_utils.find_seq_pattern(s,pattern='[N][-]*?[A-OQ-Z][-]*?[STC]',ret_start=False)
            for p in pos:
                if len(p) != 0:
                    a.append(p)
        #file all proven glycosylatin sites
        exp = []
        for seq in aln:
            s = str(seq.seq)
            pat,pos = seq_utils.find_seq_pattern(s,pattern='[n][-]*?[a-oq-z][-]*?[stc]',ret_start=False)
            for p in pos:
                if len(p) != 0:
                    exp.append(p)            
        
        #remove proven sites from potential sites
        tmp = [x for x in a if x not in exp]
        ref = sorted(set(tmp))
        for site in ref:
            count = 0
            third = []
            for seq in aln:
                s = str(seq.seq)
                pat1,pos1 = seq_utils.find_seq_pattern(s,pattern='[N][-]*?[A-OQ-Z][-]*?[STC]',ret_start=False)
                for p in pos1:
                    if abs(p[0] - site[0]) <= 2:
                        #specify what type of glycosylation site (S/T/C)
                        third.append(s[p[1]-1])
                        count +=1
            letters = set(third)
            t = "".join(letters)
            data.append((name,site[0],count/seq_count,t))
    
    df = pd.DataFrame(data,columns=['gene','position','percent','type'])
    df = df.drop_duplicates()
    df.to_csv(out_file)    
     

        
        
        
        
        
        
        

