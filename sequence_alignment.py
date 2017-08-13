# -*- coding: utf-8 -*-
"""
Created on Mon Jun  13 11:25:35 2016

@author: ewj
"""

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import seq_utils
import pandas as pd
import os, collections


def dist_mat(entries,file,out_list):
    #parse fasta file for wanted sequences
    input_iterator = SeqIO.parse(open(file,'r'),'fasta')
    filter_iterator = (x for x in input_iterator if x.id.split('|')[1] in entries)
    with open('temp.fasta','w') as f:
        SeqIO.write(filter_iterator,f,'fasta')        
    
    #run distance matrix/alignment on select sequences with clustalo    
    in_file, out_file, matrix = 'temp.fasta','out.fasta','matrix'
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file,outfile=out_file,distmat_out=matrix,force=True,distmat_full=True)
    clustalomega_cline()    
    
    #read in distance matrix     
    with open('matrix','r') as m:
        num = int(m.readline())
        l = m.read().split()
        dist = [float(x) for x in l if '0.' in x]
        rows = [dist[i:i+num] for i in range(0,len(dist),num)]
        means = []
        #find mean difference for every sequence
        for row in rows:
            means.append(sum(row)/num)
        #return sequence with lowest mean difference
        want_ind = means.index(min(means))
        out_list.append(entries[want_ind])


def two_left(df,out_list):
    #filter by length
    max_len = df['Length'].max()
    x = df.loc[df['Length'] == max_len]
    #if both sequences same length, randomly choose one
    if len(x) == 2:
        out_list.append(df.Entry.values[0])
    else:
        out_list.append(x.Entry.values)
        
        
def three_opt(df,file,out_list):
    #filter for reviewed sequence
    reviewed = df.loc[df['Status'] == 'reviewed']
    n = len(reviewed)
    if n == 1:
        out_list.append(reviewed.Entry.values)
    elif n == 2:
            two_left(reviewed,out_list)
    elif n == 0:
        dist_mat(df.Entry.values,file,out_list)
    else:
        dist_mat(reviewed.Entry.values,file,out_list)
            

def sort_bad(path):
    bad_files = os.path.join(path,'Bad')
    if not os.path.exists(bad_files):
        os.makedirs(bad_files)
        
    files = [file for file in os.listdir(path) if 'Bad' not in file]
    for file in files:  
        p = os.path.join(path,file)
        #remove flies that contain less than 5 sequences
        with open(p,'r') as f:
            s = f.read()
            count = s.count('>')
        if count < 5:
            os.rename(os.path.join(path,file),os.path.join(bad_files,file))
                

def remove_dup_org(evidence_file,table_path,fasta_path,out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    #read in evidence file, record which sequences have exp data    
    evidence = pd.read_csv(evidence_file)
    e = list(set([x.upper() for x in evidence.id]))
    files = [file for file in os.listdir(fasta_path) if 'Bad' not in file]

    for fasta in files:
        file_path = os.path.join(fasta_path,fasta)
        name = fasta[6:-6] #for lower_NAME.fasta
        df = pd.read_csv(table_path + '\%s.csv' %name)
        #get list of duplicate organisms
        dup = [item for item,count in collections.Counter(df['Organism']).items() if count > 1]        
        unique = [x for x in df['Organism'] if x not in dup]
        keep = []
        
        #keep all sequences for unique organisms
        for u in unique:
            msk = df.loc[df['Organism'] == u]
            keep.append(msk.Entry.values)
        
        #find 'most representative' sequence by evidence, reviewed, min dist, length
        for org in dup:
            msk = df.loc[df['Organism'] == org]
            entries = msk.Entry.values.tolist()

            if len(entries) == 2: #start with 2
                good = [x for x in entries if x in e]
                exp = msk.loc[msk['Entry'].isin(good)]
                if len(exp) == 1:
                    keep.append(exp.Entry.values)
                else: #none or 2
                    reviewed = msk.loc[msk['Status'] == 'reviewed']
                    if len(reviewed) == 1: 
                        keep.append(reviewed.Entry.values)
                    else:
                        two_left(msk,keep)
    
            else: #start with more than 2
                good = [x for x in entries if x in e]
                exp = msk.loc[msk['Entry'].isin(good)]
                if len(exp) == 1:
                    keep.append(exp.Entry.values)
                    
                elif len(exp) == 2:
                    reviewed = exp.loc[exp['Status'] == 'reviewed']
                    if len(reviewed) == 1:
                        keep.append(reviewed.Entry.values)
                    else:
                        two_left(exp,keep)
                
                elif len(exp) == 0: #none have exp data-check reviewed
                    three_opt(msk,file_path,keep)
        
                else: #multiple have exp data
                    three_opt(exp,file_path,keep)
        

        input_iterator = SeqIO.parse(open(file_path,'r'),'fasta')
        filter_iterator = (x for x in input_iterator if x.id.split('|')[1] in keep)
        
        out = os.path.join(out_path,'unique_%s' %fasta)
        with open(out,'w') as f:
            SeqIO.write(filter_iterator,f,'fasta')
    
    sort_bad(out_path)   


def align(in_path,out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        
    files = [file for file in os.listdir(in_path) if 'Bad' not in file]    
    for file in files:
        name = file[13:] #for unique_lower_NAME.fasta
        in_file = os.path.join(in_path,file)
        out_file = os.path.join(out_path,'aligned_%s' %name)
        #run clustal omega alignment on sequences
        clustalomega_cline = ClustalOmegaCommandline(infile=in_file,outfile=out_file,verbose=True,auto=True,force=True)
        clustalomega_cline()
        

def my_align(aln,width=60,pad=' '*2):
    
    myseq = []
    pos_all = []
    
    #find all potential and proven glycosylation sites
    for seq in aln:
        s = str(seq.seq)
        pat,pos = seq_utils.find_seq_pattern(s,pattern='[nN][-]*?[a-oA-Oq-zQ-Z][-]*?[stcSTC]',ret_start=False)
        for p in pos:        
            pos_all.append(p)
    
    unique_positions = set(pos_all)
     
    for seq in aln:
        s = str(seq.seq)
        color_list = ['yellow','cyan','green','red','magenta']
        pos_list = []
        mutations = []
        s_mut = []
        
        #check for mutated sites
        for mut in unique_positions:
            if s[mut[0]] == 'D' and s[mut[0] + 1] != 'P' and s[mut[0] + 2] in 'STC':
                mutations.append(mut)
            elif s[mut[0]] == 'S' and s[mut[0] + 1] != 'P' and s[mut[0] + 2] in 'STC':
                s_mut.append(mut)
        
        #return list of positions for NXS/T/C sites
        for x in ['Ss','Tt','Cc']:
            pat1,pos1 = seq_utils.find_seq_pattern(s,pattern='[nN][-]*?[a-oA-Oq-zQ-Z][-]*?[%s]' %x,ret_start=False)
            pos_list.append(pos1)

        pos_list.append(mutations)
        pos_list.append(s_mut)
        
        #highlight sequences based on color and position
        x = seq_utils.seq_to_pretty_list(s,pos_list,color_list)
        myseq.append(x)        
    
    #write out new sequence in new formated string
    s = ''
    for j in range(len(aln[0])//width+1):
        
        LB = j * width
        UB = LB + width
        
        for i,seq in enumerate(aln):
            q = ''.join(myseq[i][LB:UB])
            
            s += '%30s %s %s %s %i\n'%(seq.name,pad,q,pad,UB)
        s+='\n'

    return s
    

def pretty(in_path,out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    
    #write out pretty txt file for each fasta file    
    files = [file for file in os.listdir(in_path) if 'Bad' not in file]
    for file in files:
        name = file[8:-6] #for aligned_NAME.fasta
        in_file = os.path.join(in_path,file)
        out_file = os.path.join(out_path,'%s.txt' %name)       
        alignment = AlignIO.read(in_file,'fasta')
        pretty = my_align(alignment)
        with open(out_file,'w') as f:
            f.write(pretty)
        

