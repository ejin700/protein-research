# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:46:12 2016

@author: ewj
"""
from Bio import SeqIO
import myUniprotIO
from UniprotUtils import get_feature_frame
from bioservices.uniprot import UniProt
u = UniProt(verbose = True)
import pandas as pd
from tqdm import tqdm
import os, itertools


def evidence(feature,letter,out_file):
    #read query file
    x = myUniprotIO.UniprotIterator(open('query.xml','r'), return_raw_comments=True)

    #parse for wanted data (i.e. gene name, sequence id, position, etc.)
    L = []
    for rec,seqrec in tqdm(enumerate(x)):
        gene_name = seqrec.annotations.get('gene_name_primary','_none_')
        
        t= (
            get_feature_frame(seqrec,stype=feature,filter_val=letter)
            .assign(rec=rec, id=seqrec.id, gene_name=gene_name)
            .rename(columns={'start':'position'})
            )    
        L.append(t)

    msk = pd.concat(L)
    
    #filter for selected types of experimental evidence
    df = msk.loc[msk['evidence_type'].str.contains('269|305|244|213', na=False)]
    df.to_csv(out_file)  
    
    return df
    
    
def build_gene_query(gene_names):  
    #compile string for uniprot search
    L = ['gene_exact:%s' %x for x in gene_names]
    s = ' OR '.join(L)
    q = '(%s) taxonomy:mammalia fragment:no' %s
    
    return q
    
    
def chunks(l,n):
    for i in range(0, len(l), n):
        yield l[i:i+n]
        
        
def get_file_by_gene_names(gene_names,output_prefix,frmt=None,limit=20,verbose=True):
    #split list of gene names into chuncks, then search by chuncks
    for i,names in enumerate(chunks(gene_names,limit)):
        #specify output file format        
        if frmt is not None:
            assert frmt in ['tab','xml','fasta']
            ofile = '{0}_{1}_.{2}'.format(output_prefix,i,frmt)
        else:
            ofile = None
            
        query = build_gene_query(names)
        if verbose:
            print(ofile)
            print(query)
        
        #output search results in selected format
        if ofile is not None:
            with open(ofile,'w') as f:
                f.write(u.search(query,frmt=frmt))
                
                
def indiv_files(gene_names,df,gene_path,fasta_path,fasta_file):
   
    #create output folders for gene tables and fasta files
    if not os.path.exists(gene_path):
        os.makedirs(gene_path)    
    if not os.path.exists(fasta_path):
        os.makedirs(fasta_path)
    
    #iterate through list of gene name results in dataframe    
    for name in gene_names:
        options = df['Gene names'].str.upper()
        good = []
        #keep entry if first word in gene name result is same as searched gene name
        for n in options:
            l = n.split()
            if l[0] == name:
                good.append(n)
        #g = df.loc[df['Gene names'].str.upper() == name] #exact name match
        g = df.loc[df['Gene names'].isin(good)]
        #output gene table of wanted results
        g.to_csv(gene_path + '\%s.csv' %name)
        
        #parse large all_fasta file and pick out ones selected in gene table        
        input_iterator = SeqIO.parse(open(fasta_file,'r'), 'fasta')
        Entries = g.Entry.values
        filter_iterator = (x for x in input_iterator if x.id.split('|')[1] in Entries)
        
        with open(fasta_path + '\%s.fasta' %name,'w') as f:
            SeqIO.write(filter_iterator, f, 'fasta')
    

def rm_fasta_dup(in_file,out_file):
    ids = []
    desc = []
    sequences = []
    
    #compile single large fasta file of unique fasta sequences
    for seq in SeqIO.parse(in_file,'fasta'):
        #only add sequence to large file if it isn't already there
        if seq.id not in ids:
            ids.append(seq.id)
            desc.append(seq.description)
            sequences.append(str(seq.seq))
    #output single fasta file with all sequences
    with open(out_file,'w') as f:
        for x in range(len(ids)):
            f.write('>' + desc[x] + '\n' + sequences[x] + '\n')
          

def switch_letters(letters,pos,num):
    for x in range(num):
        letters[pos+x] = letters[pos+x].lower()

    
def lower(unique_ids,pos_list,seq_list,width=60):

    old_seq = []
    new_seq = []
    
    for acc,pos in itertools.zip_longest(unique_ids,pos_list):  
        for seq in seq_list:
            if acc in seq:
                #isolate wanted sequence from fasta file
                full_seq = '>' + seq
                ind = full_seq.find('\n')
                header, body = full_seq[:ind], full_seq[ind+1:]
                old_seq.append(full_seq)
                letters = list(body)
                
                #switch exp sites to lowercase (depending on location in file)
                for p in pos:
                    if p <= width - 3:
                        switch_letters(letters,p,3)
                    elif (width-3) < p < width:
                        switch_letters(letters,p,4)                        
                    elif p >= width:
                        loc = p + p//width
                        if (p % width) > (width-3):
                            switch_letters(letters,loc,4)
                        else:
                            switch_letters(letters,loc,3)
                
                #construct new string with exp data lowercase
                string = header + '\n' + ''.join(letters)
                new_seq.append(string)        
        
    return old_seq, new_seq


def mark_exp(unique_genes,evidence_file,out_path,fasta_path):
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)     
    
    df = pd.read_csv(evidence_file)   
    
    for gene in unique_genes:
        #isolate entries that have experimental data associated with it
        wanted = df[df['gene_name'].str.upper() == gene]
        unique_ids = wanted.id.unique()
        l = []
        #compile list of positions that are experimentally proven
        for x in unique_ids:
            pos_list = wanted[wanted['id'] == x]
            positions = pos_list['position'].tolist()
            index = [int(pos) for pos in positions]
            l.append(index)
        
        #retrieve necessary fasta file from fasta path
        ref_file = fasta_path + '\%s.fasta' %gene
        with open(ref_file, 'r') as f:
            s = f.read()
            z = s.split('>')
        
        #assign original sequence and new, marked lowercase sequence
        old_seq, new_seq = lower(unique_ids,l,z,width=60)
        
        #replace old sequene with new, then write out new file
        for old,new in itertools.zip_longest(old_seq,new_seq):
            s = s.replace(old,'temp').replace(new,old).replace('temp',new)
            
        with open(out_path + '\lower_%s.fasta' %gene,'w') as n:
            n.write(s)    
    
    
def sort_bad(path):
    bad_files = os.path.join(path,'Bad')
    if not os.path.exists(bad_files):
        os.makedirs(bad_files)
        
    #remove files with less than 5 sequences in them, move to Bad folder
    files = [file for file in os.listdir(path) if 'Bad' not in file]
    for file in files:  
        p = os.path.join(path,file)
        with open(p,'r') as f:
            s = f.read()
            count = s.count('>')
        if count < 5:
            os.rename(os.path.join(path,file),os.path.join(bad_files,file))    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
