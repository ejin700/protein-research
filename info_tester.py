# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:46:35 2016

@author: ewj
"""

import get_info as i
from bioservices.uniprot import UniProt
u = UniProt(verbose=True)
import pandas as pd    
import subprocess, os, glob


def main():  
    
    #specify path. VERY IMPORTANT!!
    path = r'C:\Users\ewj\Documents\Python Scripts'    
    
    #perform initial query on Uniprot
    query = 'annotation:(type:carbohyd "n linked" evidence:ECO_0000269) taxonomy:mammalia fragment:no'
    with open('query.xml','w') as q:
        q.write(u.search(query,frmt='xml'))
    
    #parse query file for wanted data and get list of unique genes
    df = i.evidence('glycosylation','N','n_linked.csv')  
    unique_genes = list(set([x.upper() for x in df.gene_name]))
    
    #mass search uniprot website then compile all tab files into one
    i.get_file_by_gene_names(unique_genes,'TAB',limit=100,frmt='tab')
    allTab = glob.glob(os.path.join(path, '*tab'))
    df_all = pd.concat((pd.read_table(f, sep = '\t', header=0) for f in allTab))
    df_all['Gene names'] = df_all['Gene names'].str.upper()    
    df_all = df_all.drop_duplicates()
    
    #mass search uniprot website then compile all fasta files into one
    i.get_file_by_gene_names(unique_genes,'FASTA',limit=100,frmt='fasta')
    subprocess.call("copy *.fasta allFasta.fasta", shell = True) #does deep copy!
    i.rm_fasta_dup('allFasta.fasta','fastaRef.fasta')

    #parse large files for idividual gene name results
    gene_path = os.path.join(path,'geneTables')
    fasta_path = os.path.join(path,'fastaFiles')
    i.indiv_files(unique_genes,df_all,gene_path,fasta_path,'fastaRef.fasta')
    
    #convert experimental data sites to lowercase
    lower_path = os.path.join(path, 'sites_lower')
    i.mark_exp(unique_genes,'n_linked.csv',lower_path,fasta_path)
    i.sort_bad(lower_path)
    
    #clean up intermediate files
    for f in glob.glob('*.fasta'):
        os.remove(f)
    for t in glob.glob('*.tab'):
        os.remove(t)
main()

