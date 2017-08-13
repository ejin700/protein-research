# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 15:20:53 2016

@author: ewj
"""

import sequence_alignment as s
import os


def main():
    
    #specify path
    path = r'C:\Users\ewj\Documents\Python Scripts'
    table_path = os.path.join(path,'geneTables')
    fasta_path = os.path.join(path,'sites_lower')
    unique_path = os.path.join(path,'toAlign')
    
    s.remove_dup_org('n_linked.csv',table_path,fasta_path,unique_path)
    
    aligned_path = os.path.join(path,'aligned')
    s.align(unique_path,aligned_path)

    pretty_path = os.path.join(path,'pretty')
    s.pretty(aligned_path,pretty_path)
    
    #remove intermediate files
    os.remove('out.fasta')
    os.remove('temp.fasta')
    os.remove('matrix')
    

main()