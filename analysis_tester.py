# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 14:56:38 2016

@author: ewj
"""

import stat_analysis as s
import os


def main():

    path = r'C:\Users\ewj\Documents\Python Scripts'
    
    aligned_path = os.path.join(path,'aligned')
    s.sort_weird(aligned_path,'weird.csv') 

    s.conserved_site(aligned_path,'sites_conserved.csv')
    s.potential_sites(aligned_path,'potential_sites_conserved.csv')
    s.comp_entire_seq(aligned_path,'n_linked.csv','total_conserved.csv')
  
    s.site_to_seq('total_conserved.csv','sites_consereved.csv','exp_site2seq.csv')  
    #s.site_to_seq('total_conserved.csv','potential_sites_conserved.csv','target_list.csv')
    
    #for n in range(1,11):
     #   s.flanking(n,'left',aligned_path)
    
    #for n in range(1,6):
     #   s.flanking_to_sites('conserved_sites.csv','right_%s.csv' %n,'right_%s.png' %n)
    #s.flanking_to_sites('conserved_sites.csv','left_1.csv','left_1.png')

main()