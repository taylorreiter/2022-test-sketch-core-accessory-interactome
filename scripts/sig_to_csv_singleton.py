#!/usr/bin/env python3
"""
Given a sourmash signature built with --singleton flag, extract the minhash integers, abund 
info, and signature name to a csv file.
"""

from sourmash import signature
from sourmash import sourmash_args
import pandas as pd
import os
import sys
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('signature')       # sourmash signature
    p.add_argument('output')          # output csv file name
    p.add_argument('ksize')           # set kmer size
    args = p.parse_args()

    # load the signature from disk
    sigfp = open(args.signature, 'rt')
    siglist = list(signature.load_signatures(sigfp, ksize = args.ksize))
    
    df_all_sigs = pd.DataFrame(columns = ['abund', 'name'])    # create an empty data frame to store all of the results                                    

    for sig in siglist:   
        mins  = sig.minhash.hashes
        df_singleton = pd.DataFrame(mins.values(), mins.keys(), columns = ['abund'])
        # construct a list that repeats the name the same number of times that we have mins
        # add it to the dataframe as a column
        name = sig.name
        name_lst = list()
        name_lst.append(name)
        name_lst.extend(name_lst*(len(mins.values()) - 1 ))
        df_singleton['name'] = name_lst
        df_all_sigs = df_all_sigs.append(df_singleton) # bind rows of one singleton to all singletons


    # write to a csv
    df_all_sigs.to_csv(args.output, index_label= "minhash")
        
if __name__ == '__main__':
    sys.exit(main())
