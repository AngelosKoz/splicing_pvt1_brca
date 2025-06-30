#!/bin/bash
import pandas as pd
import argparse


def get_splicing_eff(file, freedom=10,):
    column_names = ['chr', 'chrStart', 'chrEnd', 'transcript', 'form1', 'strand', 'splitReads', 'nonSplitReads']
    DF = pd.read_csv(str(file), sep='\t', header=None, names=column_names, index_col=False)

    # Calculate the sum of the absolute values of column 7 and column 8 (split reads, non split reads respectively)
    DF['readSum'] = DF['splitReads'] + DF['nonSplitReads']
    
    # Keep only splice sites with total reads >= freedom (set threshold: 0 to keep all)
    filtered_DF = DF.loc[DF['readSum'] >= freedom].copy()
    filtered_DF['splice_eff'] = filtered_DF['splitReads'] / filtered_DF['readSum']
    filtered_DF['transcript'] = filtered_DF['transcript'].str.replace(';', '')
    
    
    # Generate a sequence number for duplicates, starting from 1
    filtered_DF['dupe_count'] = filtered_DF.groupby('transcript').cumcount() + 1
    # Check for duplicates
    filtered_DF['is_dupe'] = filtered_DF.duplicated('transcript', keep=False)

    # Append '.N' to all, including non-duplicates if they appear more than once in the dataset (.N = sequence number)
    filtered_DF['unique_name'] = filtered_DF.apply(
        lambda row: f"{row['transcript']}.{row['dupe_count']}" 
                                 if row['is_dupe'] 
                                 else row['transcript'], axis=1
                                 )
    filtered_DF.drop(columns=['dupe_count', 'is_dupe', 'transcript'], inplace=True)
    filtered_DF.reset_index(drop=True, inplace=True)
    column_order = ['chr', 'chrStart', 'chrEnd', 'unique_name', 'form1', 'strand', 'splitReads', 'nonSplitReads', 'readSum', 'splice_eff']
    filtered_DF = filtered_DF[column_order]
    
    return filtered_DF



def main(input_file, output_file):
    filtered_DF = get_splicing_eff(input_file)
    filtered_DF.to_csv(output_file, sep='\t', index=False, header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate Splicing Efficiency',
        epilog='''Use the name of the file as input using quotes and besides the output of the file.
	Example --> "python3 get_splicing_eff.py 'initial_file_name' 'output_file_name'".
	 -- In column 7, the split reads are expected, and in column 8, the non-split reads.
	 -- Two additional columns are added to the final output: readSum (split+non-split) and splice_eff -- θ = split / (split+non-split reads)''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('input_file', type=str, help='Input file path')
    parser.add_argument('output_file', type=str, help='Output file path')
    parser.add_argument('--freedom', type=int, default=10, help='Cutoff value for readSum (default: 10)')
    args = parser.parse_args()

    main(args.input_file, args.output_file)


