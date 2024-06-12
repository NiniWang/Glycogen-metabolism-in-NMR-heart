import os
import pandas as pd
import argparse

def merge_quant_files(folders, column_name, output_file):
    # Read the quant.sf file from each folder, selecting only the specified column, and set the column name for the DataFrame
    dataframes = []
    for folder_name in folders:
        file_path = os.path.join(folder_name, 'quant.sf')
        df = pd.read_csv(file_path, sep='\t', usecols=['Name', column_name], index_col='Name')
        dataframes.append(df)

    # Merge files and set column names for each DataFrame
    merged_data = pd.concat(dataframes, axis=1, keys=folders)

    # Save the merged data to a new file
    merged_data.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge specified column from multiple quant.sf files into a single file.")
    parser.add_argument("-folders", nargs='+', help="List of folder names containing quant.sf files", required=True)
    parser.add_argument("-col", help="Column name to extract from each quant.sf file", required=True)
    parser.add_argument("-out", help="Output file name", required=True)
    args = parser.parse_args()

    merge_quant_files(args.folders, args.col, args.out)


