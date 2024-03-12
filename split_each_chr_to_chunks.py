import pandas as pd
import pickle
from pandas_plink import read_plink, read_plink1_bin
import numpy as np
import pickle5
import os
import sys


def split_X_to_chunks(input_file, output_path, chunk_size):
    from_ = 0
    to_ = chunk_size

    with open(output_file, 'wb') as gene_file_handle:
        file_name = input_file
        print(file_name)
        G = read_plink1_bin(file_name, verbose=False)  # xarray.core.dataarray.DataArray
        while True:  # while there are still samples
            G_to_pandas = G[from_:to_, :].to_pandas()  # pandas.core.frame.DataFrame
            if G_to_pandas.empty == False:  # there's more samples
                (bim, fam, bed) = read_plink(file_name, verbose=False)  # bed is dask.array.core.Array
                snp_columns = bim['snp'].to_numpy()
                G_to_pandas.set_axis(snp_columns, axis=1, inplace=True)
                G_to_pandas = G_to_pandas.fillna(0)
                G_to_pandas = G_to_pandas.astype(np.int8)
                G_to_pandas.reset_index(level=0, inplace=True)
                G_to_pandas.rename(columns={'sample': 'FID'}, inplace=True)
                pickle.dump(G_to_pandas, gene_file_handle, protocol=4)
                from_ = to_
                to_ = to_ + chunk_size
            else:
                print("finished")
                return


def main():
    file_to_split = sys.argv[1]
    output_file = sys.argv[2]
    chunk_size = 1000

    split_X_to_chunks(input_file=file_to_split,
                     output_file=output_file,
                     chunk_size=chunk_size)
]
if __name__ == "__main__":
    main()

