import pandas as pd
from sklearn.preprocessing import StandardScaler
import sys
import pickle5 as pickle


def join(X_file_name_beginning, chunk_size, gene_output_file, cov_df,X_file_name_end):
    number_samples = 0 
    with open(gene_output_file, 'wb') as gene_output_file_handle:
        this_chunk = 1
        while True:
            try:
                for chr in range(1, 23):
                    # if chr <= 4:
                    #     chunk_size=500
                    # else:
                    #     chunk_size=1000
                    print(chr)
                    dim_remove_file = X_file_name_beginning + str(chr) +  X_file_name_end
                    #dim_remove_file = dim_remove_file
                    with open(dim_remove_file, 'rb') as dim_remove_file_handle:
                        for c in range(1,this_chunk+1):
                            print(c)
                            print(pickle.load(dim_remove_file_handle))
                            batch = pickle.load(dim_remove_file_handle)
                        ID = batch['FID']
                        batch = batch.drop(['FID'], axis=1)
                        suffix = '_chr'+str(chr)
                        batch = batch.add_suffix(suffix)
                        batch['FID'] = ID
                        if chr == 1:
                            df = pd.merge(cov_df, batch, on='FID')
                        else:
                            df = pd.merge(df, batch, on='FID')
                        
                pickle.dump(df, gene_output_file_handle, protocol=4)
                this_chunk = this_chunk + 1
                number_samples = number_samples + len(df)
            except EOFError:
                print('error: '+e)
                print('number_samples after match cov =',number_samples)
                break


                
 
# handle covariate
DRM=sys.argv[1]
rep=sys.argv[2]

with open('/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/cov_matrix/rep1'
          +'/cov_matrix_MinMax_scaled_no_missing.pkl', "rb") as fh:
#with open('/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/cov_matrix/rep'+rep
#          +'/cov_matrix_MinMax_scaled_no_missing.pkl', "rb") as fh:
  cov_df = pickle.load(fh)
#cov_df = pd.read_pickle('/cs/snapless/michall/hadasak/improve_PRS/cov_matrix/cov_matrix_MinMax_scaled_no_missing.pkl')
cov_df.rename(columns={"#FID":"FID"},inplace=True)
print(cov_df.head(1))
chunk_size=1000

if DRM=="PCA":
    X_train_name_beginning = '/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/PCA/X_train_1k_chunks_PCA_dim_remove_no_missing/rep'+rep+"/chr_"
    X_test_name_beginning = '/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/PCA/X_test_1k_chunks_PCA_dim_remove_no_missing/rep'+rep+"/chr_"
    union_train_gene_output_file = '/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/PCA/X_train_1k_chunks_PCA_dim_remove_no_missing/rep' + rep + "/X_train_all_chr_MinMax_cov_MinMax.pkl"
    union_test_gene_output_file = '/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/PCA/X_test_1k_chunks_PCA_dim_remove_no_missing/rep'+rep+"/X_test_all_chr_MinMax_cov_MinMax.pkl"

if DRM=="Autoencoder":
  
    X_train_name_beginning = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/Autoencoder/X_train_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/chr_"
    X_test_name_beginning = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/Autoencoder/X_test_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/chr_"
    union_train_gene_output_file = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/Autoencoder/X_train_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/X_train_all_chr_MinMax_cov_MinMax.pkl"
    union_test_gene_output_file = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/Autoencoder/X_test_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/X_test_all_chr_MinMax_cov_MinMax.pkl"

file_name_end="_MinMax_scaled.pkl"

#train
print('train')
join(X_train_name_beginning, chunk_size, union_train_gene_output_file, cov_df,file_name_end)

#test
print('test')
join(X_test_name_beginning, chunk_size, union_test_gene_output_file, cov_df,file_name_end)
