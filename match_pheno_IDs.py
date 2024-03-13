import pandas as pd
import numpy as np
import pickle
import sys
import os

def match(x_file_name, pheno_df, x_output_file):   
    number_samples = 0
    num_batch = 1 
    with open(x_output_file, 'wb') as x_output_file_handle:
        with open(x_file_name, 'rb') as x_file_handle:
            while True:
                try:
                    X_batch = pickle.load(x_file_handle)
                    X_batch = X_batch[X_batch['FID'].isin(pheno_df['FID'])].reset_index(drop=True)
                    print('number samples in X_batch =', len(X_batch))
                    if num_batch==1:
                        print('number samples in validation =',len(X_batch))
                    num_batch = num_batch + 1
                    pickle.dump(X_batch, x_output_file_handle, protocol=4)
                    number_samples = number_samples + len(X_batch)
                except EOFError:
                    print('number_samples =',number_samples)
                    break
                
def split_y_train_to_chunks(y_output_file, pheno_df, x_chunk_file, pheno_name):
    with open(y_output_file, 'wb') as y_file_handle:
       with open(x_chunk_file, 'rb') as x_file_handle:
           while True:
               try:
                   X_batch = pickle.load(x_file_handle)
                   df = pd.merge(X_batch, pheno_df, on='FID') #in order to order samples in the same order as in genes
                   #print(df.head(1))
                   #print(df.iloc[:,[0,-1]])
                   y_batch = df.iloc[:,[0,-1]].reset_index(drop=True)
                   print('number samples in y_batch =', len(y_batch))
                   pickle.dump(y_batch, y_file_handle, protocol=4) 
               except EOFError:
                   break
                               
def split_test_y_to_chunks(input_df, x_chunk_file, pheno_name):
    pheno_df = pd.DataFrame()
    with open(x_chunk_file, 'rb') as x_file_handle:
        while True:
            try:
                X_batch = pickle.load(x_file_handle)
                df = pd.merge(X_batch, input_df, on='FID') #in order to order samples in the same order as in genes
                print(df.head(1))
                print(df.iloc[:, [0, -1]])
                y_batch =df.iloc[:,[0,-1]].reset_index(drop=True)
                pheno_df = pd.concat([pheno_df, y_batch], ignore_index=True)
            except EOFError:
                return pheno_df

pheno_name=sys.argv[1]
rep=sys.argv[2]
type=sys.argv[3]
DRM = sys.argv[4]
# handle phenotypes dataframe

pheno_file=pheno_name.lower()


pheno = pd.read_csv("/sise/nadav-group/nadavrap-group/UKBB/phenotypes/"+pheno_file, sep='\t')
if (pheno.shape[1]==1):
        pheno = pd.read_csv("/sise/nadav-group/nadavrap-group/UKBB/phenotypes/"+pheno_file, sep=' ')
    
#pheno = pd.read_csv("/sise/nadav-group/nadavrap-group/UKBB/phenotypes/"+pheno_file)
pheno.rename(columns={"#FID":"FID"},inplace=True)

print(pheno)

pheno['FID']=pheno['FID'].astype(str)
pheno = pheno.iloc[:,[0,2]]
pheno = pheno.dropna(axis=0).reset_index(drop=True)
print(pheno)

if type=='b':

    print(pheno.iloc[:, 1])
    pheno.iloc[:, 1] = pheno.iloc[:, 1].replace([1], "0")
    pheno.iloc[:, 1] = pheno.iloc[:, 1].replace([2], "1")
    print(pheno)


if DRM == "PCA":

    path = '/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/PCA/'
    #train
    union_train_gene =path+'X_train_1k_chunks_PCA_dim_remove_no_missing/rep'+rep+"/X_train_all_chr_MinMax_cov_MinMax.pkl"
    x_train_output_file = path+'X_train_1k_chunks_PCA_dim_remove_no_missing/rep'+rep+"/"+pheno_name+"_X_train_match_to_pheno_MinMax_cov_MinMax.pkl"
    # test
    union_test_gene = path+'X_test_1k_chunks_PCA_dim_remove_no_missing/rep' + rep + "/X_test_all_chr_MinMax_cov_MinMax.pkl"
    x_test_output_file = path+"X_test_1k_chunks_PCA_dim_remove_no_missing/rep" + rep + "/" + pheno_name + "_X_test_match_to_pheno_MinMax_cov_MinMax.pkl"

elif DRM == "Autoencoder":

    path = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/Autoencoder/"
    #train
    union_train_gene = path + "X_train_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/X_train_all_chr_MinMax_cov_MinMax.pkl"
    x_train_output_file = path + "X_train_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/"+pheno_name+"_X_train_match_to_pheno_MinMax_cov_MinMax.pkl"
    # test
    union_test_gene = path+ "X_test_1k_chunks_dim_remove_no_missing_500_epochs/rep" + rep + "/X_test_all_chr_MinMax_cov_MinMax.pkl"
    x_test_output_file = path+ "X_test_1k_chunks_dim_remove_no_missing_500_epochs/rep" + rep + "/" + pheno_name + "_X_test_match_to_pheno_MinMax_cov_MinMax.pkl"

y_train_output_file = path+"Y_files/rep"+rep+"/"+pheno_name+"_Y_train_1k_chunks_no_missing.pkl"
y_test_output_file = path+"Y_files/rep"+rep+"/"+pheno_name+"_Y_test_1k_chunks_no_missing.pkl"


print('train')
match(union_train_gene, pheno, x_train_output_file)

##create Y_files directory
if not os.path.exists(
    path+"Y_files/rep"+rep):
    os.makedirs(
         path+"Y_files/rep"+rep)
else:
    print( "Y files directory existed for rep"+rep)


 #--------------------y train--------------------
split_y_train_to_chunks(y_train_output_file, pheno, x_train_output_file, pheno_name)
#----------------------------------------

print('test')
match(union_test_gene, pheno, x_test_output_file)

#--------------------y test--------------------
test_pheno_df = split_test_y_to_chunks(pheno, x_test_output_file, pheno_name)
test_pheno_df = test_pheno_df.reset_index(drop=True)
test_pheno_df.to_pickle(y_test_output_file)
#----------------------------------------
