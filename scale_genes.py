import pickle
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pandas as pd
import sys

            
def StandardScaler_incremental(train_file_name_to_standard_scale, test_file_name_to_standard_scale):
    scaler = {}
    for chr in range(1, 23):
        dim_remove_file = train_file_name_to_standard_scale + str(chr) + ".pkl"
        with open(dim_remove_file, 'rb') as file_handle:
            #chr_scaler = StandardScaler()
            chr_scaler = MinMaxScaler()
            while True:
                try:
                    batch = pickle.load(file_handle)
                    batch = batch.drop(['FID'], axis=1)
                    chr_scaler.partial_fit(batch)
                except EOFError:
                    break
            scaler[chr] = chr_scaler
    scale_all_genes(train_file_name_to_standard_scale, scaler)
    scale_all_genes(test_file_name_to_standard_scale, scaler)
  
def scale_all_genes(file_name_to_standard_scale, scaler):
    for chr in range(1, 23):
        file_name = file_name_to_standard_scale + str(chr) + ".pkl"
        with open(file_name, 'rb') as file_handle:
            #output_file_name = file_name_to_standard_scale + str(chr) + "_scaled.pkl"
            output_file_name = file_name_to_standard_scale + str(chr) + "_MinMax_scaled.pkl"
            with open(output_file_name, 'wb') as gene_output_file_handle:
                while True:
                    try:
                        batch = pickle.load(file_handle)
                        ID = batch['FID']
                        batch = batch.drop(['FID'], axis=1)
                        batch_scaled = scaler[chr].transform(batch)
                        batch = pd.DataFrame(batch_scaled)
                        batch['FID'] = ID
                        pickle.dump(batch, gene_output_file_handle, protocol=4)
                    except EOFError:
                        break


rep = sys.argv[1]
DRM = sys.argv [2]

if DRM == "PCA":
    train_X_file_name = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/PCA/X_train_1k_chunks_PCA_dim_remove_no_missing/rep" + rep + '/chr_'  # for PCA
    test_X_file_name = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/PCA/X_test_1k_chunks_PCA_dim_remove_no_missing/rep"+rep+'/chr_'   # for PCA
else:
    if DRM == "Autoencoder":
        train_X_file_name = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/Autoencoder/X_train_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/chr_" #for Autoencoder
        test_X_file_name = "/sise/nadav-group/nadavrap-group/hadasa/my_storage/impoving_PRS/data/our_model/Autoencoder/X_test_1k_chunks_dim_remove_no_missing_500_epochs/rep"+rep+"/chr_"  # for Autoencoder
    else:
        print("\n\nerror! you nust chose PCA or Autoencoder as diamentialy reduction model!!!\n\n")
#test
StandardScaler_incremental(train_X_file_name, test_X_file_name)

