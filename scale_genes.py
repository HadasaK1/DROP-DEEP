import pickle
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pandas as pd
import sys

            
def StandardScaler_incremental(PCA_file_training_set,PCA_file_test_set,scaled_file_training_set,scaled_file_test_set):
    #Scaler adjustment for scaling to both files based on the training set
    with open(PCA_file_training_set, 'rb') as file_handle:
        scaler = MinMaxScaler()
        while True:
            try:
                batch = pickle.load(file_handle)
                batch = batch.drop(['FID'], axis=1)
                scaler.partial_fit(batch)
            except EOFError:
                break
    scale_all_genes(PCA_file_training_set,scaled_file_training_set, scaler)
    scale_all_genes(PCA_file_test_set,scaled_file_test_set, scaler)
  
def scale_all_genes(file_name_to_standard_scale, output_scaled_file,scaler):
    file_name = file_name_to_standard_scale
    with open(file_name, 'rb') as file_handle:
        output_file_name = output_scaled_file
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


PCA_file_training_set = sys.argv[1]
PCA_file_test_set = sys.argv[2]
scaled_file_training_set = sys.argv[3]
scaled_file_test_set = sys.argv[4]
StandardScaler_incremental(PCA_file_training_set,PCA_file_test_set,scaled_file_training_set,scaled_file_test_set)

