# DROP-DEEP
DROP-DEEP is a polygenic risk score tool that based on dimensionality reduction with principal component analysis (PCA) and a deep neural network prediction model (DNN).

steps 1-6 should be to be done for each of the chromosomes separately.
Steps 3, 5, 7 should be to be done o the training set and the test set separately. 
Step 6 is for the training and the test set together.

1.	The first step is to create the plink bed, bim, fam files, filter the samples and the SNPs according to your pre-processing parameters. Splite here your bed files to training and test files in order to train your NN. 

2.	Create covariant matrix and phenotype file.

   The covarient matrix has to be MinMax scaled. 
   The two first columns of the covarient matrix and the phenotype file has to be "FID' and "IID".

3.	If you have large data frame, you have to split your data to chanks, and to load each time around 1000 samples:
    ```
      python3 split_each_chr_to_chunks.py plink_file_name_no_suffix output_chancks_file_name
     ```

4.	Download the PCA transformers files from this link:

         https://drive.google.com/drive/folders/1oukhU_B4nM5kH9z2BxC81Kfn4kp05JAm?usp=drive_link
   
5.	 Applay our PCA transformer on your data:
   
       ```
      python3 PCA_dimension_reduction.py chanks_file output_pca_file pca_transformer
  	  ```

7.	Scale the PCA data (MinMax scale):
      ```
      python3 scale_genes.py PCA_file_training_set PCA_file_test_set scaled_file_training_set scaled_file_test_set
      ```

8.	Join the 22 chromosomes to one cohort:
      ```
      Python3 join_all_chr_after_dr.py cov_file input_files_dir merged_output_file_name
      ```

7.	Match the phenotype files (the Y files) to the feature files (the X files).
      ```
      Python3 match_pheno_IDs.py pheno_file train_X_files test_X_files output_path
      ```
9.	Run NN on the PCA features.

   for binary phenotypes:
      ```
      python3 NN_for_binary_pheno.py
      ```
    for continues phenotypes:
    ```
      python3 NN_for_linear_pheno.py
    ```
