# DROP-DEEP
DROP-DEEP is a polygenic risk score tool that based on dimensionality reduction with principal component analysis (PCA) and a deep neural network prediction model (DNN).

steps 1-7 should be to be done for each of the chromosomes separately.
Steps 3 and 5 should be to be done o the training set and the test set separately.

1.	The first step is to create the plink bed, bim, fam files, filter the samples and the SNPs according to your pre-processing parameters. Splite here your bed files to training and test files in order to train your NN. 

2.	Create covariant matrix.

  	the covarient matrix has to be MinMax scaled. 
   The two first columns has to be "FID' and "IID".

4.	If you have large data frame, you have to split your data to chanks, and to load each time around 1000 samples:

      python3 split_each_chr_to_chunks.py plink_file_name_no_suffix output_chancks_file_name

5.	Download the PCA transformers files from this link:

         https://drive.google.com/drive/folders/1oukhU_B4nM5kH9z2BxC81Kfn4kp05JAm?usp=drive_link
   
6.	 Applay our PCA transformer on your data:
   
      python3 PCA_dimension_reduction.py chanks_file output_pca_file pca_transformer

7.	Scale the PCA data (MinMax scale):
   
      python3 scale_genes.py PCA_file_training_set PCA_file_test_set scaled_file_training_set scaled_file_test_set

8.	Join the 22 chromosomes to one cohort:
   
      Python3 join_all_chr_after_dr.py cov_file input_files_dir merged_output_file_name

9.	Validate that the samples in your phenotypes file are exactly the same in the features file.

10.	Run NN on the PCA features.
   
      python3 NN.py
      python3 NN_for_binary_pheno.py
