# DROP-DEEP
DROP-DEEP is a polygenic risk score tool that based on dimensionality reduction with principal component analysis (PCA) and a deep neural network prediction model (DNN).

1.	The first step is to create the plink bed, bim, fam files, filter the samples and the SNPs according to your pre-processing parameters. Splite here your bed files to training and test files in order to train your NN. 

2.	Create covariant matrix.

3.	If you have large data frame, you have to split your data to chanks, and to load each time around 1000 samples:

python3 split_each_chr_to_chunks.py input_file_name output_file_name

4.	Applay our PCA transformer on your data:
   
python3 PCA_dimension_reduction.py

6.	Scale the PCA data (MinMax scale):
   
python3 scale_genes.py PCA

7.	Join the 22 chromosomes to one cohort:
   
Python3 join_all_chr_after_dr.py PCA

9.	Validate that the samples in your phenotypes file are exactly the same in the features file.

10.	Run NN on the PCA features.
   
python3 NN.py
python3 NN_for_binary_pheno.py
