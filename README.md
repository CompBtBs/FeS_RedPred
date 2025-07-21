# FeS_MLmodel

This repository contains:
* 2 input files:
	- Database Redox Pot Fe2S2 proteins.xlsx
	- tableAmm.txt, a utility file with amino acids' parametrization and list of other cofactors to be counted by features_calculator.py
* 2 utility modules:
	- utils.py, with functions used by both features_calculator.py and em_predict2.py
	- ML_models.py a dictionary for models with hyperparameters grid which will be tuned with GridSearch optimization
* pdb-files folder, which contains all pdb files of dataset's proteins, including in silico generated mutants
* Folder A contains the scripts for training seprate models for each specific combination of radius values r1 and r2:
	- features_calculator.py script uded to compute molecular descriptors values. These descriptors are saved in a dataset_features_r1_r2.xlsx file which serves as input for model training.
	- em_predict.py, the main script used to launch models training and to test their performance
* Folder B includes the code for training a single model that simultaneously considers all features calculated for every r1 and r2:
	- features_calculator.py
	- total.py to merge all features in one single file total.xlsx, avoiding repetitions
	-em_predict.py

The remaining models were constructed using the scripts in folder A and modifying the features_calculators.py output files removing the selected features.

Warning: we run all codes in linux, when running features_calculator.py in windows a modification on PDBParser library is needed (l.192 resname = line[17:20].replace(' ',''))
