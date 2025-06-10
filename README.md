# Design-of-Surfactant-Molecules-Under-Performance-Constraints
This repository includes the data and models developed to ensure the reproducibility of the results presented in the paper. The "property_predictor" script allows the application of these models to any molecule represented by its SMILES notation.

The heads and tails are provided as pickle files, available both as bond matrices and in SMILES notation. The final surfactant structures are also shared in SMILES format.
Each property prediction model is stored in two files:
  .joblib – Contains the trained model; and
  .joblibparameters – Lists the molecular descriptors required for model usage.
Additionally,  a summary of the methods used to develop predictive models that estimate relevant properties based on Mordred molecular descriptors is provided.
-	Synthetic Accessibility Score (SAscore): Partial least squares regression (PLS regression) is applied with 4 principal components (PCs).   
-	log(CMC): PLS regression is applied with 7 principal components (PCs). 
-	Surface tension: A PLS regression model is built, using 5 PCs.
-	Krafft point: Lasso regression with 5-fold cross validation is used, the number of non-zero parameters (i.e., descriptors) decreased to 22.
-	log(LC50): Lasso regression with 5-fold cross validation is used with 68 non-zero parameters.
-	Biodegradability: Lasso regression with 5-fold cross validation is used with 28 non-zero parameters.
