#%%
# -*- coding: utf-8 -*-
import joblib
from mordred import Calculator, descriptors
from rdkit import Chem
import warnings
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer

def find_mordred_descriptors(smiles_list):
    clean_list_of_smiles = [smi for smi in smiles_list if Chem.MolFromSmiles(smi,sanitize=False) != None]
    df = []
    calc = Calculator(descriptors, ignore_3D=True)
    mols = [Chem.MolFromSmiles(i) for i in clean_list_of_smiles] 
    df=calc.pandas(mols)
    return df,mols,clean_list_of_smiles

def run(models, smiles_list):
    warnings.simplefilter(action='ignore', category = FutureWarning)
    mordred_database = find_mordred_descriptors(smiles_list)
    all_results = pd.DataFrame()
    for mt in models.keys():
        if mt == "sascore":
            name= "linear_SAscore"
        elif mt == "LC50":
            name = "logarithmic_LC50"
        elif mt == "CMC":
            name = "logarithmic_CMC"
        elif mt == "biodegradability":
            name = "linear_Biodeg"
        elif mt == "SurfaceTension":
            name= "linear_SurfaceTension"
        elif mt == "KrafftPoint":
            name = "linear_KrafftPoint"
        else:
            print("Please pick a model from the given database. There is no model "+mt+".")
            break
 
        model = joblib.load(name+".joblib")[0]

        mordred_descriptors, scale_min, scale_max = joblib.load(name+".joblib_parameters") 
        
        filter_md = mordred_database[mordred_descriptors]
        filter_md = filter_md.append(scale_max, ignore_index=True)
        filter_md = filter_md.append(scale_min, ignore_index=True)
        filter_md_array=np.array(filter_md).T
        for i, var in enumerate(filter_md_array): 
            var = (var-var[-1])/(var[-2]-var[-1])
            filter_md_array[i] = var
        x_array = filter_md_array.T
        df=pd.DataFrame(x_array)
        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        x_array2=imp.fit_transform(df)
        x_predict = x_array2[:-2,:]
        if (type(model) is tuple) == True:
            model=model[0]
            results = model.predict(x_predict)        
        elif "pls" in name:
            results = model.predict(x_predict)
            results=np.array([i[0] for i in results])            
        else:
            results = model.predict(x_predict)       
        if "logarithmic" in name:
            results = 10**results
        all_results[mt]=results       
    return all_results

models=["sascore"," LC50","CMC","biodegradability","SurfaceTension","KrafftPoint"]
smiles=["CCCCCCCCCCCCCOC(=O)O", "CCCCCCCCCCCCOOC(=O)O", "CCCCCCCCCCCC(=O)NC(C)=O", "CCCCCCCCCCCCNC(=O)O"]
properties=run(models,smiles)