#this script can be used to merge all features in one single file total.xlsx, avoiding repetitions

import os
import logging
from tqdm import tqdm
from collections import OrderedDict
from Bio.PDB import PDBList, PDBParser
#from PyBioMed.PyProtein import CTD
import numpy as np
import pandas as pd


list_bar_radii= list (np.arange(8,17))
list_CofAtom_radii = list (np.arange(3,6))
list_CofAtom_names = ['FE1', 'FE2', 'S1', 'S2', 'FE']
dir = 'dataset_features/'
df_total = pd.DataFrame()


df_protein = pd.read_excel(f'{dir}dataset_protein_8_3.xlsx')
df_ini = df_protein[["PDB Code(s)","mV","pH","n_chains"]]
df_total = pd.concat([df_total, df_ini], axis=1)
protein_columns = [col for col in df_protein.columns if 'Protein.' in col]
filtered_df_protein = df_protein[protein_columns]
df_total = pd.concat([df_total, filtered_df_protein], axis=1)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)




for CofAtom in list_CofAtom_radii:
    dataset_CofAtom = pd.read_excel(f'{dir}dataset_protein_8_{str(CofAtom)}.xlsx')

    logger.info(f'Processing dataset_protein_8_{str(CofAtom)}.xlsx')

    CofAtom_columns = [col for col in dataset_CofAtom.columns if 'CofAtom.' in col]
    filtered_dataset_CofAtom = dataset_CofAtom[CofAtom_columns].copy()
    filtered_dataset_CofAtom.columns = [f'{str(CofAtom)}_{col}' for col in filtered_dataset_CofAtom.columns]
    df_total = pd.concat([df_total, filtered_dataset_CofAtom], axis=1)

    for CofAtom_name in list_CofAtom_names:

        Cof_Oxigen_columns = [col for col in dataset_CofAtom.columns if f'{CofAtom_name}.Oxigen' in col]
        filtered_Cof_oxigen = dataset_CofAtom[Cof_Oxigen_columns].copy()
        filtered_Cof_oxigen.columns = [f'{str(CofAtom)}_{col}' for col in filtered_Cof_oxigen.columns]
        df_total = pd.concat([df_total, filtered_Cof_oxigen], axis=1)

        Cof_nitrogen_columns = [col for col in dataset_CofAtom.columns if f'{CofAtom_name}.Nitrogen' in col]
        filtered_Cof_nitrogen = dataset_CofAtom[Cof_nitrogen_columns].copy()
        filtered_Cof_nitrogen.columns = [f'{str(CofAtom)}_{col}' for col in filtered_Cof_nitrogen.columns]
        df_total = pd.concat([df_total, filtered_Cof_nitrogen], axis=1)

        Cof_carbon_columns = [col for col in dataset_CofAtom.columns if f'{CofAtom_name}.Carbon' in col]
        filtered_Cof_carbon = dataset_CofAtom[Cof_carbon_columns].copy()
        filtered_Cof_carbon.columns = [f'{str(CofAtom)}_{col}' for col in filtered_Cof_carbon.columns]
        df_total = pd.concat([df_total, filtered_Cof_carbon], axis=1)

        Cof_sulfur_columns = [col for col in dataset_CofAtom.columns if f'{CofAtom_name}.Sulfur' in col]
        filtered_Cof_sulfur = dataset_CofAtom[Cof_sulfur_columns].copy()
        filtered_Cof_sulfur.columns = [f'{str(CofAtom)}_{col}' for col in filtered_Cof_sulfur.columns]
        df_total = pd.concat([df_total, filtered_Cof_sulfur], axis=1)

        Cof_nearest_columns = [col for col in dataset_CofAtom.columns if f'{CofAtom_name}_nearest.' in col]
        filtered_Cof_nearest = dataset_CofAtom[Cof_nearest_columns].copy()
        filtered_Cof_nearest.columns = [f'{str(CofAtom)}_{col}' for col in filtered_Cof_nearest.columns]
        df_total = pd.concat([df_total, filtered_Cof_nearest], axis=1)

        Cof_around_columns = [col for col in dataset_CofAtom.columns if f'Around_{CofAtom_name}' in col]
        filtered_Cof_around= dataset_CofAtom[Cof_around_columns].copy()
        filtered_Cof_around.columns = [f'{str(CofAtom)}_{col}' for col in filtered_Cof_around.columns]
        df_total = pd.concat([df_total, filtered_Cof_around], axis=1)



for bar in list_bar_radii:
    dataset_bar = pd.read_excel(f'{dir}dataset_protein_{str(bar)}_3.xlsx')
    logger.info(f'Processing dataset_protein_{str(bar)}_3.xlsx')
    bar_columns = [col for col in dataset_bar.columns if 'Bar.' in col]
    filtered_dataset_bar = dataset_bar[bar_columns].copy()
    filtered_dataset_bar.columns = [f'{str(bar)}_{col}' for col in filtered_dataset_bar.columns]
    df_total = pd.concat([df_total,filtered_dataset_bar], axis=1)


logger.info('Saving total.xlsx')

df_total.to_excel(f'{dir}total.xlsx')
    
