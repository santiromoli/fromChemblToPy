# This script gets information from the ChEMBL API. First, we collect all compounds that binds (or not) to CHEMBL1824, and its IC50.
# Then, we collect some properties of the ligands, such as "qed_weighted", "cx_most_apka", "cx_most_bpka", "hba", "hbd", "psa".
# original code from: https://github.com/chembl/notebooks/blob/main/ChEMBL_API_example_for_webinar.ipynb
# >> pip install chembl-webresource-client

from chembl_webresource_client.new_client import new_client
import pandas as pd
import scipy.io
import numpy as np

tasks = [('CHEMBL1868', 'VEGFR1'), ('CHEMBL3712907', 'TMIGD3'), ('CHEMBL278', 'ITGA4'), ('CHEMBL4718', 'MNKI'),
         ('CHEMBL3130', 'PIK3CD'), ('CHEMBL5568', 'ROS1'),      ('CHEMBL3594', 'CA9'),  ('CHEMBL268', 'CTSK'), 
         ('CHEMBL1824', 'erb2B')]

for chembl,nam in tasks:
    # query activity API
    activities = new_client.activity.filter(target_chembl_id__in = [chembl], pchembl_value__isnull = False, 
                                            standard_type = "IC50", standard_units = 'nM', IC50_value__lte = 10000, 
                                            standar_relation__iexact = '=', assay_type = 'B'
                                            ).only(['molecule_chembl_id', 'ic50_value'])
    act_df = pd.DataFrame(activities)
    ic50_values_df = act_df[['molecule_chembl_id','value']]

    # find the list of compounds that are within the act_df dataframe
    cmpd_chembl_ids = list(act_df['molecule_chembl_id'])

    # molecule API
    molecules = new_client.molecule.filter(molecule_chembl_id__in = cmpd_chembl_ids  
                                           ).only([ 'molecule_chembl_id', 'molecule_properties'])
    mol_df = pd.DataFrame(molecules)
    #re-arranging the length of mol_df
    mol_df = ic50_values_df.merge(mol_df,how='right', left_on='molecule_chembl_id', right_on='molecule_chembl_id')
    # convert nested cells (ie those containing a dictionary) to individual columns in the dataframe
    ligands = ['qed_weighted', 'hba', 'hbd', 'psa', 'cx_most_apka', 'cx_most_bpka']
    mat = {}
    num = lambda s: np.nan if s is None else float(s)
    for lig in ligands:
        mol_df[lig] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x[lig])
        mat[lig] = np.array([num(e) for e in mol_df[lig]])
    #add the IC50 value to the array
    mat["IC50_value"] = np.array([num(e) for e in ic50_values_df['value']])
    # post-processing: as you wish
    scipy.io.savemat('XeYt_%s.mat'%nam,mat)
