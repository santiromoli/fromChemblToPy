# This script gets information from the ChEMBL API. In this second implementation we will add to the training datasets the information
# of the test ligand, TAS-120, in order to compare the predictions. So, we only have to add four targets (FGFR1-4) that are binded by TAS-120,
# taking care not to include among them the TAS-120.
# Please, take notice that you have to install the API client first:
# >> pip install chembl-webresource-client

from chembl_webresource_client.new_client import new_client
import pandas as pd
import scipy.io
import numpy as np

tasks = [('CHEMBL1868', 'VEGFR1'), ('CHEMBL3712907', 'TMIGD3'), ('CHEMBL278', 'ITGA4'), ('CHEMBL4718', 'MNKI'),
         ('CHEMBL3130', 'PIK3CD'), ('CHEMBL5568', 'ROS1'), ('CHEMBL3594', 'CA9'), ('CHEMBL268', 'CTSK'), 
         ('CHEMBL1824', 'erb2B'), ('CHEMBL3650', 'FGFR1'), ('CHEMBL4142', 'FGFR2'), ('CHEMBL2742', 'FGFR3'), ('CHEMBL3973', 'FGFR4')]

for chembl,nam in tasks:
    # query activity API
    activities = new_client.activity.filter(target_chembl_id__in = [chembl], pchembl_value__isnull = False, 
                                            standard_type = "IC50", standard_units = 'nM', IC50_value__lte = 10000, 
                                            standar_relation__iexact = '=', assay_type = 'B'
                                            ).only(['molecule_chembl_id', 'ic50_value'])
    act_df = pd.DataFrame(activities)
    act_df = act_df.query("molecule_chembl_id == 'CHEMBL3701238'") #Since the API does not support the attribute 'exclude'
    ic50_values_df = act_df[['molecule_chembl_id','value']]

    # find the list of compounds that are within the act_df dataframe
    cmpd_chembl_ids = list(act_df['molecule_chembl_id'])

    # molecule API
    molecules = new_client.molecule.filter(molecule_chembl_id__in = cmpd_chembl_ids  
                                           ).only([ 'molecule_chembl_id', 'molecule_properties'])
    mol_df = pd.DataFrame(molecules)

    # convert nested cells (ie those containing a dictionary) to individual columns in the dataframe
    ligands = ['qed_weighted', 'hba', 'hbd', 'psa', 'cx_most_apka', 'cx_most_bpka']
    mat = {}
    num = lambda s: np.nan if s is None else float(s)
    for lig in ligands:
        mol_df[lig] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x[lig])
        mat[lig] = np.array([num(e) for e in mol_df[lig]])
    
    # post-processing: as you wish
    scipy.io.savemat('XeYt_%s.mat'%nam,mat)
