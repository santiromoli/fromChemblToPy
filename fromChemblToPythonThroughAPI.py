# This script gets information from the ChEMBL API. First, we collect all compounds that binds (or not) to CHEMBL1824, and its IC50.
# Then, we collect some properties of the ligands, such as "qed_weighted", "cx_most_apka", "cx_most_bpka", "hba", "hbd", "psa".

from chembl_webresource_client.new_client import new_client
import pandas as pd
#activity API:
activities = new_client.activity.filter(target_chembl_id__in = ['CHEMBL1824'], pchembl_value__isnull = False, 
                                        standard_type = "IC50", standard_units = 'nM', IC50_value__lte = 10000, 
                                        standar_relation__iexact = '=', assay_type = 'B'
                                        ).only(['molecule_chembl_id', 'ic50_value'])
act_df = pd.DataFrame(activities)
ic50_values_df = act_df[['molecule_chembl_id','value']]
#find the list of compounds that are within the act_df dataframe:
cmpd_chembl_ids = list(act_df['molecule_chembl_id'])
#molecule API
molecules = new_client.molecule.filter(molecule_chembl_id__in = cmpd_chembl_ids  
                                       ).only([ 'molecule_chembl_id', 'molecule_properties'])
mol_df = pd.DataFrame(molecules)
# Convert nested cells (ie those containing a dictionary) to individual columns in the dataframe
mol_df['qed_weighted'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['qed_weighted'])
mol_df['cx_most_apka'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['cx_most_apka'])
mol_df['cx_most_bpka'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['cx_most_bpka'])
mol_df['hba'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['hba'])
mol_df['hbd'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['hbd'])
mol_df['psa'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['psa'])
mol_df['cx_logd'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['cx_logd'])
mol_df['cx_logp'] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x['cx_logp'])
# post-processing: as you wish

# original code from: https://github.com/chembl/notebooks/blob/main/ChEMBL_API_example_for_webinar.ipynb
