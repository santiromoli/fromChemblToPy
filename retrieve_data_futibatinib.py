""" Query the ChEMBL API.
database: ChEMBL
drug: futibatinib aka TAS-120 sold as Lytgobi (anti bile duct cancer drug)
receptors (proteins of the membrane that interact with the surroundings)
  target receptors: FGFR1, FGFR2, FGFR3, FGFR4
  off-target receptors: VEGFR1, TMIGD3, ITGA4, MNKI, PIK3CD, ROS1, CA9, CTSK, erbB2
  ligands: qed_weighted, cx_most_apka, cx_most_bpka, hba, hbd, psa

# In this second implementation we will add to the training datasets the information
# of the test ligand, TAS-120, in order to compare the predictions. So, we only have to add four targets (FGFR1-4) that are 
binded by TAS-120,
# taking care not to include among them the TAS-120.

strength: IC50 (how much is needed for inhibition)

Background
The most common reason for a drug not to reach the market is toxicity. Toxicity has many issues. I wanted  to abord only one: selectivity. It is a common fact that many drugs interact with more than one target (a particular protein with which the drug interacts). In some cases, this is positive: the same illness is treated through different ways (polyphamacology). But, almost all drugs that interact with more than one target produce undesirable side effects. So, it is very useful to know (or at least, to have an insight) in the early stages of drug design, if the drug candidate will interact with more than one target (and with whom). This information is of utmost importance to leave the process of drug design (for the high number of undesirable side effects) or try to improve the drug candidate in order to improve its selectivity.


So, I thought the problem of regression like this: we will train the GPML with data of nine off-targets of TAS-120 (erbB-2, ROS1, CTSK, etc) {they are all names of receptors: proteins of the membrane that interact with the surrondings}. Then, we use ChEMBL to gather all the information of the ligands of each of this TAS-120 off-targets -that information will be the Xtraining: "qed_weighted", "cx_most_apka", "cx_most_bpka", "hba", "hbd", "psa"-. The Ytraining will be the IC50.

Then, the Xtest will be the data of TAS-120 and the Ytest the IC50 (this is the value I'm interested in: if they are something like the ones published in ChEMBL, then they are good news!).

Please install the API client first via:
>> pip install chembl-webresource-client
"""

from chembl_webresource_client.new_client import new_client
import pandas as pd
import scipy.io
import numpy as np

chem = ('CHEMBL3701238', 'FUTIBATINIB')

# 9 off-targets and 4 targets
targets = [('CHEMBL1868', 'VEGFR1'), ('CHEMBL3712907', 'TMIGD3'),
           ('CHEMBL278', 'ITGA4'),   ('CHEMBL4718', 'MNKI'),
           ('CHEMBL3130', 'PIK3CD'), ('CHEMBL5568', 'ROS1'),
           ('CHEMBL3594', 'CA9'),    ('CHEMBL268', 'CTSK'),
           ('CHEMBL1824', 'erbB2'),
           ('CHEMBL3650', 'FGFR1'),  ('CHEMBL4142', 'FGFR2'),
           ('CHEMBL2742', 'FGFR3'),  ('CHEMBL3973', 'FGFR4')
           ]
ligands = ['qed_weighted', 'hba', 'hbd', 'psa', 'cx_most_apka', 'cx_most_bpka']
ligand_names = ['QED Weighted', 'HBA', 'HBD', 'PSA', 'CX Acidic pKa', 'CX Basic pKa']

num = lambda s: np.nan if s is None else float(s)  # conversion
mat = {'ligands': np.array(ligand_names)}

for mcid,name in targets:

    print("Query target %s"%name)
    # 1) query activity API
    activities = new_client.activity.filter(target_chembl_id__in = [mcid], pchembl_value__isnull = False, 
                                            standard_type = "IC50", standard_units = 'nM', IC50_value__lte = 10000, 
                                            standar_relation__iexact = '=', assay_type = 'B'
                                            ).only(['molecule_chembl_id', 'ic50_value'])
    act_df = pd.DataFrame(activities)

    # the API does not support the attribute 'exclude'
    act_df = act_df.query("molecule_chembl_id != '%s'"%chem[0])
    ic50_values_df = act_df[['molecule_chembl_id','value']]

    # find the list of compounds that are within the act_df dataframe
    cmpd_chembl_ids = list(act_df['molecule_chembl_id'])

    # 2) query molecule API
    molecules = new_client.molecule.filter(molecule_chembl_id__in = cmpd_chembl_ids  
                                           ).only([ 'molecule_chembl_id', 'molecule_properties'])
    mol_df = pd.DataFrame(molecules)
    # re-arrange the length of mol_df
    mol_df = ic50_values_df.merge(mol_df, how='right',
                                  left_on='molecule_chembl_id',
                                  right_on='molecule_chembl_id')

    # convert nested cells (i.e. those containing a dictionary) to individual columns in the dataframe
    data = []
    for i,lig in enumerate(ligands):
        mol_df[lig] = mol_df.loc[ mol_df['molecule_properties'].notnull(), 'molecule_properties'].apply(lambda x: x[lig])
        data.append( [num(e) for e in mol_df[lig]] )
    mat[name] = np.array(data)

    # add the IC50 value to the array
    mat[name+"_IC50"] = np.array([num(e) for e in ic50_values_df['value']])
    print(ic50_values_df)
    
# save data to disk
scipy.io.savemat('XeYt.mat',mat)
