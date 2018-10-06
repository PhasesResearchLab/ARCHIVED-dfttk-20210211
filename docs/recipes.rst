

=======
Recipes
=======

Make ESPEI datasets from a QHA database, optionally writing the (nicely named) files to dict.
The QHA database requires the following metadata schema:

.. code-block:: json

   {
     'metadata': {
       'phase_name': 'FCC_A1',
       'tag': 'ed447049-ad67-4090-ba99-378188d3416b',
       'sublattice': {
         'configuration': [['Cr', 'Ni']],
         'occupancies': [[0.03125, 0.96875]]
        },
     }
   }

.. code-block:: python

   ########
   # EDIT #
   ########

   phase_name = 'BCC_A2'
   configuration_to_find = [['Ni', 'V']]
   sublattice_site_ratios = [1.0]
   db_username = 'BrandonResultsRO'
   db_password = 'piqhg38hap3'
   db_uri = 'mongodb://206.189.190.225:27018'
   WRITE_FILES = True

   temperature_index = 59  # index of 300 K temperature (close to 298 K), found by hand

   refstate_tags = {
       'Fe': '4ac77fce-0e43-4c07-8418-4843a2cd5723',
       'Ni': '0059ee69-4a8f-4e86-9895-9b40cf67dd96',
       'Cr': '8cceb186-2796-4488-ba8c-2380c5278f62',
       'V': 'fba46b6b-1699-419f-b5e1-da9533530701',
   }

   ########################
   ### SCRIPT
   ########################

   from dfttk.analysis.formation_energies import get_formation_energy, get_thermal_props
   from dfttk.espei_compat import make_dataset, dfttk_config_to_espei, dfttk_occupancies_to_espei
   from pymatgen import Structure
   import numpy as np
   from dfttk.utils import recursive_flatten
   import json
   from pymongo import MongoClient

   # create the MongoClient
   cli = MongoClient(db_uri)
   db = cli.results
   db.authenticate(name=db_username, password=db_password)
   coll = db.qha

   # construct the energies, assupmtion of same temperature grid
   # energies are J/mol-atom
   refstate_energies = {}
   for el, tag in refstate_tags.items():
       qha_result = coll.find_one({'metadata.tag': tag})
       refstate_energies[el] = get_thermal_props(qha_result)

   # calculate all the dataset values
   configs     = []
   occupancies = []
   hm_values   = []
   sm_values   = []
   cpm_values  = []

   # we'll change the T to the right temperatures later
   fixed_conds = {'P': 101325, 'T': 0}
   temp_conds = {'P': 101325, 'T': 0}


   for qha_res in coll.find({'metadata.sublattice.configuration': configuration_to_find, 'metadata.phase_name': phase_name}):
       configs.append(qha_res['metadata']['sublattice']['configuration'])
       occupancies.append(qha_res['metadata']['sublattice']['occupancies'])

       tprops = get_thermal_props(qha_res)
       struct = Structure.from_dict(qha_res['structure'])
       hm_form = get_formation_energy(tprops, struct, refstate_energies, 'HM', idx=temperature_index)
       sm_form = get_formation_energy(tprops, struct, refstate_energies, 'SM', idx=temperature_index)
       cpm_form = get_formation_energy(tprops, struct, refstate_energies, 'CPM', thin=10)[:-2]
       fixed_temp = tprops['T'][temperature_index]
       cpm_temps = tprops['T'][::10][:-2]

       hm_values.append(hm_form)
       sm_values.append(sm_form)
       cpm_values.append(cpm_form)

   fixed_conds['T'] = fixed_temp.tolist()
   temp_conds['T'] = cpm_temps.tolist()

   # make the HM, SM, CPM values arrays of the proper shape
   hm_values = np.array([[hm_values]])
   sm_values = np.array([[sm_values]])
   cpm_values = np.array(cpm_values).T[np.newaxis, ...]

   if WRITE_FILES:
       # write JSON files
       comps = [c.upper() for c in sorted(recursive_flatten(configuration_to_find))]
       for prop, vals, conds in [('HM_FORM', hm_values, fixed_conds), ('SM_FORM', sm_values, fixed_conds), ('CPM_FORM', cpm_values, temp_conds)]:
           ds = make_dataset(phase_name, prop, sublattice_site_ratios, configs, conds, vals, occupancies=occupancies, tag=tag)
           with open('{}-{}-{}-DFTTK.json'.format('-'.join(comps), phase_name, prop), 'w') as fp:
               json.dump(ds, fp)

