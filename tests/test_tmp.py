#!python
#

#Temp test

import os
import re
import json
from pathlib import Path
import itertools
'''
import dfttk.structure_builders.sqs_db as sqs_db

atat_sqsdb_path = "/storage/home/mjl6505/package/atat/atat/data/sqsdb"

'''
syms = ['Fe', 'Fe', 'Cr', 'Fe', 'Cr', 'Cr', 'Ni', 'Fe', 'Ni', 'Ni']

syms_group = [a[0] for a in itertools.groupby(syms)]
print(syms_group)