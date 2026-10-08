# -*- coding: utf-8 -*-
"""
Created on Mon Feb  2 15:22:26 2026

@author: hli490
"""

import os
import re
import sys
import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime

path = r'C:\Users\hli490\OneDrive - Vrije Universiteit Amsterdam\Desktop\CDs-to-DOs\Global paper'
df0= pd.read_csv(os.path.join(path,'data/public_emdat_custom_request_floods_diseases.csv'))
df_admin1 = pd.read_csv(os.path.join(path,'data/llmGeodis_nogeom.csv')) # adding missing admin 1; mannually checked (not all), the name is the place where the disaster was reported
df_admin1['admin1'] = df_admin1['admin1'].fillna(df_admin1['name']) 

df = pd.merge(df0,df_admin1,on='DisNo.',how="left")

df.to_csv(os.path.join(path,'data/public_emdat_custom_request_admin1.csv'))
