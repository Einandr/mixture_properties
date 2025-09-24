import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy import interpolate
import math
import ast
import utils.terra_component as tc
import utils.mixture_component as dc
import utils.material as material
import os
import yaml
import re
from chempy import balance_stoichiometry
from chempy import mass_fractions
from chempy import chemistry
from pprint import pprint
from enum import Enum, auto
from pathlib import Path
from configparser import ConfigParser


path = r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS'
dir_run = 'RUN'
dir_props_out = 'props_out'
dir_soulution = 'solution'

# path_prop = r'D:\YASIM\VORON\2024_01_All_Discrete_Phase_Development\kosh\new_version'
path_prop = os.getcwd()
file_props = 'props_full_database.txt'
file_config = 'run_config.ini'
file_MPL = 'MPL_constants.txt'


path_run = ''.join((path, '/', dir_run))
Path(path_run).mkdir(parents=True, exist_ok=True)
os.chdir(path_run)


df_props = pd.read_csv(''.join((path_prop, '\\', file_props)), delimiter=' ', index_col='name')
df_props['T_range'] = df_props['T_range'].apply(lambda x: np.fromstring(x.replace('\'', '').replace('\"', '').replace('\n', ''), dtype=float, sep=','))
df_props['f_ranges'] = df_props['f_ranges'].apply(lambda x: np.array(ast.literal_eval(x.replace('\'', '').replace('\"', '').replace('\n', ''))))


g_O2 = {'O2': 1}
g_N2 = {'N2': 1}
g_CO2 = {'CO2': 1}
g_CO = {'CO': 1}
g_H2O = {'H2O': 1}
g_C6H14 = {'C6H14': 1}
g_C10H22 = {'C10H22': 1}
g_C6H6 = {'C6H6': 1}
g_C7H16 = {'C7H16': 1}
g_C9H12 = {'C9H12': 1}
g_C9H18 = {'C9H18': 1}
g_BHD = {'C6H14': 0.091,
         'C10H22': 0.727,
         'C6H6': 0.182}
g_KERO = {'C9H12': 0.132,
          'C10H22': 0.767,
          'C9H18': 0.101}



gas_mixture_reference = {'O2': g_O2,
                         'N2': g_N2,
                         'CO2': g_CO2,
                         'CO': g_CO,
                         'H2O': g_H2O,
                         'C6H14': g_C6H14,
                         'C10H22': g_C10H22,
                         'C6H6': g_C6H6,
                         'C9H12': g_C9H12,
                         'C9H18': g_C9H18,
                         'BHD': g_BHD,
                         'KERO': g_KERO,
                         'C7H16': g_C7H16}




show_plots = True


# Universal gas constant [J/mol-K]
R0 = 8.31
# Temperature for Formatiom Enthalpy
T0 = 298.15
# Low temperature point in piecewise-linear Cp data - needed to calculate H for correctly fit Terra values
T_base = 100
# temperature delta to fit phase transition latent heat
dT_phase_transition = 200
# first temperature value in Cp range
T_first = 0
# last temperature value in Cp range
T_last = 6000
# temperature delta for Cp range
dT = 100


# dispersed material initialization
gas_material = material.Material(gas_mixture_reference, 'mat_gas', df_props, show_plots, True, T0, T_base, dT_phase_transition, T_first, T_last, dT)




























