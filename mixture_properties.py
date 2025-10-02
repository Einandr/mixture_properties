import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy import interpolate
import math
import ast
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
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Set, Optional

# USER IMPORTS
from utils.generate_data_base import *
from utils import initiate_kinetics as kin
from utils import component_terra as tc
from utils import component_mixture as mc
from utils import material


path = r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS'
dir_run = 'RUN_mixture_properties'
dir_props_out = 'props_out'
dir_solution = 'solution'

path_data = ''.join((os.getcwd(), '\\', 'data'))
terra_props = 'props_TERRA.txt'
chemkin_path = 'ker_mixture'
chemkin_thermo = 'ker_mixture.dat'
chemkin_trans = 'ker_mixture.trans'
chemkin_reactions = 'ker_mixture.inp'
generate_new_chemkin_db = False

file_config = 'run_config.ini'
file_MPL = 'MPL_constants.txt'


path_run = ''.join((path, '/', dir_run))
Path(path_run).mkdir(parents=True, exist_ok=True)
os.chdir(path_run)


df_terra = pd.read_csv(''.join((path_data, '\\', terra_props)), delimiter=' ')
df_terra['T_range'] = df_terra['T_range'].apply(lambda x: np.fromstring(x.replace('\'', '').replace('\"', '').replace('\n', ''), dtype=float, sep=','))
df_terra['f_ranges'] = df_terra['f_ranges'].apply(lambda x: np.array(ast.literal_eval(x.replace('\'', '').replace('\"', '').replace('\n', ''))))

df_terra['name_isomer'] = df_terra['name_isomer'].fillna('')
df_terra['name'] = df_terra['name_brutto'] + ' ' + df_terra['name_isomer']
df_terra['name'] = df_terra['name'].str.strip()
df_terra.set_index('name', inplace=True)



# Создаем файл базы данных в той же директории, где находится файл CHEMKIN THERMO
if generate_new_chemkin_db:
    generate_data_base(''.join((path_data, '\\', chemkin_path)), chemkin_thermo, ''.join((path_data, '\\', chemkin_path)), chemkin_thermo.replace('.dat', '.db'))


g_O2 = {'O2': material.Component(1, material.Source.C)}
g_N2 = {'N2': material.Component(1, material.Source.C)}
g_CO2 = {'CO2': material.Component(1, material.Source.C)}
g_CO = {'CO': material.Component(1, material.Source.T)}
g_H2O = {'H2O': material.Component(1, material.Source.T)}
g_C6H14 = {'C6H14': material.Component(1, material.Source.T)}
g_C10H22 = {'C10H22': material.Component(1, material.Source.T)}
g_C6H6 = {'C6H6': material.Component(1, material.Source.T)}
g_C7H16 = {'C7H16': material.Component(1, material.Source.T)}
g_C9H12 = {'C9H12': material.Component(1, material.Source.T)}
g_C9H18 = {'C9H18': material.Component(1, material.Source.T)}

g_BHD = {'C6H14': material.Component(0.091, material.Source.T),
         'C10H22': material.Component(0.727, material.Source.T),
         'C6H6': material.Component(0.182, material.Source.T)}


g_KERO = {
    'C9H12 unknown': material.Component(0.132, material.Source.T),
    'C10H22 n-Decane': material.Component(0.767, material.Source.C),
    'C9H18 1-Nonene': material.Component(0.101, material.Source.C)
}





# gas_mixture_reference = {'O2': g_O2,
#                          'N2': g_N2,
#                          'CO2': g_CO2,
#                          'CO': g_CO,
#                          'H2O': g_H2O,
#                          'C6H14': g_C6H14,
#                          'C10H22': g_C10H22,
#                          'C6H6': g_C6H6,
#                          'C9H12': g_C9H12,
#                          'C9H18': g_C9H18,
#                          'BHD': g_BHD,
#                          'KERO': g_KERO,
#                          'C7H16': g_C7H16}

gas_mixture_reference = {'KERO': g_KERO}



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


def check_source_conflicts(mixture: Dict) -> Dict[str, Set[material.Source]]:
    component_sources = {}
    for components_dict in mixture.values():
        for component_name, component in components_dict.items():
            if component_name not in component_sources:
                component_sources[component_name] = set()
            component_sources[component_name].add(component.source)
    return component_sources


def get_components_by_source(mixture: Dict, source: material.Source) -> List[str]:
    components = []
    component_sources = check_source_conflicts(mixture)
    for component_name, sources in component_sources.items():
        if source in sources and len(sources) == 1:
            components.append(component_name)
    return components


# Проверка конфликтов
component_sources = check_source_conflicts(gas_mixture_reference)
conflicts = [
    component_name
    for component_name, sources in component_sources.items()
    if len(sources) > 1
]

if conflicts:
    print("Обнаружены конфликты источников для следующих компонентов:")
    for conflict in conflicts:
        print(f"- {conflict}: {component_sources[conflict]}")
    raise ValueError("Обнаружены конфликты источников. Исправьте данные и повторите попытку.")

# Составление списков
components_with_source_C = get_components_by_source(gas_mixture_reference, material.Source.C)
components_with_source_T = get_components_by_source(gas_mixture_reference, material.Source.T)

print("Компоненты с источником CHEMKIN:")
print(components_with_source_C)
print("\nКомпоненты с источником TERRA:")
print(components_with_source_T)






# Инициализация химической кинетики
components_chemkin = kin.initiate_kinetics(''.join((path_data, '\\', chemkin_path, '\\', chemkin_thermo.replace('.dat', '.db'))), components_with_source_C)




# dispersed material initialization
gas_material = material.Material(gas_mixture_reference, 'mat_gas', df_terra, components_chemkin, show_plots, True, T0, T_base, dT_phase_transition, T_first, T_last, dT)




























