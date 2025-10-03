import yaml
from pathlib import Path
import os
import numpy as np

# USER IMPORTS
from utils.generate_db_thermo import *
from utils import initiate_kinetics as kin




path_data = ''.join((os.getcwd(), '\\', 'data', '\\', 'ker_mixture'))

chemkin_thermo = 'ker_mixture.dat'
chemkin_reactions = 'ker_air_china_rp3_19sp.inp'

dir_output = r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS'
file_db_output = chemkin_thermo.replace('.dat', '.db')





dir_run = 'RUN_thermo_converter'
generate_new_db_thermo = True


# Создаем файл базы данных в той же директории, где находится файл CHEMKIN THERMO
if generate_new_db_thermo:
    generate_db_thermo(path_data, chemkin_thermo, path_data, file_db_output)

# Инициализация химической кинетики
get_species_names = kin.read_chemkin_file(''.join((path_data, '/', chemkin_reactions)))
list_components = kin.initiate_kinetics(''.join((path_data, '/', file_db_output)), get_species_names)

# test = list_components[0].formula

path_run = ''.join((dir_output, '/', dir_run))
Path(path_run).mkdir(parents=True, exist_ok=True)
os.chdir(path_run)

def numpy_to_list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.int64, np.float64)):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


yaml.add_representer(np.ndarray, lambda dumper, data: dumper.represent_list(numpy_to_list(data)))
yaml.add_representer(np.int64, lambda dumper, data: dumper.represent_int(numpy_to_list(data)))
yaml.add_representer(np.float64, lambda dumper, data: dumper.represent_float(numpy_to_list(data)))

for TC in list_components:
    T_grid_list = numpy_to_list(TC.T_grid)
    Cp_grid_list = numpy_to_list(TC.Cp_grid)
    heat_capacity_dict = {T: Cp for T, Cp in zip(T_grid_list, Cp_grid_list)}
    props = {'formation_heat': {TC.T0: TC.H0_298_J()},
             'gas_constant': TC.R,
             'heat_capacity': heat_capacity_dict}
    with open(f'props_{TC.name}.yaml', 'w') as f:
        yaml.dump(props, f)

print('THE END')


