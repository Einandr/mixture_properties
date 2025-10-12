import pandas as pd
import os
from pathlib import Path
import ast

# USER IMPORTS
from utils.generate_db_thermo import *
from utils.generate_db_transport import *
from utils import initiate_kinetics as kin
from utils import material
from utils import output as ou

path = r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS'

# Пути к результатам расчета в CHEMKIN EQUILIBRIUM
dir_chemkin_equilibrium_solution = 'CHEMKIN'
file_chemkin_equilibrium_solution = 'Dagaut_original.xlsx'
dir_run = 'RUN_chemkin_equilibrium'

# Пути к кинетике
chemkin_path = r'D:\YASIM\!Chemical_Kinetics\ker_Dagaut'
chemkin_thermo = 'thermDagaut.dat'
chemkin_transport = 'transpDagaut.dat'
chemkin_reactions = 'Dagaut_Ori.inp.txt'
generate_new_db_thermo = False
generate_new_db_transport = False

path_data = ''.join((os.getcwd(), '\\', 'data'))
terra_props = 'props_TERRA.txt'

path_run = ''.join((path, '/', dir_run))
Path(path_run).mkdir(parents=True, exist_ok=True)
os.chdir(path_run)

# НЕ НУЖНО ТАК КАК ВСЕ КОМПОНЕНТЫ ЧЕМКИН
df_terra = pd.read_csv(''.join((path_data, '\\', terra_props)), delimiter=' ')
df_terra['T_range'] = df_terra['T_range'].apply(lambda x: np.fromstring(x.replace('\'', '').replace('\"', '').replace('\n', ''), dtype=float, sep=','))
df_terra['f_ranges'] = df_terra['f_ranges'].apply(lambda x: np.array(ast.literal_eval(x.replace('\'', '').replace('\"', '').replace('\n', ''))))
df_terra['name_isomer'] = df_terra['name_isomer'].fillna('')
df_terra['name'] = df_terra['name_brutto'] + ' ' + df_terra['name_isomer']
df_terra['name'] = df_terra['name'].str.strip()
df_terra.set_index('name', inplace=True)


df = pd.read_excel(''.join((path, '/', dir_chemkin_equilibrium_solution, '/', file_chemkin_equilibrium_solution)), sheet_name='1.soln_vs_parameter')
df.columns = df.columns.str.strip()
df.columns = df.columns.str.replace('Equilibrium_Mass_fraction_', 'Y_', regex=False)
df.columns = df.columns.str.replace('Initial_Mass_fraction_', 'Y_in_', regex=False)
df.columns = df.columns.str.replace('_\\(\\)', '', regex=True)
df = df.sort_values(by='Equilibrium_Temperature_(K)').reset_index(drop=True)
df['Cp [J/kg-K]'] = df['Equilibrium_Enthalpy_(J/kg)'].diff() / df['Equilibrium_Temperature_(K)'].diff()
df.loc[0, 'Cp [J/kg-K]'] = df.loc[1, 'Cp [J/kg-K]']


equilibrium_columns = [col for col in df.columns if col.startswith('Y_') and not col.startswith('Y_in_')]
gas_mixture_reference = {}
for component in equilibrium_columns:
    component_name = component[2:]  # Убираем префикс 'Y_'
    var_name = f'g_{component_name}'
    component_dict = {component_name: material.Component(1, material.Source.C)}
    gas_mixture_reference[component_name] = component_dict

# Компоненты в равновесной смеси:
print('Компоненты в равновесной смеси:\n', equilibrium_columns)


# Создаем файлы БД в той же директории, где находится файлы CHEMKIN THERMO, CHEMKIN TRANSPORT
if generate_new_db_thermo:
    generate_db_thermo(chemkin_path, chemkin_thermo, chemkin_path, chemkin_thermo.replace('.dat', '.db'))
if generate_new_db_transport:
    generate_db_transport(chemkin_path, chemkin_transport, chemkin_path, chemkin_transport.replace('.dat', '.db'))


show_plots = True


# Universal gas constant [J/mol-K]
R0 = 8.31
# Temperature for Formation Enthalpy
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


# Проверка конфликтов
component_sources = material.check_source_conflicts(gas_mixture_reference)
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
components_with_source_C = material.get_components_by_source(gas_mixture_reference, material.Source.C)


# Инициализация химической кинетики
os.chdir(path_run)
components_chemkin = kin.initiate_kinetics(''.join((chemkin_path, '\\', chemkin_thermo.replace('.dat', '.db'))), components_with_source_C, show_plots=show_plots, T_last=T_last, path_db_transport=''.join((chemkin_path, '\\', chemkin_transport.replace('.dat', '.db'))))

# material initialization
gas_material = material.Material(gas_mixture_reference, 'mat_gas', df_terra, components_chemkin, show_plots, True, T0, T_base, dT_phase_transition, T_first, T_last, dT)


# вычисляем новые данные
df = df.rename(columns={'Equilibrium_Temperature_(K)': 'T [K]', 'Equilibrium_Enthalpy_(J/kg)': 'H [J/kg]'})
df['Mu_visc [kg/m-s]'] = 0.0
df['Lambda [W/m-K]'] = 0.0
df['D [m^2/s]'] = 0.0
df['M [kg/mol]'] = 0.0
df['R [J/kg-K]'] = 0.0
df['eps_dk [Kelvins]'] = 0.0
df['sigma [angstroms]'] = 0.0


for index, row in df.iterrows():
    # Создаём словарь с массовыми долями компонентов для текущей строки
    mass_fractions_dict = {col[2:]: row[col] for col in df.columns if col.startswith('Y_') and not col.startswith('Y_in')}
    T = row['T [K]']
    df.loc[index, 'Mu_visc [kg/m-s]'] = gas_material.get_mixture_viscosity(mass_fractions_dict, T)
    df.loc[index, 'Lambda [W/m-K]'] = gas_material.get_mixture_heat_conductivity(mass_fractions_dict, T)
    df.loc[index, 'D [m^2/s]'] = gas_material.get_mixture_diffusivity(mass_fractions_dict, T)
    df.loc[index, 'M [g/mol]'] = gas_material.get_mixture_mu(mass_fractions_dict)
    df.loc[index, 'R [J/kg-K]'] = gas_material.get_mixture_R(mass_fractions_dict)
    df.loc[index, 'eps_dk [Kelvins]'] = gas_material.get_mixture_eps_dk(mass_fractions_dict)
    df.loc[index, 'sigma [angstroms]'] = gas_material.get_mixture_sigma(mass_fractions_dict)



df.to_excel('material_with_properties.xlsx', index=False)

# Определяем верхнее ограничение по температуре (по умолчанию None — не использовать)
max_temperature = None

if max_temperature is not None:
    filtered_df = df[df['T [K]'] <= max_temperature]
else:
    filtered_df = df

average_values = {
    'M [g/mol]': filtered_df['M [g/mol]'].mean(),
    'R [J/kg-K]': filtered_df['R [J/kg-K]'].mean(),
    'eps_dk [Kelvins]': filtered_df['eps_dk [Kelvins]'].mean(),
    'sigma [angstroms]': filtered_df['sigma [angstroms]'].mean(),
}

dH0_value = df.loc[df['T [K]'] == T0, 'H [J/kg]'].values
if len(dH0_value) > 0:
    average_values['dH0 [J/kg]'] = dH0_value[0]
else:
    average_values['dH0 [J/kg]'] = None

df_average_properties = pd.DataFrame({
    'M [g/mol]': [average_values['M [g/mol]']],
    'R [J/kg-K]': [average_values['R [J/kg-K]']],
    'dH0 [J/kg]': [average_values['dH0 [J/kg]']],
    'eps_dk [Kelvins]': [average_values['eps_dk [Kelvins]']],
    'sigma [angstroms]': [average_values['sigma [angstroms]']],
})

df_average_properties.to_excel('material_average_properties.xlsx', index=False)
selected_df = df[['T [K]', 'Cp [J/kg-K]', 'H [J/kg]', 'Mu_visc [kg/m-s]', 'Lambda [W/m-K]', 'D [m^2/s]', 'M [g/mol]', 'R [J/kg-K]']]
selected_df.to_excel('material_selected_properties.xlsx', index=False)

path_output_material = 'material'
os.makedirs(path_output_material, exist_ok=True)
os.chdir(path_output_material)

ou.output_props('material_gas', True, T_base, T0, selected_df, df_average_properties)
os.chdir('..')


print('debug')