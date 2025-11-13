import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import yaml
from matplotlib.patches import Polygon
from decimal import Decimal
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap
from collections import OrderedDict


def plot_property(T_grid, property_grid, mixture_name, is_gas, save_path,
                  y_label, title_label, file_name, legend_label=r'$C_p$'):
    """
    Функция для построения и сохранения графика одного свойства смеси.

    :param T_grid: Массив температур.
    :param property_grid: Массив значений свойства.
    :param mixture_name: Имя смеси для названия файла.
    :param is_gas: Флаг, указывающий, является ли смесь газом.
    :param save_path: Путь для сохранения графика.
    :param y_label: Метка оси Y.
    :param title_label: Заголовок графика.
    :param file_name: Имя файла для сохранения.
    :param legend_label: Метка для легенды.
    """
    plt.figure()
    plt.grid(True)
    plt.plot(T_grid, property_grid, label=legend_label)
    plt.legend(fontsize='medium', loc='best')
    plt.title(title_label)
    plt.xlabel(r'$Т\ (К)$')
    plt.ylabel(y_label)
    plt.xlim(0)
    plt.ylim(0)
    plt.savefig(f'{save_path}/{"gas" if is_gas else "disp"}_{mixture_name}_{file_name}.jpeg',
                dpi=400, bbox_inches='tight')
    plt.close()


def plot_mixture_properties(data, mixture_name, is_gas, save_path='.'):
    """
    Функция для построения и сохранения графиков свойств смеси.

    :param data: DataFrame с данными свойств смеси.
    :param mixture_name: Имя смеси для названия файлов.
    :param is_gas: Флаг, указывающий, является ли смесь газом.
    :param save_path: Путь для сохранения графиков.
    """
    T_grid = data['T [K]']
    Cp_grid = data['Cp [J/kg-K]']
    H_grid = data['H [J/kg]']

    # График теплоёмкости
    plot_property(T_grid, Cp_grid, mixture_name, is_gas, save_path,
                  y_label=r'$C_p\ (\frac{Дж}{кг\cdot К})$',
                  title_label=r'$Теплоёмкость$',
                  file_name='Cp',
                  legend_label=r'$C_p$')

    # График энтальпии
    plot_property(T_grid, H_grid, mixture_name, is_gas, save_path,
                  y_label=r'$H\ (\frac{Дж}{кг})$',
                  title_label=r'$Энтальпия$',
                  file_name='H',
                  legend_label=r'$H$')

    if is_gas:
        viscosity_grid = data['Mu_visc [kg/m-s]']
        heat_conductivity_grid = data['Lambda [W/m-K]']
        diffusivity_grid = data['D [m^2/s]']

        # График вязкости
        plot_property(T_grid, viscosity_grid, mixture_name, is_gas, save_path,
                      y_label=r'$\mu\ (Па \cdot с)$',
                      title_label=r'$Вязкость$',
                      file_name='viscosity',
                      legend_label=r'$\mu$')

        # График теплопроводности
        plot_property(T_grid, heat_conductivity_grid, mixture_name, is_gas, save_path,
                      y_label=r'$\lambda\ (\frac{Вт}{м \cdot K})$',
                      title_label=r'$Теплопроводность$',
                      file_name='heat_conductivity',
                      legend_label=r'$\lambda$')

        # График диффузии
        plot_property(T_grid, diffusivity_grid, mixture_name, is_gas, save_path,
                      y_label=r'$D\ (\frac{м^2}{c})$',
                      title_label=r'$Диффузия$',
                      file_name='diffusivity',
                      legend_label=r'$D$')


def get_prop(df, str_start, x, y, quant, str_end):
    print('data from get props: ', df)
    prop = str_start
    prop_points = ''
    n_points = 0
    for ind, item in enumerate(df[x]):
        # print('check index: ', ind, item, y)
        if df[y][item] is not None:
            n_points += 1
            prop_points += ''.join((' ', str(df[x][item]), ' ', str(Decimal(str(df[y][item])).quantize(Decimal(quant)))))
            # print(prop_points)
    prop += str(n_points)
    prop += prop_points
    prop += str_end
    return prop


def get_prop_from_yaml(df, str_start, str_end, quant):
    print('data from get props: ', df)
    prop = str_start
    prop_points = ''
    n_points = 0
    for ind, item in enumerate(df):
        # print('check index: ', ind, item)
        if df[item] is not None:
            n_points += 1
            prop_points += ''.join((' ', str(item), ' ', str(Decimal(str(df[item])).quantize(Decimal(quant)))))
            # print(prop_points)
    prop += str(n_points)
    prop += prop_points
    prop += str_end
    return prop


def output_props(name, is_gas, T_base, T0, df, df_constants):
    print('from output_props T0 and df:', T0, df)


    # Make YAML file for QUBIC calculations
    if is_gas:
        file_yaml_out = ''.join(('g_qubiq_', name, '.yaml'))
    else:
        file_yaml_out = ''.join(('p_qubiq_', name, '.yaml'))

    T_header = 'T [K]'

    # Определяем температурно-зависимые свойства
    if not is_gas:
        props = OrderedDict([
            ('heat_capacity', 'Cp [J/kg-K]')
        ])
    else:
        props = OrderedDict([
            ('heat_capacity', 'Cp [J/kg-K]'),
            ('viscosity', 'Mu_visc [kg/m-s]'),
            ('heat_conductivity', 'Lambda [W/m-K]'),
            ('mass_diffusivity', 'D [m^2/s]')
        ])

    props_for_yaml = CommentedMap()

    # Извлечение значений eps_dk и sigma из df_constants
    eps_dk_value = None
    sigma_value = None
    if 'eps_dk [Kelvins]' in df_constants.columns:
        eps_dk_value = df_constants['eps_dk [Kelvins]'].iloc[0]
    if 'sigma [angstroms]' in df_constants.columns:
        sigma_value = df_constants['sigma [angstroms]'].iloc[0]

    # Создание комментария для первой строки YAML
    comment_header = f"# Transport properties calculated basen on kinetic theory formulas for Lennard-Jones potential well depth eps_dk = {eps_dk_value} Kelvins, Lennard-Jones collision diameter sigma = {sigma_value} angstroms\n"

    # Добавляем константные свойства в начало YAML
    constant_keys = OrderedDict([
            ('formation_heat', 'dH0 [J/kg]'),
            ('molar_mass', 'M [kg/mol]'),
            ('gas_constant', 'R [J/kg-K]'),
            ('density', 'rho [kg/m3]'),
            ('eps_dk', 'eps_dk [Kelvins]'),
            ('sigma', 'sigma [angstroms]')
        ])

    print('df_constants test:', df_constants)

    for key, value in constant_keys.items():
        if value in df_constants.columns:
            if key == 'formation_heat':
                formation_heat_value = df_constants[value].iloc[0]
                props_for_yaml.yaml_set_comment_before_after_key(key, before=f'\n{value}, {T0} K\n')
                formation_heat_map = CommentedMap()
                formation_heat_map[float(T0)] = float(formation_heat_value)
                props_for_yaml[key] = formation_heat_map
            elif key in ['molar_mass', 'gas_constant', 'density']:
                if (is_gas and key in ['molar_mass', 'gas_constant']) or (not is_gas and key in ['molar_mass', 'gas_constant', 'density']):
                    constant_value = df_constants[value].iloc[0]
                    props_for_yaml.yaml_set_comment_before_after_key(key, before=f'\n{value}\n')
                    props_for_yaml[key] = float(constant_value)

    for key, value in props.items():
        df_for_dict = pd.DataFrame(columns=[T_header, value])
        rarefy_output = False

        if value == 'Cp [J/kg-K]':
            # if we want float precision:
            df_for_dict[T_header] = df[T_header].map('{:.2f}'.format).astype(float)
            df_for_dict[value] = df[value].map('{:.2f}'.format).astype(float)
            # if integer precision is enough:
            # df_for_dict[T_header] = df[T_header].astype(int)
            # df_for_dict[value] = df[value].astype(int)
            df_for_dict.set_index(T_header, drop=True, inplace=True)
            # !!! Remove leading repeating values !!!
            df_for_dict.drop(df_for_dict[df_for_dict.index < T_base].index, inplace=True)
            last_value = df_for_dict[value].iloc[-1]
            print('last value: ', last_value)

            # !!! Remove trailing repeating values - if they appear anywhere before, not repeat !!!
            ind_for_drop = []
            for ind, row in df_for_dict[::-1].itertuples():
                # print(ind, row)
                if row == last_value and df_for_dict[value][ind] == df_for_dict.shift()[value][ind]:
                    ind_for_drop.append(ind)
                if row == last_value and df_for_dict[value][ind] != df_for_dict.shift()[value][ind]:
                    print('reached last trailing repeating value, stopped searching to avoid deletion of correctly twice repeating values in midrange')
                    break
            df_for_dict.drop(ind_for_drop, inplace=True)
            print('ind_for_drop:', ind_for_drop)
            if rarefy_output and is_gas:
                ind_for_drop2 = []
                for ind, row in df_for_dict[::-1].itertuples():
                    # rarefy - only for gas - BE CAREFUL NOT TO DELETE PHASE TRANSITION DATA
                    print('check:', is_gas, ind, ind % 200)
                    if ind % 100 != 0:
                        print('here true')
                        ind_for_drop2.append(ind)
                print('ind_for_drop2:', ind_for_drop2)
                df_for_dict.drop(ind_for_drop2, inplace=True)
            print('rarefied Cp df for dict: ', df_for_dict)
            return_df_for_dict = df_for_dict.copy()
        else:
            # if we want float precision:
            df_for_dict[T_header] = df[T_header].map('{:.2f}'.format).astype(float)
            df_for_dict[value] = df[value].map('{:.3e}'.format).astype(float)
            df_for_dict.set_index(T_header, drop=True, inplace=True)
            df_for_dict.drop(df_for_dict[df_for_dict.index < T_base].index, inplace=True)
            # RAREFY heat_capacity viscosity heat_conductivity diffusivity
            ind_for_drop = []
            if rarefy_output:
                for ind, row in df_for_dict[::-1].itertuples():
                    # print(ind, row)
                    if ind % 100 != 0:
                        ind_for_drop.append(ind)
            df_for_dict.drop(ind_for_drop, inplace=True)
            # return_df_for_dict = df_for_dict.copy()
            print('rarefied df for dict: ', df_for_dict)

        data_map = CommentedMap()
        for index, row in df_for_dict.iterrows():
            print('making props_for_yaml index value:', index, value)
            data_map[float(index)] = float(row[value])
            print('data map:', data_map)
        props_for_yaml.yaml_set_comment_before_after_key(key, before=f'\n{value}\n')
        props_for_yaml[key] = data_map
    yaml = YAML()
    yaml.width = 4096  # Увеличиваем ширину для длинных строк
    with open(file_yaml_out, 'w') as f:
        f.write(comment_header)
        yaml.dump(props_for_yaml, f)


    # Make txt file for UDF for FLUENT calculations - write Cp only
    if is_gas:
        file_txt_out = ''.join(('g_fluent_UDF_', name, '.txt'))
    else:
        file_txt_out = ''.join(('p_fluent_UDF_', name, '.txt'))
    with open(file_txt_out, 'w') as f:
        for key, value in props.items():
            if value == 'Cp [J/kg-K]':
                f.write(T_header)
                f.write('\n')
                strlist = return_df_for_dict.index.tolist()
                # if we want float precision:
                # f.write(str([float('{:.1f}'.format(x)) for x in strlist]))
                # if integer precision is enough:
                f.write(str([int(x) for x in strlist]))
                f.write('\n')
                f.write(value)
                f.write('\n')
                print('inside output utils')
                print(key, value)
                strlist = return_df_for_dict[value].values.tolist()
                # if we want float precision:
                f.write(str([float('{:.2f}'.format(x)) for x in strlist]))
                # if integer precision is enough:
                # f.write(str([int(x) for x in strlist]))


    if is_gas:
        # Make txt file for Fluent journal to define gas phase properties for FLUENT calculations
        # write Cp, heat_conductivity, viscosity, diffusivity

        print('before Cp journal call')
        str_start_Cp = ''.join(('/define/materials/change-create/', name, ' mixture-template y piecewise-linear '))
        str_end_Cp = ' n n n n\n'
        prop_Cp = get_prop_from_yaml(props_for_yaml['heat_capacity'], str_start_Cp, str_end_Cp, '.01')

        str_start_Lambda = ''.join(('/define/materials/change-create/', name, ' mixture-template n y piecewise-linear '))
        str_end_Lambda = ' n n n\n'
        prop_Lambda = get_prop_from_yaml(props_for_yaml['heat_conductivity'], str_start_Lambda, str_end_Lambda, '.000001')

        str_start_viscosity = ''.join(('/define/materials/change-create/', name, ' mixture-template n n y piecewise-linear '))
        str_end_viscosity = ' n n\n'
        prop_viscosity = get_prop_from_yaml(props_for_yaml['viscosity'], str_start_viscosity, str_end_viscosity, '.000001')

        str_start_diffusivity = ''.join(('/define/materials/change-create/', name, ' mixture-template n n n y piecewise-linear '))
        str_end_diffusivity = ' n\n'
        prop_diffusivity = get_prop_from_yaml(props_for_yaml['mass_diffusivity'], str_start_diffusivity, str_end_diffusivity, '.000001')


        with open(''.join(('g_fluent_', name, '.jou')), 'w') as f:
            f.write(';Cp\n')
            f.write(prop_Cp)
            f.write(';Lambda\n')
            f.write(prop_Lambda)
            f.write(';Viscosity\n')
            f.write(prop_viscosity)
            f.write(';Diffusivity\n')
            f.write(prop_diffusivity)
    else:
        str_start_Cp = ''.join(('/define/materials/change-create/', name, ' mixture-template n y piecewise-linear '))
        str_end_Cp = ' n n n n n n\n'
        prop_Cp = get_prop_from_yaml(props_for_yaml['heat_capacity'], str_start_Cp, str_end_Cp, '.01')
        with open(''.join(('p_fluent_', name, '.jou')), 'w') as f:
            f.write(';Cp\n')
            f.write(prop_Cp)

    return return_df_for_dict.index.tolist()






