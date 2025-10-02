import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import yaml
from matplotlib.patches import Polygon
from decimal import Decimal



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


def output_props(name, is_gas, T_base, df):
    # Make YAML file for QUBIC calculations
    if is_gas:
        file_yaml_out = ''.join(('g_qubiq_', name, '.yaml'))
    else:
        file_yaml_out = ''.join(('p_qubiq_', name, '.yaml'))
    T_header = 'T [K]'

    if not is_gas:
        props = {
            'heat_capacity': 'Cp [Дж/кг-К]',
        }
    else:
        props = {
            'heat_capacity': 'Cp [Дж/кг-К]',
            'viscosity': 'Mu_visc [kg/m-s]',
            'heat_conductivity': 'Lambda [W/m-K]',
            'diffusivity': 'D [m^2/s]'
        }

    print('check data from output_props: name, df:\n', name, df)

    # define dict of dataframe and fill with dataframes
    props_for_yaml = props.copy()
    for key, value in props_for_yaml.items():
        df_for_dict = pd.DataFrame(columns=[T_header, value])
        rarefy_output = False

        if value == 'Cp [Дж/кг-К]':
            # if we want float precision:
            df_for_dict[T_header] = df[T_header].map('{:.1f}'.format).astype(float)
            df_for_dict[value] = df[value].map('{:.1f}'.format).astype(float)
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
            df_for_dict[T_header] = df[T_header].map('{:.3e}'.format).astype(float)
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
        props_for_yaml[key] = df_for_dict
    for key, value in props_for_yaml.items():
        props_for_yaml[key] = list(value.to_dict(orient='dict').values())[0]
    with open(file_yaml_out, 'w') as f:
        yaml.dump(props_for_yaml, f)




    # Make txt file for UDF for FLUENT calculations - write Cp only
    if is_gas:
        file_txt_out = ''.join(('g_fluent_UDF_', name, '.txt'))
    else:
        file_txt_out = ''.join(('p_fluent_UDF_', name, '.txt'))
    with open(file_txt_out, 'w') as f:
        for key, value in props.items():
            if value == 'Cp [Дж/кг-К]':
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
                # f.write(str([float('{:.1f}'.format(x)) for x in strlist]))
                # if integer precision is enough:
                f.write(str([int(x) for x in strlist]))




    if is_gas:
        # Make txt file for Fluent journal to define gas phase properties for FLUENT calculations
        # write Cp, heat_conductivity, viscosity, diffusivity

        print('before Cp journal call')
        str_start_Cp = ''.join(('/define/materials/change-create/', name, ' mixture-template y piecewise-linear '))
        str_end_Cp = ' n n n n\n'
        prop_Cp = get_prop_from_yaml(props_for_yaml['heat_capacity'], str_start_Cp, str_end_Cp, '1')

        str_start_Lambda = ''.join(('/define/materials/change-create/', name, ' mixture-template n y piecewise-linear '))
        str_end_Lambda = ' n n n\n'
        prop_Lambda = get_prop_from_yaml(props_for_yaml['heat_conductivity'], str_start_Lambda, str_end_Lambda, '.000001')

        str_start_viscosity = ''.join(('/define/materials/change-create/', name, ' mixture-template n n y piecewise-linear '))
        str_end_viscosity = ' n n\n'
        prop_viscosity = get_prop_from_yaml(props_for_yaml['viscosity'], str_start_viscosity, str_end_viscosity, '.000001')

        str_start_diffusivity = ''.join(('/define/materials/change-create/', name, ' mixture-template n n n y piecewise-linear '))
        str_end_diffusivity = ' n\n'
        prop_diffusivity = get_prop_from_yaml(props_for_yaml['diffusivity'], str_start_diffusivity, str_end_diffusivity, '.000001')


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
        prop_Cp = get_prop_from_yaml(props_for_yaml['heat_capacity'], str_start_Cp, str_end_Cp, '1')
        with open(''.join(('p_fluent_', name, '.jou')), 'w') as f:
            f.write(';Cp\n')
            f.write(prop_Cp)



    return return_df_for_dict.index.tolist()






