import pandas as pd
import fileinput, glob, os
import re
import math
import matplotlib.pyplot as plt
import numpy as np
import pylab
from pathlib import Path
import mendeleev


verify_CHEMKIN = True
verify_TERRA = True
GOST_plots = False
set_x_limit = False

style_C = '-'       # CHEMKIN
style_T = '--'      # TERRA
line_width = 2
x_limit = 3000
terra_Y_norm_factor = 15.7596

verify_sources = [verify_CHEMKIN, verify_TERRA]
number_of_sources = sum(verify_sources)


path = r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS'
dir_run = 'RUN_materials_verification'


# Тут просто добавляем или удаляем необходимые материалы
if verify_CHEMKIN:
    df_chemkin_01 = pd.read_excel(r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS\RUN_chemkin_equilibrium\material_selected_properties.xlsx', index_col='T [K]')
    df_chemkin_02 = pd.read_excel(r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS\RUN_chemkin_equilibrium\material_with_properties.xlsx', index_col='T [K]')


def get_molar_mass(formula):
    # изменением массы в радикалах ионах +- пренебрегается
    elements = re.findall('([A-Z][a-z]*)(\d*)', formula)
    total_mass = 0.0
    for element, count in elements:
        count = int(count) if count else 1  # Если количество не указано, считаем 1
        try:
            atomic_mass = mendeleev.element(element).atomic_weight
            total_mass += atomic_mass * count
            print(f"Ннайден элемент: {element} из формулы {formula} с атомной массой {atomic_mass}")
        except:
            print(f"Не удалось найти элемент: {element} из формулы {formula}")
            return None
    print(f"Для компонента {formula} молярная масса {total_mass / 1000}")
    return total_mass / 1000


# Тут просто добавляем или удаляем необходимые материалы
if verify_TERRA:
    df_terra_01 = pd.read_excel(r'D:\YASIM\VORON\2025_08_KEROSENE_PROPS\terra_results_SI.xlsx', index_col='   T', header=1)
    df_terra_01.columns = df_terra_01.columns.str.strip()

    # Множество столбцов, которые не являются компонентами
    names_not_components = {'m_total', 'p', 'T', 'v', 'S', 'I', 'U', 'M', 'Cp', 'k', 'Cp\'', 'k\'', 'Ap', 'Bv', 'Gt', 'MMg', 'Rg', 'Cpg', 'kg', 'Cp\'g', 'k\'g', 'Mu', 'Lt', 'Lt\'', 'Pr', 'Pr\'', 'A', 'z'}

    # Этап 1: Переименование столбцов (добавление префикса Y_)
    new_columns = []
    for col in df_terra_01.columns:
        if col not in names_not_components:
            new_columns.append(f'Y_{col}')
        else:
            new_columns.append(col)
    df_terra_01.columns = new_columns

    # Этап 2: Расчёт масс компонентов и общей массы
    mass_columns = []
    for col in df_terra_01.columns:
        if col.startswith('Y_'):
            component_name = col[2:]  # Убираем префикс Y_
            M = get_molar_mass(component_name)  # Молярная масса компонента
            if M is not None:
                # print(f'Molar mass for {component_name} is {M:.6f} kg/mol')
                mass_col = f'm_{component_name}'
                df_terra_01[mass_col] = df_terra_01[col] * M  # Масса компонента
                mass_columns.append(mass_col)

    # Суммируем массы всех компонентов по строкам
    df_terra_01['m_total'] = df_terra_01[mass_columns].sum(axis=1)

    # Этап 3: Расчёт массовых долей
    for col in df_terra_01.columns:
        if col.startswith('Y_'):
            component_name = col[2:]
            mass_col = f'm_{component_name}'
            df_terra_01[col] = df_terra_01[mass_col] / df_terra_01['m_total']

    # Удаляем временные столбцы с массами
    df_terra_01 = df_terra_01.drop(columns=mass_columns)

    # Переименование столбцов
    df_terra_01 = df_terra_01.rename(columns={
        'Cp\'': 'Cp [J/kg-K]',
        'I': 'H [J/kg]',
        'MMg': 'M [g/mol]',
        'Rg': 'R [J/kg-K]',
        'Mu': 'Mu_visc [kg/m-s]',
        'Lt': 'Lambda [W/m-K]'
    })

    # Умножаем значения
    df_terra_01['Cp [J/kg-K]'] = df_terra_01['Cp [J/kg-K]'] * 1000
    df_terra_01['H [J/kg]'] = df_terra_01['H [J/kg]'] * 1000


path_run = ''.join((path, '/', dir_run))
Path(path_run).mkdir(parents=True, exist_ok=True)
os.chdir(path_run)

xlabel = r'$T\ (K)$'
y_Cp_label = r'$Cp\ (\frac{Дж}{кг \cdot К})$'
y_H_label = r'$H\ (\frac{Дж}{кг})$'


def plot_result(pic_name, ylabel, tuple_of_data):
    annotation_ratios_list = []                             # привязка линий выносок
    annotation_offsets_x_list = []                          # отступы линий выносок X
    annotation_offsets_y_list = []                          # отступы линий выносок Y
    number_of_plots_from_one_source_list = []               # число графиков из одного источника для верификации

    for ind, item in enumerate(tuple_of_data):
        df = item[0]
        y = item[1]
        style = item[2]
        label = item[3]
        if GOST_plots:
            ylim = item[4]
            annotation_ratios_list.append(item[5])
            annotation_offsets_x_list.append(item[6])
            annotation_offsets_y_list.append(item[7])
            number_of_plots_from_one_source_list.append(len(item[1]))
        if ind == 0:
            ax = df.plot(y=y, use_index=True, style=style, grid=True, label=label, xlabel=xlabel, ylabel=ylabel, xlim=0, markersize=2.5, markevery=200, linewidth=line_width)
        else:
            df.plot(ax=ax, y=y, use_index=True, style=style, grid=True, label=label, xlabel=xlabel, ylabel=ylabel, xlim=0, markersize=2.5, markevery=200, linewidth=line_width)

        if set_x_limit:
            ax.set_xlim(0.0, x_limit)

        # Проверка на 0 и удаление бесполезных линий
        for line in ax.lines:
            if max(line.get_ydata()) == 0:
                line.remove()



    # НАСТРОЙКИ АЛЯ ГОСТ
    if GOST_plots:
        cumulative_lines = sum(number_of_plots_from_one_source_list)
        # Убираем крайние черточки на шкалах
        ax.set_xlim(0, df.index.max())
        # if ind == 1:
        # xticks = ax.get_xticks()
        # xticks = xticks[:-1]
        # ax.set_xticks(xticks)
        # ax.xaxis.set_label_coords(0.98, -0.02)

        ax.set_ylim(ylim[0], ylim[1])
        yticks = ax.get_yticks()
        yticks = yticks[:-1]
        ax.set_yticks(yticks)

        # Временно установить Х лимиты если нужно для какого-то конкретного графика
        if set_x_limit:
            ax.set_xlim(0.0, x_limit)

        # Поворачиваем подпись OY горизонтально
        ax.yaxis.label.set_rotation(0)
        ax.yaxis.set_label_coords(-0.06, 0.95)

        # Стрелочки на осях
        plt.annotate('', xy=(1.02, 0), xytext=(0.0, 0),
             arrowprops=dict(facecolor='black', shrink=0.0, width=0.01, headlength=8, headwidth=5),
             xycoords='axes fraction', textcoords='axes fraction')
        plt.annotate('', xy=(0, 1.02), xytext=(0, 0),
             arrowprops=dict(facecolor='black', shrink=0.0, width=0.01, headlength=8, headwidth=5),
             xycoords='axes fraction', textcoords='axes fraction')

        # Правая и верхняя границы графика
        plt.gca().spines[['top', 'right']].set_color('gray')



        # print('number_of_plots_from_one_source', number_of_plots_from_one_source)
        print(annotation_offsets_x_list)
        print(annotation_offsets_y_list)
        print(number_of_plots_from_one_source_list)

        label_counter = 0
        ind_source = 0
        cumulative_number_of_plots = number_of_plots_from_one_source_list[ind_source]
        for i, line in enumerate(ax.lines):
            if i >= cumulative_number_of_plots:
                ind_source += 1
                cumulative_number_of_plots += number_of_plots_from_one_source_list[ind_source]
            x_data = line.get_xdata()
            y_data = line.get_ydata()


            ind_line_inside_source = i + number_of_plots_from_one_source_list[ind_source] - cumulative_number_of_plots
            # print('i cumulative_number_of_plots number_of_plots_from_one_source_list[ind_source] ind_source ind_line_inside_source', i, cumulative_number_of_plots, number_of_plots_from_one_source_list[ind_source], ind_source, ind_line_inside_source)

            annotation_ratio = annotation_ratios_list[ind_source][ind_line_inside_source]
            annotation_offset_x = annotation_offsets_x_list[ind_source][ind_line_inside_source]
            annotation_offset_y = annotation_offsets_y_list[ind_source][ind_line_inside_source]
            ind = round(annotation_ratio * len(x_data))

            # if (i+1) % 2 != 0:
            # print(line._linestyle, style_S, max_value)
            # if (line._linestyle == style_S) and (cumulative_lines > number_of_sources):
            #     label_counter += 1
            #     label = str(label_counter)
            #     ax.annotate(label, xy=(x_data[ind], y_data[ind]), xytext=(x_data[ind] + annotation_offset_x*max(x_data), y_data[ind] + annotation_offset_y*max(y_data)),
            #                 arrowprops=dict(facecolor='black', shrink=0.00001, headlength=5, headwidth=0.1, width=0.1),
            #                 color='black', fontname='Times New Roman')
            # print('i cumulative_number_of_plots number_of_plots_from_one_source_list[ind_source] ind_source ind_line_inside_source', i, cumulative_number_of_plots, number_of_plots_from_one_source_list[ind_source], ind_source, ind_line_inside_source)

        # легенда
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    if not GOST_plots:
        plt.legend(loc='best', fontsize='small')
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.savefig(''.join(('pic_', pic_name, '.jpeg')), dpi=400, bbox_inches='tight')
    plt.close()


# CHEMKIN
if verify_CHEMKIN and verify_TERRA:
    plot_result('01_Cp', r'$Cp\ (\frac{Дж}{кг \cdot К})$', ((df_chemkin_01, ['Cp [J/kg-K]'], [style_C], ['CHEMKIN - Cp'], [0, 80], [0.9], [-0.3], [-0.1]),
                                                            (df_terra_01, ['Cp [J/kg-K]'], [style_T], ['TERRA - Cp'], [0, 80], [0.9], [-0.3], [-0.1])))
    plot_result('02_H', r'$H\ (\frac{Дж}{кг})$', ((df_chemkin_01, ['H [J/kg]'], [style_C], ['CHEMKIN - H'], [0, 80], [0.9], [-0.3], [-0.1]),
                                                  (df_terra_01, ['H [J/kg]'], [style_T], ['TERRA - H'], [0, 80], [0.9], [-0.3], [-0.1])))
    plot_result('03_mu', r'$\mu (\frac{кг}{м \cdot с})$', ((df_chemkin_01, ['Mu_visc [kg/m-s]'], [style_C], ['CHEMKIN - $\mu$'], [0, 80], [0.9], [-0.3], [-0.1]),
                                                           (df_terra_01, ['Mu_visc [kg/m-s]'], [style_T], ['TERRA - $\mu$'], [0, 80], [0.9], [-0.3], [-0.1])))
    plot_result('04_lambda', r'$\lambda (\frac{Вт}{м \cdot К})$', ((df_chemkin_01, ['Lambda [W/m-K]'], [style_C], ['CHEMKIN - $\lambda$'], [0, 80], [0.9], [-0.3], [-0.1]),
                                                                   (df_terra_01, ['Lambda [W/m-K]'], [style_T], ['TERRA - $\lambda$'], [0, 80], [0.9], [-0.3], [-0.1])))
    plot_result('05_D', r'$D (\frac{м^2}{с})$', ((df_chemkin_01, ['D [m^2/s]'], [style_C], ['CHEMKIN - D'], [0, 80], [0.9], [-0.3], [-0.1]),))
    plot_result('06_M', r'$M (\frac{кг}{моль})$', ((df_chemkin_01, ['M [g/mol]'], [style_C], ['CHEMKIN - M'], [0, 80], [0.9], [-0.3], [-0.1]),
                                                   (df_terra_01, ['M [g/mol]'], [style_T], ['TERRA - M'], [0, 80], [0.9], [-0.3], [-0.1])))
    plot_result('07_R', r'$R (\frac{Дж}{кг \cdot К})$', ((df_chemkin_01, ['R [J/kg-K]'], [style_C], ['CHEMKIN - R'], [0, 80], [0.9], [-0.3], [-0.1]),
                                                         (df_terra_01, ['R [J/kg-K]'], [style_T], ['TERRA - R'], [0, 80], [0.9], [-0.3], [-0.1])))
    plot_result('08_Y', r'$Y$', ((df_chemkin_02, ['Y_H2O', 'Y_CO2', 'Y_CO'], [style_C, style_C, style_C], ['CHEMKIN - массовая доля H2O', 'CHEMKIN - массовая доля CO2', 'CHEMKIN - массовая доля CO'], [0, 1], [0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [-0.1, -0.1, 0.1]),
                                 (df_terra_01, ['Y_H2O', 'Y_CO2', 'Y_CO'], [style_T, style_T, style_T], ['TERRA - массовая доля H2O', 'TERRA - массовая доля CO2', 'TERRA - массовая доля CO'], [0, 1], [0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [-0.1, -0.1, 0.1])))





    # plot_result('01_p_diameter', r'$d\ (мкм)$', ((df_fluent_01, ['p_diameter'], [style_F], ['Fluent D0 G0 - диаметр частицы'], [0, 80], [0.9], [-0.3], [-0.1]),
    #                                              (df_fluent_02, ['p_diameter'], [style_F], ['Fluent D0.333 G0 - диаметр частицы'], [0, 80], [0.8], [-0.3], [-0.1]),
    #                                              (df_fluent_03, ['p_diameter'], [style_F], ['Fluent D1 G0 - диаметр частицы'], [0, 80], [0.7], [0.2], [0.1]),
    #                                              (df_qubiq_particles_01, ['diameter'], [style_Q], ['QUBIQ D0 G0 - диаметр частицы'], [0, 80], [0.4], [0.1], [-0.1]),
    #                                              (df_qubiq_particles_02, ['diameter'], [style_Q], ['QUBIQ D0.333 G0 - диаметр частицы'], [0, 80], [0.4], [0.2], [-0.1]),
    #                                              (df_qubiq_particles_03, ['diameter'], [style_Q], ['QUBIQ D1 G0 - диаметр частицы'], [0, 80], [0.4], [0.3], [-0.1])))
    # plot_result('02_p_mass', r'$m\ (кг)$', ((df_fluent_01, ['p_mass'], [style_F], ['Fluent D0 G0 - масса частицы'], [0, 1.6e-10], [0.6], [0.1], [0.1]),
    #                                         (df_fluent_02, ['p_mass'], [style_F], ['Fluent D0.333 G0 - масса частицы'], [0, 1.6e-10], [0.3], [-0.1], [-0.1]),
    #                                         (df_fluent_03, ['p_mass'], [style_F], ['Fluent D1 G0 - масса частицы'], [0, 1.6e-10], [0.2], [-0.1], [-0.1]),
    #                                         (df_qubiq_particles_01, ['mass'], [style_Q], ['QUBIQ D0 G0 - масса частицы'], [0, 1.6e-10], [0.4], [0.1], [0.1]),
    #                                         (df_qubiq_particles_02, ['mass'], [style_Q], ['QUBIQ D0.333 G0 - масса частицы'], [0, 1.6e-10], [0.4], [0.1], [0.1]),
    #                                         (df_qubiq_particles_03, ['mass'], [style_Q], ['QUBIQ D1 G0 - масса частицы'], [0, 1.6e-10], [0.4], [0.1], [0.1])))
    # plot_result('03_p_mass_all', r'$m\ (кг)$', ((df_script, ['p_mass[kg]', 'p_mass_volatile[kg]', 'p_mass_combustible[kg]', 'p_mass_inert[kg]'], ['.', '-', '-', '-'], ['Скрипт - масса частицы', 'Скрипт - масса летучего компонента частицы', 'Скрипт - масса горючего компонента частицы', 'Скрипт - масса инертного компонента частицы']),
    #                                             (df_fluent, ['p_mass', 'p_mass_volatile', 'p_mass_combustible', 'p_mass_inert'], ['.', '-', '-', '-'], ['Fluent - масса частицы', 'Fluent - масса горючего летучего частицы', 'Fluent - масса горючего компонента частицы', 'Fluent - масса инертного компонента частицы'])))
    # plot_result('03_p_density', r'$\rho\ \left(\frac{кг}{м^3}\right)$', ((df_fluent_01, ['p_density'], [style_F], ['Fluent D0 G0 - плотность частицы'], [550, 690], [0.6], [0.1], [0.02]),
    #                                                                      (df_fluent_02, ['p_density'], [style_F], ['Fluent D0.333 G0 - плотность частицы'], [550, 690], [0.6], [0.1], [0.04]),
    #                                                                      (df_fluent_03, ['p_density'], [style_F], ['Fluent D1 G0 - плотность частицы'], [550, 690], [0.4], [0.1], [0.04]),
    #                                                                      (df_qubiq_particles_01, ['density'], [style_Q], ['QUBIQ D0 G0 - плотность частицы'], [550, 690], [0.8], [0.1], [0.1]),
    #                                                                      (df_qubiq_particles_02, ['density'], [style_Q], ['QUBIQ D0.333 G0 - плотность частицы'], [550, 690], [0.8], [0.1], [0.1]),
    #                                                                      (df_qubiq_particles_03, ['density'], [style_Q], ['QUBIQ D1 G0 - плотность частицы'], [550, 690], [0.8], [0.1], [0.1])))
    # plot_result('04_pT', r'$T\ (K)$', ((df_fluent_01, ['p_temperature'], [style_F], [r'Fluent D0 G0 - температура частицы'], [290, 440], [0.4], [0.1], [-0.1]),
    #                                    (df_fluent_02, ['p_temperature'], [style_F], [r'Fluent D0.333 G0 - температура частицы'], [290, 440], [0.6], [0.1], [-0.1]),
    #                                    (df_fluent_03, ['p_temperature'], [style_F], [r'Fluent D1 G0 - температура частицы'], [290, 440], [0.8], [0.1], [-0.1]),
    #                                    (df_qubiq_particles_01, ['temperature'], [style_Q], ['QUBIQ D0 G0 - температура частицы'], [290, 440], [0.4], [0.1], [-0.1]),
    #                                    (df_qubiq_particles_02, ['temperature'], [style_Q], ['QUBIQ D0.333 G0 - температура частицы'], [290, 440], [0.4], [0.1], [-0.1]),
    #                                    (df_qubiq_particles_03, ['temperature'], [style_Q], ['QUBIQ D1 G0 - температура частицы'], [290, 440], [0.4], [0.1], [-0.1])))
    # plot_result('06_gT', r'$T\ (K)$', ((df_script, ['g_temperature[K]'], ['-'], ['QUBIQ - температура газа']),
    #                                    (df_fluent, ['g_temperature'], ['-'], [r'Fluent - температура газа'])))
    # plot_result('06_pT_gT', r'$T\ (K)$', ((df_script, ['p_temperature[K]', 'g_temperature[K]'], ['-', '-'], ['Скрипт - температура частицы', 'Скрипт - температура газа']),
    #                                       (df_fluent, ['p_temperature', 'g_temperature'], ['-', '-'], ['Fluent - температура частицы', 'Fluent - температура газа'])))
    # plot_result('07_pY', r'$Y$', ((df_script, ['p_Y_volatile', 'p_Y_combustible', 'p_Y_inert'], ['-', '-', '-'], ['Скрипт - массовая доля летучих', 'Скрипт - массовая доля горючих', 'Скрипт - массовая доля инертных']),
    #                               (df_fluent, ['p_Y_volatile', 'p_Y_combustible', 'p_Y_inert'], ['-', '-', '-'], ['Fluent - массовая доля летучих', 'Fluent - массовая доля горючих', 'Fluent - массовая доля инертных'])))
    # plot_result('08_gY', r'$Y$', ((df_script, ['g_Y_oxidizer', 'g_Y_volatile', 'g_Y_nitrogen', 'g_Y_CPCF'], ['o-', 'o-', 'o-', 'o-'], ['Скрипт - массовая доля O2', 'Скрипт - массовая доля летучего', 'Скрипт - массовая доля N2', 'Скрипт - массовая доля конечных продуктов окисления']),
    #                               (df_fluent, ['g_Y_O2', 'g_Y_GV', 'g_Y_N2', 'g_Y_CPCF'], ['-', '-', '-', '-'], ['Fluent - массовая доля O2', 'Fluent - массовая доля летучего', 'Fluent - массовая доля N2', 'Fluent - массовая доля конечных продуктов окисления'])))
    # plot_result('09_Gvol_Gcomb', r'$G\ \left(\frac{kg}{m^2 \cdot c}\right)$', ((df_script, ['G_vol[kg/m^2/s]', 'G_comb[kg/m^2/s]'], ['o-', 'o-'], ['Скрипт - массовый поток летучего', 'Скрипт - массовый поток горючего']),
    #                                                                            (df_fluent, ['G_vol', 'G_comb'], ['-', '-'], ['Fluent - массовый поток летучего', 'Fluent - массовый поток горючего'])))
    # plot_result('10_Theta', r'$\theta$', ((df_script, ['theta_vol', 'theta_comb'], ['-', '-'], ['Скрипт - к-т вдува летучего', 'Скрипт - к-т вдува горючего']),
    #                                       (df_fluent, ['theta_vol', 'theta_comb'], ['-', '-'], ['Fluent - к-т вдува летучего', 'Fluent - к-т вдува горючего'])))
    # plot_result('11_Y_surf', r'$Y_surf$', ((df_script, ['Y_surf_vol', 'Y_surf_comb'], ['o-', 'o-'], ['Скрипт - концентрация летучего GV у поверхности', 'Скрипт - концентрация О2 у поверхности']),
    #                                        (df_fluent, ['Y_surf_vol', 'Y_surf_comb'], ['-', '-'], ['Fluent - концентрация летучего GV у поверхности', 'Fluent - концентрация О2 у поверхности'])))
    # plot_result('12_pQ', r'$Q\ \left(Дж\right)$', ((df_script, ['p_dQ_vol_evap[J]', 'p_dQ_comb_mass[J]', 'p_dQ_comb_prelim[J]', 'p_dQ_comb_final[J]', 'p_dQ_convection[J]'], ['.', '.', '.', '.', '.'], ['Скрипт - p_dQ_vol_evap', 'Скрипт - p_dQ_comb_mass', 'Скрипт - p_dQ_comb_prelim', 'Скрипт - p_dQ_comb_final', 'Скрипт - p_dQ_convection']),
    #                                                (df_fluent, ['p_dQ_vol_evap', 'p_dQ_comb_mass', 'p_dQ_comb_prelim', 'p_dQ_comb_final', 'p_dQ_convection'], ['-', '-', '-', '-', '-'], ['Fluent - p_dQ_vol_evap', 'Fluent - p_dQ_comb_mass', 'Fluent - p_dQ_comb_prelim', 'Fluent - p_dQ_comb_final', 'Fluent - p_dQ_convection'])))
    # plot_result('13_gQ', r'$Q\ \left(Дж\right)$', ((df_script, ['g_dQ[J]', 'g_dQ_GV_mass[J]', 'g_dQ_CPCF_mass[J]', 'g_dQ_oxidizer_mass[J]', 'g_dQ_convection[J]', 'g_dQ_final_combustion[J]', 'g_dQ_from_particle_limits[J]'], ['.', '.', '.', '.', '.', '.', '.'], ['Скрипт - g_dQ', 'Скрипт - g_dQ_GV_mass', 'Скрипт - g_dQ_CPCF_mass', 'Скрипт - g_dQ_oxidizer_mass', 'Скрипт - g_dQ_convection', 'Скрипт - g_dQ_final_combustion', 'Скрипт - g_dQ_from_particle_limits']),
    #                                                (df_fluent, ['g_dQ', 'g_dQ_GV_mass', 'g_dQ_CPCF_mass', 'g_dQ_oxidizer_mass', 'g_dQ_convection', 'g_dQ_final_combustion', 'g_dQ_from_particle_limits'], ['-', '-', '-', '-', '-', '-', '-'], ['Fluent - g_dQ', 'Fluent - g_dQ_GV_mass', 'Fluent - g_dQ_CPCF_mass', 'Fluent - g_dQ_oxidizer_mass', 'Fluent - g_dQ_convection', 'Fluent - g_dQ_final_combustion', 'Fluent - g_dQ_from_particle_limits'])))
    # plot_result('14_pH', r'$H\ \left(\frac{Дж}{кг}\right)$', ((df_script, ['p_enthalpy[J/kg]'], ['-'], ['Скрипт - энтальпия частицы']),
    #                                                           (df_fluent, ['p_enthalpy'], ['-'], ['Fluent - энтальпия частицы'])))
    # plot_result('15_gH', r'$H\ \left(\frac{Дж}{кг}\right)$', ((df_script, ['g_enthalpy[J/kg]'], ['-'], ['Скрипт - энтальпия газа']),
    #                                                           (df_fluent, ['g_enthalpy'], ['-'], ['Fluent - энтальпия газа'])))
    # plot_result('16_dH', r'$dH\ \left(\frac{Дж}{кг}\right)$', ((df_script, ['dH_vol_evap[J/kg]', 'dH_comb_preliminary[J/kg]', 'dH_comb_final[J/kg]'], ['.', '.', '.'], ['Скрипт - тепловой эффект испарения', 'Скрипт - тепловой эффект сгорания первичных', 'Скрипт - тепловой эффект сгорания конечных']),
    #                                                           (df_fluent, ['dH_vol_evap', 'dH_comb_preliminary', 'dH_comb_final'], ['-', '-', '-'], ['Fluent - тепловой эффект испарения', 'Fluent - тепловой эффект сгорания первичных', 'Fluent - тепловой эффект сгорания конечных'])))
    # plot_result('17_k_comb', r'$H\ \left(\frac{с}{м}\right)$', ((df_script, ['k_comb[s/m]'], ['.'], ['Скрипт - k_comb']),
    #                                                             (df_fluent, ['k_comb'], ['-'], ['Fluent - k_comb'])))
    # plot_result('18_Y_comb', r'$Y$', ((df_script, ['Y_surf_comb', 'g_Y_oxidizer'], ['.', '.'], ['Скрипт - массовая доля O2 на поверхности', 'Скрипт - массовая доля O2']),
    #                                   (df_fluent, ['Y_surf_comb', 'g_Y_O2'], ['-', '-'], ['Fluent - массовая доля O2 на поверхности', 'Fluent - массовая доля O2'])))



#
# # SCRIPT - QUBIQ VERIFICATION
# if verify_SCRIPT and verify_TERRA and not verify_CHEMKIN:
#     plot_result('01_p_diameter', r'$d,\ мкм$', ((df_script, ['p_diameter[mkm]'], [style_S], ['Скрипт - диаметр частицы'], [0, 20], [0.2], [0.1], [0.1]),
#                                                  (df_qubiq_particles, ['diameter'], [style_Q], ['QUBIQ - диаметр частицы'], [0, 20], [0.4], [0.1], [0.1])))
#     plot_result('04_pT', r'$T,\ K$', ((df_script, ['p_temperature[K]'], [style_S], ['Скрипт - температура частицы'], [1800, 2000], [0.2], [0.1], [-0.1]),
#                                        (df_qubiq_particles, ['temperature'], [style_Q], [r'QUBIQ - температура частицы'], [1800, 2000], [0.4], [0.1], [-0.1])))
#     plot_result('07_pY', r'$Y$', ((df_script, ['p_Y_volatile', 'p_Y_combustible', 'p_Y_inert'], [style_S, style_S, style_S], ['Скрипт - массовая доля летучих', 'Скрипт - массовая доля горючих', 'Скрипт - массовая доля инертных'], [0, 1], [0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [-0.1, -0.1, 0.1]),
#                                   (df_qubiq_particles, ['Y1', 'Y2', 'Y3'], [style_Q, style_Q, style_Q], ['QUBIQ - массовая доля летучих', 'QUBIQ - массовая доля горючих', 'QUBIQ - массовая доля инертных'], [0, 1], [0.15, 0.35, 0.55], [0.1, 0.1, 0.1], [0.1, 0.1, -0.1])))
#     plot_result('02_p_mass', r'$m,\ кг$', ((df_script, ['p_mass[kg]'], [style_S], ['Скрипт - масса частицы'], [0, 1e-11], [0.2], [0.1], [0.1]),
#                                            (df_qubiq_particles, ['mass'], [style_Q], ['QUBIQ - масса частицы'], [0, 1e-11], [0.4], [0.1], [0.1])))
#     plot_result('03_p_mass_all', r'$m,\ кг$', ((df_script, ['p_mass_volatile[kg]', 'p_mass_combustible[kg]', 'p_mass_inert[kg]', 'p_mass[kg]'], [style_S, style_S, style_S, style_S], ['Скрипт - масса летучего компонента частицы', 'Скрипт - масса горючего компонента частицы', 'Скрипт - масса инертного компонента частицы', 'Скрипт - масса частицы'], [0, 1e-11], [0.1, 0.1, 0.05, 0.1], [0.1, 0.1, 0.1, 0.1], [0.2, 0.2, 0.1, 0.2]),
#                                                (df_qubiq_particles, ['mass_volatile', 'mass_combustible', 'mass_inert', 'mass'], [style_Q, style_Q, style_Q, style_Q], ['QUBIQ - масса горючего летучего частицы', 'QUBIQ - масса горючего компонента частицы', 'QUBIQ - масса инертного компонента частицы', 'QUBIQ - масса частицы'], [0, 1e-11], [0.15, 0.15, 0.1, 0.15], [0.1, 0.1, 0.1, 0.1], [0.2, 0.2, 0.1, 0.2])))
#
# if verify_TERRA_gas:
#         plot_result('05_gT', r'$T\ (K)$', ((df_script, ['g_temperature[K]'], ['-'], ['Скрипт - температура газа'], [1800, 2000], [0.2], [0.1], [0.1]),
#                                            (df_qubiq_gas, ['temperature'], [style_Q], ['QUBIQ - температура газа'], [1800, 2000], [0.4], [0.1], [0.1])))
#         plot_result('06_pT_gT', r'$T\ (K)$', ((df_script, ['p_temperature[K]', 'g_temperature[K]'], [style_S, style_S], ['Скрипт - температура частицы', 'Скрипт - температура газа'], [1950, 1975], [0.2, 0.4], [0.1, 0.1], [-0.05, 0.01]),
#                                               (df_qubiq_particles, ['temperature'], [style_Q], [r'QUBIQ - температура частицы'], [1950, 1975], [0.1], [0.1], [-0.05]),
#                                               (df_qubiq_gas, ['temperature'], [style_Q], ['QUBIQ - температура газа'], [1950, 1975], [0.3], [0.1], [0.1])))
#         plot_result('08_gY', r'$Y$', ((df_script, ['g_Y_oxidizer', 'g_Y_volatile', 'g_Y_nitrogen', 'g_Y_CPCF'], [style_S, style_S, style_S, style_S], ['Скрипт - массовая доля O2', 'Скрипт - массовая доля летучего', 'Скрипт - массовая доля N2', 'Скрипт - массовая доля конечных продуктов окисления'], [0, 1], [0.1, 0.3, 0.5, 0.7], [0.1, 0.1, 0.1, 0.1], [0.2, 2, -0.2, 2]),
#                                       (df_qubiq_gas, ['mass_fractions_O2', 'mass_fractions_GV', 'mass_fractions_N2', 'mass_fractions_CPCF'], [style_Q, style_Q, style_Q, style_Q], ['QUBIQ - массовая доля O2', 'QUBIQ - массовая доля летучего', 'QUBIQ - массовая доля N2', 'QUBIQ - массовая доля конечных продуктов окисления'], [0, 1], [0.15, 0.35, 0.55, 0.75], [0.1, 0.1, 0.1, 0.1], [0.2, 30, -0.2, 2])))
#         plot_result('06_gP', r'$P\ (Па)$', ((df_script, ['g_pressure[Pa]'], [style_S], ['Скрипт - давление газа'], [499500, 500500], [0.2], [0.1], [0.0005]),
#                                             (df_qubiq_gas, ['pressure'], [style_Q], ['QUBIQ - давление газа'], [499500, 500500], [0.3], [0.1], [0.0005])))
#
#


print('debug')







