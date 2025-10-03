import sqlite3
import csv
from . component_chemkin import *



def read_chemkin_file(path):
    species_names_strings = []
    write_species = False
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            index_comment = line.find('!')
            if index_comment != -1:
                line = line[:index_comment]
            print(line)
            if line == 'SPECIES':
                write_species = True
            if line == 'END':
                write_species = False
            if write_species:
                species_names_strings.append(line.strip())
    species_names_strings.remove('SPECIES')
    species_names_strings = [x for x in species_names_strings if x]
    species_names = []
    for item in species_names_strings:
        species_names.extend(item.split())
    print('Species order form Chemkin file: ', species_names, '\n')
    with open('species_list.csv', 'w', newline='') as f:
        cw = csv.writer(f, delimiter=' ')
        cw.writerow(species_names)
    return species_names


def initiate_kinetics(path_db_thermo, components_names, T_last=None, path_db_transport=None):
    # Подключение к базе данных термодинамики
    conn_thermo = sqlite3.connect(path_db_thermo)
    cursor_thermo = conn_thermo.cursor()

    # Получаем список доступных в базе имен компонентов
    cursor_thermo.execute('SELECT name FROM thermo')
    names = [row[0] for row in cursor_thermo.fetchall()]

    # Проверяем, все ли компоненты из components_names есть в базе
    missing_components = [
        name for name in components_names
        if name not in names
    ]

    if missing_components:
        raise ValueError(
            f"Следующие компоненты отсутствуют в БД THERMO: {', '.join(missing_components)}"
        )

    # Подключение к базе данных транспортных свойств, если путь передан
    conn_transport = None
    cursor_transport = None
    if path_db_transport is not None:
        conn_transport = sqlite3.connect(path_db_transport)
        cursor_transport = conn_transport.cursor()


    # Если все компоненты найдены, продолжаем
    list_components = []
    for name in components_names:
        cursor_thermo.execute('SELECT * FROM thermo WHERE name=?', (name,))
        row = cursor_thermo.fetchone()
        if not row:
            raise ValueError(f"Компонент {name} не найден в базе данных")

        # Получаем описание столбцов
        columns = [desc[0] for desc in cursor_thermo.description]

        # Создаем словарь для сопоставления имен столбцов с их индексами
        column_indices = {column: index for index, column in enumerate(columns)}

        # Извлекаем данные из строки по именам столбцов
        name = row[column_indices['name']]
        date = row[column_indices['date']]
        formula = row[column_indices['formula']]
        phase = row[column_indices['phase']]
        t_low = float(row[column_indices['t_low']])
        t_mid = float(row[column_indices['t_mid']])
        t_high = float(row[column_indices['t_high']])
        a1_low = float(row[column_indices['a1_low']])
        a2_low = float(row[column_indices['a2_low']])
        a3_low = float(row[column_indices['a3_low']])
        a4_low = float(row[column_indices['a4_low']])
        a5_low = float(row[column_indices['a5_low']])
        a6_low = float(row[column_indices['a6_low']])
        a7_low = float(row[column_indices['a7_low']])
        a1_high = float(row[column_indices['a1_high']])
        a2_high = float(row[column_indices['a2_high']])
        a3_high = float(row[column_indices['a3_high']])
        a4_high = float(row[column_indices['a4_high']])
        a5_high = float(row[column_indices['a5_high']])
        a6_high = float(row[column_indices['a6_high']])
        a7_high = float(row[column_indices['a7_high']])
        atomic = row[column_indices['atomic']]

        print('check kinetics initiate: name, date, formula, phase, atomic, t_low, t_mid, t_high, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high',
              name, date, formula, phase, atomic, t_low, t_mid, t_high, a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low, a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high)

        # Извлечение данных из базы данных транспортных свойств, если она подключена
        geom_index = None                               # not used yet
        L_J_potential_well_depth = 100                  # default value
        L_J_collision_diameter = 4                      # default value
        dipole_momentum = None                          # not used yet
        polarizability = None                           # not used yet
        rotational_relaxation_collision_number = None   # not used yet
        comment = None                                  # not used yet

        if conn_transport is not None:
            cursor_transport.execute('SELECT * FROM transport WHERE name=?', (name,))
            transport_row = cursor_transport.fetchone()
            if transport_row:
                transport_columns = [desc[0] for desc in cursor_transport.description]
                transport_column_indices = {column: index for index, column in enumerate(transport_columns)}

                geom_index = transport_row[transport_column_indices['geom_index']]
                L_J_potential_well_depth = transport_row[transport_column_indices['L_J_potential_well_depth']]
                L_J_collision_diameter = transport_row[transport_column_indices['L_J_collision_diameter']]
                dipole_momentum = transport_row[transport_column_indices['dipole_momentum']]
                polarizability = transport_row[transport_column_indices['polarizability']]
                rotational_relaxation_collision_number = transport_row[transport_column_indices['rotational_relaxation_collision_number']]
                comment = transport_row[transport_column_indices['comment']]

        # Создаем объект Component
        list_components.append(
            Component(
                name=name,
                date=date,
                formula=formula,
                phase=phase,
                T_low=t_low,
                T_mid=t_mid,
                T_high=t_high,
                atomic=atomic,
                a1_low=a1_low,
                a2_low=a2_low,
                a3_low=a3_low,
                a4_low=a4_low,
                a5_low=a5_low,
                a6_low=a6_low,
                a7_low=a7_low,
                a1_high=a1_high,
                a2_high=a2_high,
                a3_high=a3_high,
                a4_high=a4_high,
                a5_high=a5_high,
                a6_high=a6_high,
                a7_high=a7_high,
                dT=200,         # [K] температурная дельта
                T_base=100,     # [K] температура начала таблицы
                T0=298.15,       # [K] температура для вычисления энтальпии образования
                eps_dk=L_J_potential_well_depth,
                sigma=L_J_collision_diameter,
                T_last=T_last
            )
        )











    # # вывод полной информации об используемой химической кинетики из базы
    # cursor.execute('select * FROM thermo')
    # columns = list(map(lambda x: x[0], cursor.description))
    # print('Database column names:')
    # print(columns, '\n')
    # rows = cursor.fetchall()
    # for row in rows:
    #     print(row)
    #
    # # список доступных в базе имен компонентов
    # cursor.execute('select name FROM thermo')
    # try:
    #     names = list(map(lambda x: x[0], cursor.fetchall()))
    # except IndexError:
    #     names = []
    # print('Species names available in DB: ', names)
    #
    # # список доступных в базе химических формул компонентов
    # cursor.execute('select formula FROM thermo')
    # try:
    #     species_formulas_db = list(map(lambda x: x[0], cursor.fetchall()))
    # except IndexError:
    #     species_formulas_db = []
    # print('Species formulas available in DB: ', species_formulas_db)
    #
    # class GetComponentFromDB:
    #     def __init__(self, name):
    #         self.name = name
    #         self.query_select_formula = ''.join(('SELECT formula FROM thermo WHERE name=\'', name, '\''))
    #         self.query_select_trange = ''.join(('SELECT t_low, t_mid, t_high FROM thermo WHERE name=\'', name, '\''))
    #         self.query_select_a_low = ''.join(('SELECT a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low FROM thermo WHERE name=\'', name, '\''))
    #         self.query_select_a_high = ''.join(('SELECT a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high FROM thermo WHERE name=\'', name, '\''))
    #         # self.conn = sqlite3.connect(path_db_thermo)
    #         # cursor = conn.cursor()
    #
    #         # date
    #         cursor.execute(''.join(('SELECT date FROM thermo WHERE name=\'', self.name, '\'')))
    #         try:
    #             self.date = cursor.fetchall()[0][0]
    #         except IndexError:
    #             self.date = []
    #
    #         # Formula
    #         cursor.execute(''.join(('SELECT formula FROM thermo WHERE name=\'', self.name, '\'')))
    #         try:
    #             self.formula = cursor.fetchall()[0][0]
    #         except IndexError:
    #             self.formula = []
    #
    #         # Phase
    #         cursor.execute(''.join(('SELECT phase FROM thermo WHERE name=\'', self.name, '\'')))
    #         try:
    #             self.phase = cursor.fetchall()[0][0]
    #         except IndexError:
    #             self.phase = []
    #
    #         # Phase
    #         cursor.execute(''.join(('SELECT atomic FROM thermo WHERE name=\'', self.name, '\'')))
    #         try:
    #             self.atomic = cursor.fetchall()[0][0]
    #         except IndexError:
    #             self.atomic = []
    #
    #         # T_range
    #         cursor.execute(''.join(('SELECT t_low, t_mid, t_high FROM thermo WHERE name=\'', self.name, '\'')))
    #         try:
    #             self.T_range = list(cursor.fetchall()[0])
    #         except IndexError:
    #             self.T_range = []
    #
    #         # a_low
    #         cursor.execute(''.join(('SELECT a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low FROM thermo WHERE name=\'', self.name, '\'')))
    #         try:
    #             self.a_low = list(cursor.fetchall()[0])
    #         except IndexError:
    #             self.a_low = []
    #
    #         # a_high
    #         cursor.execute(''.join(('SELECT a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high FROM thermo WHERE name=\'', self.name, '\'')))
    #         try:
    #             self.a_high = list(cursor.fetchall()[0])
    #         except IndexError:
    #             self.a_high = []


    # dT = 200
    # T_base = 100
    # T0 = 298.15  # [K] температура для вычисления энтальпии образования

    # DEBUG TEST BEGIN
    # test_name = 'N2'
    # # test = get_species_data_from_db(test_name)
    # TC = GetComponentFromDB(test_name) # Test Component
    # print('name and T_range:', TC.name, TC.T_range)

    # TC2 = Component(TC.name, TC.date, TC.formula, TC.phase, TC.T_range[0], TC.T_range[1], TC.T_range[2], TC.atomic, TC.a_low[0], TC.a_low[1], TC.a_low[2], TC.a_low[3], TC.a_low[4], TC.a_low[5], TC.a_low[6], TC.a_high[0], TC.a_high[1], TC.a_high[2], TC.a_high[3], TC.a_high[4], TC.a_high[5], TC.a_high[6], dT, T_base, T0)
    # print('TC2 T_range: ', TC2.T_range)
    # print('TC2 a_low: ', TC2.a[0])
    # print('TC2 a_high: ', TC2.a[1])
    # print('TC2 T_grid: ', TC2.T_grid)
    # print('TC2 Cp_grid: ', TC2.Cp_grid)
    # print('TC2 T0: ', TC2.T0)
    # print('TC2 H0_298: ', TC2.H0_298_J())
    # print('TC2 mu: ', TC2.mu)
    # print('TC2 R: ', TC2.R)

    # print('Cp low region: ', TC2.cp(100))
    # print('Cp high region: ', TC2.cp(6000))
    # print('H0 298: ', TC2.h0_298())
    # DEBUG TEST END

    # Создаем список экземпляров класса Component, соблюдая порядок компонент как в CHEMKIN!!!:
    # list_components = []
    # for name in components_names:
    #     TC = GetComponentFromDB(name)           # Temporary Component - объект со свойствами из базы данных
    #     list_components.append(Component(TC.name, TC.date, TC.formula, TC.phase, TC.T_range[0], TC.T_range[1], TC.T_range[2], TC.atomic, TC.a_low[0], TC.a_low[1], TC.a_low[2], TC.a_low[3], TC.a_low[4], TC.a_low[5], TC.a_low[6], TC.a_high[0], TC.a_high[1], TC.a_high[2], TC.a_high[3], TC.a_high[4], TC.a_high[5], TC.a_high[6], dT, T_base, T0))


    conn_thermo.close()
    if conn_transport is not None:
        conn_transport.close()

    print('Components names in CHEMKIN order: ', components_names)
    return list_components
