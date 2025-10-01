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


def initiate_kinetics(path_db, components_names):
    # Подключение к базе данных
    conn = sqlite3.connect(path_db)
    cursor = conn.cursor()

    # Получаем список доступных в базе имен компонентов
    cursor.execute('SELECT name FROM thermo')
    species_names_db = [row[0] for row in cursor.fetchall()]

    # Проверяем, все ли компоненты из components_names есть в базе
    missing_components = [
        name for name in components_names
        if name not in species_names_db
    ]

    if missing_components:
        raise ValueError(
            f"Следующие компоненты отсутствуют в базе данных: {', '.join(missing_components)}"
        )

    # Если все компоненты найдены, продолжаем
    list_components = []
    for name in components_names:
        cursor.execute('SELECT * FROM thermo WHERE name=?', (name,))
        row = cursor.fetchone()
        if not row:
            raise ValueError(f"Компонент {name} не найден в базе данных")

        # Извлекаем данные из строки
        (
            name, date, formula, phase, atomic,
            t_low, t_mid, t_high,
            a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low,
            a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high
        ) = row

        # Преобразуем числовые значения из строк в float
        t_low = float(t_low)
        t_mid = float(t_mid)
        t_high = float(t_high)

        a_low = [float(a1_low), float(a2_low), float(a3_low), float(a4_low),
                 float(a5_low), float(a6_low), float(a7_low)]
        a_high = [float(a1_high), float(a2_high), float(a3_high), float(a4_high),
                  float(a5_high), float(a6_high), float(a7_high)]

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
                a1_low=a_low[0],
                a2_low=a_low[1],
                a3_low=a_low[2],
                a4_low=a_low[3],
                a5_low=a_low[4],
                a6_low=a_low[5],
                a7_low=a_low[6],
                a1_high=a_high[0],
                a2_high=a_high[1],
                a3_high=a_high[2],
                a4_high=a_high[3],
                a5_high=a_high[4],
                a6_high=a_high[5],
                a7_high=a_high[6],
                dT=200,         # [K] температурная дельта
                T_base=100,     # [K] температура начала таблицы
                T0=298.15       # [K] температура для вычисления энтальпии образования
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
    #     species_names_db = list(map(lambda x: x[0], cursor.fetchall()))
    # except IndexError:
    #     species_names_db = []
    # print('Species names available in DB: ', species_names_db)
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
    #         # self.conn = sqlite3.connect(path_db)
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

    print('Components names in CHEMKIN order: ', components_names)
    return list_components
