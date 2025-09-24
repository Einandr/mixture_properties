import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *
import time
from datetime import datetime
import os
import csv
import re
import numpy as np
import sqlite3
#import logging
#logging.basicConfig(level=logging.DEBUG)



#work_path = 'D:\\Simulation\\!Python_CIAM\\2_Comb_DES\\'                       # рабочая директория проекта
#f_thermo = 'Ya_Kerosene_Air.dat'

def generate_data_base(dir_thermo, file_thermo, dir_output, file_output):
    '''
    :param dir_thermo: Директория с файлом THERMO CHEMKIN
    :param file_thermo: Имя файла THERMO CHEMKIN
    :param dir_output: Директория где будет находиться выходной файл базы данных (можно оставить ту же что dir_thermo)
    :param file_output: Имя выходного файла базы данных
    :return:
    '''

    time_script_start = datetime.now()
    print('time start:', time_script_start)

    os.chdir(dir_output)
    conn = sqlite3.connect(file_output)
    cursor = conn.cursor()
    try:
        # создание таблицы SQL
        cursor.execute('''CREATE TABLE IF NOT EXISTS thermo (name TEXT, date TEXT, formula TEXT, phase TEXT, t_low REAL, t_mid REAL, t_high REAL,  
                          a1_low REAL, a2_low REAL, a3_low REAL, a4_low REAL, a5_low REAL, a6_low REAL, a7_low REAL,
                          a1_high REAL, a2_high REAL, a3_high REAL, a4_high REAL, a5_high REAL, a6_high REAL, a7_high REAL, atomic TEXT)
                       ''')
        conn.commit()
        cursor.execute('DELETE FROM thermo')
        conn.commit()
    except ValueError:
        print('table already exist')


    def _safe_str_to_Nan_convertion(value):
        return 'Nan' if value == '' else value

    def _parse_formula(symbol_value_string_pair):
        # должен принимать только строку из 5 букв, где буквы 1,2 - символ элемента, буквы 3,4,5 - количество элементов
        symbol = symbol_value_string_pair[:2].strip()
        value = symbol_value_string_pair[2:].strip()
        return symbol, value


    with open(''.join((dir_thermo, '/', file_thermo)), "r") as f:
        # reader = csv.reader(f, delimiter=' ')
        # row_pressure = [row for idx, row in enumerate(reader) if idx == 2]
        temp_range_list = []
        temp_range = []
        a1_low = None
        a2_low = None
        a3_low = None
        a4_low = None
        a5_low = None
        a6_low = None
        a7_low = None
        a1_high = None
        a2_high = None
        a3_high = None
        a4_high = None
        a5_high = None
        a6_high = None
        a7_high = None
        sp_name = None
        sp_date = None
        sp_formula = None
        sp_phase = None
        sp_t_low = None
        sp_t_high = None
        sp_t_mid = None
        atomic = None

        sp_previous_number = 0

        for ind, row in enumerate(f):
            if ind == 1:
                temp_range_list = row.split()
            sp_current_number = (ind - 2) // 4

            if (ind - 2) % 4 == 0 and row.strip() != 'END':
                sp_name = row[:16].strip()
                sp_date_string = ''.join((row[18:20], ' ', row[20:22], ' ', row[22:24]))
                try:
                    sp_date = datetime.strptime(sp_date_string, '%d %m %y').date()
                except ValueError:
                    sp_date = 'None'
                    print('No date information for current species')
                if '&' not in row:
                    print('ROW IS', row)
                    sp_formula_01 = _parse_formula(row[24:29])
                    sp_formula_02 = _parse_formula(row[29:34])
                    sp_formula_03 = _parse_formula(row[34:39])
                    sp_formula_04 = _parse_formula(row[39:44])
                    formula_list = [sp_formula_01, sp_formula_02, sp_formula_03, sp_formula_04]
                    print(formula_list)

                    # sp_formula_01 = row[24:29].strip()
                    # sp_formula_02 = row[29:34].strip()
                    # sp_formula_03 = row[34:39].strip()
                    # sp_formula_04 = row[39:44].strip()
                    sp_formula = ''
                    for formula_tuple in formula_list:
                        symbol, value = formula_tuple
                        print('outside - symbol:', symbol, 'value: ', value)
                        if symbol not in ('00', '0', '') and value not in ('00', '0', ''):
                            print('inside - symbol:', symbol, 'value: ', value)
                            sp_formula += symbol + value

                    # sp_formula = row[24:44].strip()
                else:
                    sp_formula = row.split('&')[1].strip()
                sp_formula = re.sub(r'\s+', '', sp_formula)
                # далее сделано если нет строки выше
                sp_list = (re.sub(r'(\d+[.,]?\d+)', r'\1 ', sp_formula))  # выделяем границу молей предыдущего компонента и имени следующего
                sp_list = re.split(r'\s', sp_list)  # разбиваем строку по пробелам
                sp_list = [x for x in sp_list if x]  # удаляем лишние пробелы
                # sp_components = [x for ind, x in enumerate(sp_list) if ((ind % 2) != 1)]
                # sp_moles = [x for ind, x in enumerate(sp_list) if ((ind % 2) != 0)]
                sp_phase = row[44]
                sp_t_low = float(row[45:55].strip())
                sp_t_high = float(row[55:65].strip())
                sp_t_mid = float(row[65:73].strip())
                atomic = _safe_str_to_Nan_convertion(row[73:78].strip())
            if (ind - 2) % 4 == 1:
                a1_high = float(row[0:15])
                a2_high = float(row[15:30])
                a3_high = float(row[30:45])
                a4_high = float(row[45:60])
                a5_high = float(row[60:75])
            if (ind - 2) % 4 == 2 and ind > 2:
                a6_high = float(row[0:15])
                a7_high = float(row[15:30])
                a1_low = float(row[30:45])
                a2_low = float(row[45:60])
                a3_low = float(row[60:75])
            if (ind - 2) % 4 == 3 and ind > 2:
                a4_low = float(row[0:15])
                a5_low = float(row[15:30])
                a6_low = float(row[30:45])
                a7_low = float(row[45:60])
                try:
                    test_string = ''.join(('INSERT INTO thermo VALUES ( \'', sp_name, '\', \'', str(sp_date), '\', \'', sp_formula, '\', \'', sp_phase, '\', ', str(sp_t_low), ', ', str(sp_t_mid), ', ', str(sp_t_high), ', ',
                                           str(a1_low), ', ', str(a2_low), ', ', str(a3_low), ', ', str(a4_low), ', ', str(a5_low), ', ', str(a6_low), ', ', str(a7_low), ', ',
                                           str(a1_high), ', ', str(a2_high), ', ', str(a3_high), ', ', str(a4_high), ', ', str(a5_high), ', ', str(a6_high), ', ', str(a7_high), ', \'', str(atomic), '\')'))
                    cursor.execute(test_string)
                    conn.commit()
                except ValueError:
                    print('not reached full definition')

            # создание таблицы SQL температурных диапазонов

            try:
                a_all = np.array([[a1_low, a2_low, a3_low, a4_low, a5_low, a6_low, a7_low],
                                  [a1_high, a2_high, a3_high, a4_high, a5_high, a6_high, a7_high]])

            except ValueError:
                print('a not defined yet')

            print(a_all)
            print(row.strip(), len(row), ind)

    for item in temp_range_list:
        temp_range.append(float(item))

    print('\nDatabase SQL table is:\n')

    cursor.execute('select * FROM thermo')

    columns = list(map(lambda x: x[0], cursor.description))
    print(columns, '\n')

    rows = cursor.fetchall()
    for row in rows:
        print(row)

    conn.close()

    time_script_finish = datetime.now()
    print('time finish:', time_script_finish)
    delta = time_script_finish - time_script_start
    print('script execution time:', delta.seconds, 'seconds')
    print('script execution time:', "{:4.2}".format(delta.seconds / 60), 'minutes')
    print('script execution time:', "{:4.2}".format(delta.seconds / 3600), 'hours')
    print('THE END')

    return None
