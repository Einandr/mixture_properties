from datetime import datetime
import os
import sqlite3


def generate_db_transport(dir_thransport, file_transport, dir_output, file_output):
    """
    :param dir_thransport: Директория с файлом TRANSPORT CHEMKIN
    :param file_transport: Имя файла TRANSPORT CHEMKIN
    :param dir_output: Директория где будет находиться выходной файл базы данных (можно оставить ту же что dir_thransport)
    :param file_output: Имя выходного файла базы данных
    :return:
    """

    time_start = datetime.now()
    os.chdir(dir_output)
    conn = sqlite3.connect(file_output)
    cursor = conn.cursor()
    try:
        # создание таблицы SQL
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS transport (
                name TEXT,
                geom_index INT,
                L_J_potential_well_depth  REAL,
                L_J_collision_diameter REAL,  
                dipole_momentum REAL,
                polarizability REAL,
                rotational_relaxation_collision_number REAL,
                comment TEXT
            )
        ''')
        conn.commit()
        cursor.execute('DELETE FROM transport')
        conn.commit()
    except ValueError:
        print('table already exist')

    with open(''.join((dir_thransport, '/', file_transport)), "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('!', '/', '.')):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue

            name = parts[0]
            try:
                geom_index = int(parts[1])                                  # If the index is 0, the molecule is a single atom. If the index is 1 the molecule is linear, and if it is 2, the molecule is nonlinear.
                L_J_potential_well_depth = float(parts[2])                  # Kelvins
                L_J_collision_diameter = float(parts[3])                    # angstroms
                dipole_momentum = float(parts[4])                           # Debye -- 10^-18 cm^3/2 ergs^1/2
                polarizability = float(parts[5])                            # cubic angstroms
                rotational_relaxation_collision_number = float(parts[6])    # at 298 K
            except (ValueError, IndexError):
                continue

            comment = None
            if '!' in line:
                comment_pos = line.find('!')
                comment = line[comment_pos + 1:].strip()

            cursor.execute('''
                            INSERT INTO transport (
                                name, geom_index, L_J_potential_well_depth,
                                L_J_collision_diameter, dipole_momentum,
                                polarizability, rotational_relaxation_collision_number,
                                comment
                            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (
                name, geom_index, L_J_potential_well_depth,
                L_J_collision_diameter, dipole_momentum,
                polarizability, rotational_relaxation_collision_number,
                comment
            ))

    conn.commit()

    print('\nDB SQL transport table is:\n')
    cursor.execute('select * FROM transport')
    columns = list(map(lambda x: x[0], cursor.description))
    print(columns, '\n')
    rows = cursor.fetchall()
    for row in rows:
        print(row)
    conn.close()

    time_finish = datetime.now()
    delta = time_finish - time_start
    total_seconds = delta.total_seconds()
    print(f'generate transport DB execution time: {total_seconds:.2f} seconds, {total_seconds / 60:.2f} minutes, {total_seconds / 3600:.2f} hours')
    return None
