# your_script.py
import pandas as pd
import psycopg2
from psycopg2 import Error
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import py3Dmol
import json
from django.conf import settings
connection_params = {
        'user': settings.DATABASES['default']['USER'],
        'password': settings.DATABASES['default']['PASSWORD'],
        'host': settings.DATABASES['default']['HOST'],
        'port': settings.DATABASES['default']['PORT'],
        'database': settings.DATABASES['default']['NAME'],
    }

def create_molecule_image(connection_params, value):
    connection = None
    atom_numbers = []

    try:
        # Создаем новое соединение на основе параметров
        connection = psycopg2.connect(**connection_params)
        cursor = connection.cursor()

        # Получаем данные для молекулы из базы данных
        cursor.execute(f'select * from connections where id = {value}')
        columns = [col.name for col in cursor.description]
        connections_current = pd.DataFrame(cursor.fetchall(), columns=columns)

        cursor.execute(f'select * from coordinates where id = {value}')
        columns = [col.name for col in cursor.description]
        coordinates_current = pd.DataFrame(cursor.fetchall(), columns=columns)

        cursor.execute(f'SELECT * FROM cycle WHERE id = {value}')
        columns = [col.name for col in cursor.description]
        cycle_data = pd.DataFrame(cursor.fetchall(), columns=columns)

        cursor.execute("""
                    SELECT id, number_cycle, AVG(x) AS center_x, AVG(y) AS center_y, AVG(z) AS center_z
                    FROM (
                        SELECT cycle.id, number_cycle, atom_in_cycle, x, y, z 
                        FROM cycle
                        JOIN coordinates ON cycle.id = coordinates.id
                        WHERE cycle.id = %s AND atom_in_cycle = number_atoms
                    ) AS a
                    GROUP BY id, number_cycle;
                """, (value,))
        cycle_centers_columns = [col.name for col in cursor.description]
        cycle_centers_data = pd.DataFrame(cursor.fetchall(), columns=cycle_centers_columns)

        # Преобразование объектов Decimal во float
        cycle_centers_data['center_x'] = cycle_centers_data['center_x'].astype(float)
        cycle_centers_data['center_y'] = cycle_centers_data['center_y'].astype(float)
        cycle_centers_data['center_z'] = cycle_centers_data['center_z'].astype(float)

        # Создаем молекулу RDKit
        mol = Chem.RWMol()

        atom_indices = {}  # Словарь для отслеживания индексов атомов
        atom_colors = []  # Список для хранения цветов атомов
        atom_cycle_info = {}  # Dictionary to store cycle information for each atom

        for index, row in coordinates_current.iterrows():
            atom_num = row['number_atoms']
            element = row['element']
            x = float(row['x'])
            y = float(row['y'])
            z = float(row['z'])

            atom = Chem.Atom(element)
            atom.SetAtomMapNum(atom_num)
            atom.SetDoubleProp('x', x)
            atom.SetDoubleProp('y', y)
            atom.SetDoubleProp('z', z)
            atom_numbers.append(atom_num)
            atom_idx = mol.AddAtom(atom)
            atom_indices[atom_num] = atom_idx

            # Получаем цвет атома
            color = get_atom_color(element)
            atom_colors.append(color)

        for index, row in cycle_data.iterrows():
            atom_num = int(row['atom_in_cycle']-1)  # Convert to int
            cycle_num = int(row['number_cycle'])  # Convert to int
            atom_cycle_info[atom_num] = cycle_num

        bond_cylinders = []  # List to store bond cylinder data

            # Добавляем номер атома в список


        for index, row in connections_current.iterrows():
            atom_1 = row['atom_1']
            atom_2 = row['atom_2']
            bond_type = int(row['type'])

            if atom_1 in atom_indices and atom_2 in atom_indices and bond_type is not None:
                atom_idx_1 = atom_indices[atom_1]
                atom_idx_2 = atom_indices[atom_2]

                # Convert bond type
                if bond_type == 1:
                    rdkit_bond_type = Chem.BondType.SINGLE
                elif bond_type == 2:
                    rdkit_bond_type = Chem.BondType.DOUBLE
                elif bond_type == 3:
                    rdkit_bond_type = Chem.BondType.TRIPLE

                # Add bond to molecule
                if not mol.GetBondBetweenAtoms(atom_idx_1, atom_idx_2):
                    mol.AddBond(atom_idx_1, atom_idx_2, rdkit_bond_type)

                # Add bond data for cylinders
                bond_cylinders.append({"start": atom_idx_1, "end": atom_idx_2, "type": bond_type})

        Chem.SanitizeMol(mol)
        AllChem.EmbedMolecule(mol)

        coords = np.array([(atom.GetDoubleProp('x'), atom.GetDoubleProp('y'), atom.GetDoubleProp('z'))
                            for atom in mol.GetAtoms()])

        # Создаем объект py3Dmol.view
        view = py3Dmol.view(width=1800, height=1800)



        # Выводим данные молекулы в формате JSON
        molecule_data = {
            'coords': coords.tolist(),
            'bonds': bond_cylinders,
            'atom_colors': atom_colors,
            'atom_numbers': atom_numbers,
            'cycle_info': atom_cycle_info,
            'cycle_centers': cycle_centers_data.to_dict('records')  # Convert DataFrame to list of dictionaries

        }
        return molecule_data

    except (Exception, Error) as error:
        print(f"Ошибка при работе с PostgreSQL для молекулы {value}: {error}")

    finally:
        if connection:
            connection.close()



def get_atom_color(element):
    # Задаем цвет для каждого элемента по номеру
    element_colors  = {
        1: '#ffffff',
        2: '#d9ffff',
        3: '#cc80ff',
        4: '#c2ff00',
        5: '#ffb5b5',
        6: '#505050',
        7: '#3050f8',
        8: '#ff0d0d',
        9: '#90e050',
        10: '#b3e3f5',
        11: '#ab5cf2',
        12: '#8aff00',
        13: '#bfa6a6',
        14: '#f0c8a0',
        15: '#ff8000',
        16: '#ffff30',
        17: '#1ff01f',
        18: '#80d1e3',
        19: '#8f40d4',
        20: '#3dff00',
        21: '#e6e6e6',
        22: '#bfc2c7',
        23: '#a6a6ab',
        24: '#8a99c7',
        25: '#9c7ac7',
        26: '#e06633',
        27: '#f090a0',
        28: '#50d050',
        29: '#c88033',
        30: '#7d80b0',
        31: '#c28f8f',
        32: '#668f8f',
        33: '#bd80e3',
        34: '#ffa100',
        35: '#a62929',
        36: '#5cb8d1',
        37: '#702eb0',
        38: '#00ff00',
        39: '#94ffff',
        40: '#94e0e0',
        41: '#73c2c9',
        42: '#54b5b5',
        43: '#3b9e9e',
        44: '#248f8f',
        45: '#0a7d8c',
        46: '#006985',
        47: '#c0c0c0',
        48: '#ffd98f',
        49: '#a67573',
        50: '#668080',
        51: '#9e63b5',
        52: '#d47a00',
        53: '#940094',
        54: '#429eb0',
        55: '#57178f',
        56: '#00c900',
        57: '#70d4ff',
        58: '#ffffc7',
        59: '#d9ffc7',
        60: '#c7ffc7',
        61: '#a3ffc7',
        62: '#8fffc7',
        63: '#61ffc7',
        64: '#45ffc7',
        65: '#30ffc7',
        66: '#1fffc7',
        67: '#00ff9c',
        68: '#00e675',
        69: '#00d452',
        70: '#00bf38',
        71: '#00ab24',
        72: '#4dc2ff',
        73: '#4da6ff',
        74: '#2194d6',
        75: '#267dab',
        76: '#266696',
        77: '#175487',
        78: '#d0d0e0',
        79: '#ffd123',
        80: '#b8b8d0',
        81: '#a6544d',
        82: '#575961',
        83: '#9e4fb5',
        84: '#ab5c00',
        85: '#754f45',
        86: '#428296',
        87: '#420066',
        88: '#007d00',
        89: '#70abfa',
        90: '#00baff',
        91: '#00a1ff',
        92: '#008fff',
        93: '#0080ff',
        94: '#006bff',
        95: '#545cf2',
        96: '#785ce3',
        97: '#8a4fe3',
        98: '#a136d4',
        99: '#b31fd4',
        100: '#b31fba',
        101: '#b30da6',
        102: '#bd0d87',
        103: '#c70066',
        104: '#cc0059',
        105: '#d1004f',
        106: '#d90045',
        107: '#e00038',
        108: '#e6002e',
        109: '#eb0026',
        110: '#eb0026',
        111: '#eb0026',
        112: '#bce0eb',
        113: '#f7d8c2',
        114: '#f7d8c2',
        115: '#f6d9c4',
        116: '#f6dac0',
        117: '#f6d9c4',
        118: '#f7d8c2'
     }

    return element_colors.get(element, '#000000')  # Возвращаем черный цвет по умолчанию, если элемент не найден
def get_unique_elements(connection_params, value):
    connection = None
    try:
        connection = psycopg2.connect(**connection_params)
        cursor = connection.cursor()

        cursor.execute(f'''
            SELECT DISTINCT elements.element, elements.name_element
            FROM coordinates
            JOIN elements ON coordinates.element = elements.element
            WHERE coordinates.id = {value}
        ''')
        unique_elements = [{'element': row[0], 'name': row[1], 'color': get_atom_color(row[0])} for row in cursor.fetchall()]

        return unique_elements

    except (Exception, Error) as error:
        print(f"Error retrieving unique elements for molecule {value}: {error}")

    finally:
        if connection:
            connection.close()


def get_groups_data(connection_params, value):
    connection = None
    try:
        connection = psycopg2.connect(**connection_params)
        cursor = connection.cursor()

        cursor.execute(f'''
            SELECT number_atom, groups
            FROM groups
            WHERE id = {value}
        ''')
        groups_data = [{'number_atom': row[0], 'groups': row[1]} for row in cursor.fetchall()]

        return groups_data

    except (Exception, Error) as error:
        print(f"Error retrieving groups data for molecule {value}: {error}")

    finally:
        if connection:
            connection.close()