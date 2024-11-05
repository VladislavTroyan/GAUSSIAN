#!/usr/bin/env python
# coding: utf-8

# In[1]:


import psycopg2
from psycopg2 import Error
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from sqlalchemy import create_engine
import pandas as pd
import numpy as np
from itertools import combinations


def db(id):
    try:
        connection = psycopg2.connect(user="postgres",
                                      password="11111111",
                                      host="127.0.0.1",
                                      port="5432",
                                      database="GAUSSIAN")
        connection.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        cursor = connection.cursor()

        engine = create_engine('postgresql://postgres:11111111@localhost:5432/GAUSSIAN')

        # Возьмем все существующие группы из БД
        cursor.execute(f'select * from connections where id = {id}')
        columns = []
        for idx, col in enumerate(cursor.description):
            columns.append(col.name)

        connections_current = pd.DataFrame(cursor.fetchall(), columns=columns)
        cursor.execute(f'select * from coordinates where id = {id}')
        columns = []

        for idx, col in enumerate(cursor.description):
            columns.append(col.name)
        coordinates_current = pd.DataFrame(cursor.fetchall(), columns=columns)

        cursor.execute(f'select * from groups where id = {id}')
        columns = []
        for idx, col in enumerate(cursor.description):
            columns.append(col.name)

        groups_current = pd.DataFrame(cursor.fetchall(), columns=columns)

        cursor.execute(f'select cas from raschety where id = {id}')
        cas = cursor.fetchall()
        cas = ''.join([''.join(i) for i in cas])

    except (Exception, Error) as error:
        print("Ошибка при работе с PostgreSQL", error)

    finally:
        if connection:
            cursor.close()
            connection.close()
    return coordinates_current, connections_current, groups_current, cas


def R_COOH(coordinates_current, connections_current, groups_current):
    flag_COOH = False

    # Проверяем наличие нужных групп
    if ('C-C,O,OH' in list(groups_current['groups'])) & ('O-C' in list(groups_current['groups'])):
        # Определяем номера элементов, у которых эти группы
        num_C = list(groups_current[groups_current['groups'] == 'C-C,O,OH']['number_atom'])
        num_O = list(groups_current[groups_current['groups'] == 'O-C']['number_atom'])

        for Ci in num_C:
            for O in num_O:
                # Если есть связь между углеродом и кислородом
                if len(connections_current[
                           ((connections_current['atom_1'] == Ci) & (connections_current['atom_2'] == O)) |
                           ((connections_current['atom_2'] == Ci) & (connections_current['atom_1'] == O))]) == 1:
                    flag_COOH = True
                    break
    return flag_COOH


def R_CHO(coordinates_current, connections_current, groups_current):
    flag_CHO = False
    # Проверяем наличие нужных групп
    if 'C-C,H,O' in list(groups_current['groups']):
        flag_CHO = True
        # Восстановим номера атомов, из которых формируется альдегидная группа
        num_C = list(groups_current[groups_current['groups'] == 'C-C,H,O']['number_atom'])
        for C in num_C:
            con = connections_current[(connections_current['atom_1'] == C) | (connections_current['atom_2'] == C)]
            atom = set(list(con['atom_1']) + list(con['atom_2']))
            for N in atom:
                if int(coordinates_current[coordinates_current['number_atoms'] == N]['element']) == 8:
                    O = N
                elif int(coordinates_current[coordinates_current['number_atoms'] == N]['element']) == 1:
                    H = N
    return flag_CHO


def ester(coordinates_current, connections_current, groups_current):
    flag_ester = False
    if 'C-C,2O' in list(groups_current['groups']):
        # Определяем номера элементов, у которых эти группы
        num_C = list(groups_current[groups_current['groups'] == 'C-C,2O']['number_atom'])

        # Проходимся по всем молекулам с такой группой
        for C in num_C:
            # Определяем, с кем этот атом связан
            con = connections_current[(connections_current['atom_1'] == C) | (connections_current['atom_2'] == C)]
            atom = set(list(con['atom_1']) + list(con['atom_2']))
            # Проходимся по каждому соседу и находим кислороды
            for O in atom:
                if int(coordinates_current[coordinates_current['number_atoms'] == O]['element'].iloc[0]) == 8:
                    # Если у этого кислорода есть соседи помимо углерода, то он нам и нужен
                    if len(connections_current[
                               (connections_current['atom_1'] == O) | (connections_current['atom_2'] == O)]) > 1:
                        # Смотрим, кто второй сосед у этого кислорода
                        con_O = connections_current[
                            (connections_current['atom_1'] == O) | (connections_current['atom_2'] == O)]
                        for i in set(list(con_O['atom_1']) + list(con_O['atom_2'])):
                            if i != C and i != O:
                                # Если сосед - углерод, то говорим, что это радикальная группа, соответственно данная моделула - сложный эфир
                                if int(coordinates_current[coordinates_current['number_atoms'] == i]['element'].iloc[
                                           0]) == 6:
                                    flag_ester = True
                        break

    return flag_ester


def R_O_R(coordinates_current, connections_current, groups_current):
    flag_R_O = False
    # Проверяем наличие нужных групп
    if 'O-2C' in list(groups_current['groups']):
        flag_R_O = True
        # Восстановим номера атомов, из которых формируется альдегидная группа
        O = list(groups_current[groups_current['groups'] == 'O-2C']['number_atom'])

    return flag_R_O


def R_OH(coordinates_current, connections_current, groups_current):
    flag_OH = False
    for i, row in groups_current.iterrows():
        if ('OH' in row['groups']):
            row_not_OH = row['groups'].replace('OH', '')
            if not 'O' in row_not_OH:
                flag_OH = True
    return flag_OH


def R_N(coordinates_current, connections_current, groups_current, cycle, arom_cycle):
    flag_N_ar = False
    flag_N = False
    # Проверим наличие аминогруппы ароматического кольца
    if 'N-2C' in list(groups_current['groups']):
        num_N = list(groups_current[groups_current['groups'] == 'N-2C']['number_atom'])
        for N in num_N:
            con = connections_current[(connections_current['atom_1'] == N) | (connections_current['atom_2'] == N)]
            atom = list(set(list(con['atom_1']) + list(con['atom_2'])))
            ind_N = atom.index(N)
            atom.pop(ind_N)
            for ind, row in arom_cycle.iterrows():
                num_cycle = row['number_cycle']
                atom_in_cycle_list = list(cycle[cycle['number_cycle'] == num_cycle]['atom_in_cycle'])
                if (atom[0] in atom_in_cycle_list) and (atom[1] in atom_in_cycle_list) and (row['arom']):
                    flag_N_ar = True
                    break
            # Если группа есть, но она не входит в состав ароматического цикла
            if not flag_N_ar:
                flag_N = True

    # Проверим наличие остальных групп с азотом
    for element in list(set(list(groups_current['groups']))):
        if element.startswith('N-') and element != 'N-2C' and element != 'N-C,2O':
            flag_N = True
            break

    return flag_N_ar, flag_N


def R_NO2(coordinates_current, connections_current, groups_current):
    flag_NO2 = False
    if 'N-C,2O' in list(groups_current['groups']):
        flag_NO2 = True
    return flag_NO2


def R_Cl(coordinates_current, connections_current, groups_current):
    flag_R_Cl = False
    if ('Cl-C' in list(groups_current['groups'])):
        flag_R_Cl = True
    return flag_R_Cl


def R_Br(coordinates_current, connections_current, groups_current):
    flag_R_Br = False
    if ('Br-C' in list(groups_current['groups'])):
        flag_R_Br = True
    return flag_R_Br


def R_F(coordinates_current, connections_current, groups_current):
    flag_R_F = False
    if ('F-C' in list(groups_current['groups'])):
        flag_R_F = True
    return flag_R_F


def R_I(coordinates_current, connections_current, groups_current):
    flag_R_I = False
    if ('I-C' in list(groups_current['groups'])):
        flag_R_I = True
    return flag_R_I


def R_SH(coordinates_current, connections_current, groups_current):
    flag_R_SH = False
    if ('S-C,H' in list(groups_current['groups'])):
        flag_R_SH = True
    return flag_R_SH


def R_SOOOH(coordinates_current, connections_current, groups_current):
    flag_R_SOOOH = False
    if ('S-C,3O' in list(groups_current['groups'])) and ('O-H' in list(groups_current['groups'])) and (
            'O-S' in list(groups_current['groups'])):
        num_S = list(groups_current[groups_current['groups'] == 'O-C,H']['number_atom'])
        for S in num_S:
            flag_O_H = False
            flag_OS = False
            con = connections_current[(connections_current['atom_1'] == S) | (connections_current['atom_2'] == S)]
            atom = set(list(con['atom_1']) + list(con['atom_2']))
            # Удаляем саму серу, чтобы делать меньше иттераций
            atom.remove(S)
            OS = 0
            for O in atom:
                if int(coordinates_current[coordinates_current['number_atoms'] == O]['element'].iloc[0]) == 8:
                    if groups_current[groups_current['number_atom'] == O]['groups'] == 'O-H':
                        flag_O_H = True
                    elif groups_current[groups_current['number_atom'] == O]['groups'] == 'O-S':
                        OS += 1
            if OS == 2:
                flag_OS = True

            if flag_O_H and flag_OS:
                flag_R_SOOOH = True
                break

    return flag_R_SOOOH


def R_CO_R(coordinates_current, connections_current, groups_current):
    flag_R_CO_R = False
    if ('C-2C,O' in list(groups_current['groups'])) and ('O-C' in list(groups_current['groups'])):
        # Определяем номера элементов, у которых эти группы
        num_C = list(groups_current[groups_current['groups'] == 'C-2C,O']['number_atom'])
        num_O = list(groups_current[groups_current['groups'] == 'O-C']['number_atom'])
        for C in num_C:
            for O in num_O:
                if len(connections_current[
                           ((connections_current['atom_1'] == C) & (connections_current['atom_2'] == O)) |
                           ((connections_current['atom_1'] == O) & (connections_current['atom_2'] == C))]) > 0:
                    flag_R_CO_R = True
                    break

    return flag_R_CO_R


def aromatic_cycle(coor):
    # Преобразуем координаты в массив
    arr = []
    for i, row in coor.iterrows():
        arr.append([float(row['x']), float(row['y']), float(row['z'])])

    ar = False

    keyList = [0, 1, 2]
    mean_residual = {key: None for key in keyList}
    coef = {key: None for key in keyList}

    # Проверяем только те циклы, которые состоят из более 5 точек
    if len(arr) > 4:
        points = np.array(arr)
        # Если значения в одном из столбцов равны - атомы уже лежат в одной плоскости
        all_equal = False
        for column in points.T:
            if all(x == column[0] for x in column):
                all_equal = True
                break
        # Если не равны, пытаемся построить плоскость МНК
        if all_equal:
            ar = True
        else:
            #  1 случай:
            A = np.c_[points[:, 0], points[:, 1], np.ones(points.shape[0])]
            b = -points[:, 2]
            x_coef, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
            if rank > 2:
                coef[0] = x_coef
                mean_residual[0] = residuals
                if mean_residual[0] < 0.089:
                    ar = True

            # 2 случай:
            A = np.c_[points[:, 0], points[:, 2], np.ones(points.shape[0])]
            b = -points[:, 1]
            x_coef, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
            if rank > 2:
                coef[1] = x_coef
                mean_residual[1] = residuals
                if mean_residual[1] < 0.089:
                    ar = True

            # 3 случай:
            A = np.c_[points[:, 1], points[:, 2], np.ones(points.shape[0])]
            b = -points[:, 0]
            x_coef, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
            if rank > 2:
                coef[2] = x_coef
                mean_residual[2] = residuals
                if mean_residual[2] < 0.089:
                    ar = True
    #     print(mean_residual)
    return ar


def find_all_paths(adj_matrix, start_node, end_node):
    visited = [False] * len(adj_matrix)
    current_path = []
    all_paths = []
    find_paths(adj_matrix, start_node, end_node, visited, current_path, all_paths)
    return all_paths


def find_paths(adj_matrix, current_node, end_node, visited, current_path, all_paths):
    visited[current_node] = True
    current_path.append(current_node + 1)
    if current_node == end_node:
        all_paths.append(list(current_path))
    else:
        for neighbor, is_connected in enumerate(adj_matrix[current_node]):
            if is_connected and not visited[neighbor]:
                find_paths(adj_matrix, neighbor, end_node, visited, current_path, all_paths)
    current_path.pop()
    visited[current_node] = False


def find_cycle(coordinates_current, connections_current, ID):
    df = connections_current.merge(coordinates_current, right_on='number_atoms', left_on='atom_1')[
        ['atom_1', 'atom_2', 'element']].rename(columns={'element': 'element_1'})
    df = df.merge(coordinates_current, right_on='number_atoms', left_on='atom_2')[
        ['atom_1', 'atom_2', 'element_1', 'element']].rename(columns={'element': 'element_2'})
    # Сформировать граф
    n = max(coordinates_current['number_atoms'])
    grahf = np.zeros((n, n))
    for i in range(len(connections_current)):
        atom_1 = df.loc[i, 'atom_1']
        atom_2 = df.loc[i, 'atom_2']
        el_1 = df.loc[i, 'element_1']
        el_2 = df.loc[i, 'element_2']
        flag = False
        for index, row in enumerate(grahf):
            if flag:
                break
            else:
                if index == atom_1 - 1:
                    for column, value in enumerate(row):
                        if (column == atom_2 - 1):
                            grahf[index][column] = 1
                            grahf[column][index] = 1
                            flag = True
                            break

    # Запустить поиск всех путей
    list_itog = []
    for i in range(len(coordinates_current)):
        # Смотрим для кажлого атома его соседей
        s = np.nonzero(grahf[i])[0]
        for j in s:
            # Принимаем текущую вершину и его соседа - за начальную и конечную вершину
            start_node = i + 1
            end_node = j + 1
            all_paths = find_all_paths(grahf, start_node - 1, end_node - 1)
            # Если путей больше, чем один - есть цикл, через который может пройти алгоритм
            if len(all_paths) > 1:
                for path in all_paths:
                    # Уберем ненужные данные
                    if len(path) > 2:
                        # Отсортируем все массивы
                        path.sort()
                        # Если такого пути нет, запишем его в общий список
                        if not path in list_itog:
                            list_itog.append(path)
    list_cycle = list_itog.copy()
    # Удалим те циклы, которые включают в себя 2 малых цикла
    for i in list_itog:
        num_i = list_itog.index(i)
        list_i = list(range(len(list_itog)))
        list_i.pop(num_i)
        comb = list(combinations(list_i, 2))
        for c in comb:
            m = list_itog[c[0]]
            n = list_itog[c[1]]
            merge = m + n
            merge = list(set(merge))
            merge.sort()
            if list_itog[num_i] == merge:
                num_2 = list_cycle.index(merge)
                list_cycle.pop(num_2)
                break
    # Удалим циклы, которые включают в себя 3 малых цикла
    all_c = list_cycle.copy()
    for cycle in list_cycle:
        comb = list(combinations(list_cycle, 3))
        for c in comb:
            merge = c[0] + c[1] + c[2]
            merge = list(set(merge))
            m1 = set(c[0])
            m2 = set(c[1])
            m3 = set(c[2])
            common_elements = m1.intersection(m2, m3)
            cycle_new = [x for x in merge if x not in common_elements]
            cycle_new.sort()
            if cycle == cycle_new:
                num_2 = all_c.index(cycle)
                all_c.pop(num_2)
                break
    cycle = pd.DataFrame(columns=['id', 'number_cycle', 'atom_in_cycle'])
    arom_cycle = pd.DataFrame(columns=['id', 'number_cycle', 'arom'])
    n = 1
    print(all_c)
    for c in all_c:
        coor = coordinates_current[coordinates_current['number_atoms'].isin(c)]
        number_atom_in_cycle = list(coor['number_atoms'])
        count_neighbors = 0
        for atom_in_cycle in list(coor['number_atoms']):
            for i, row in connections_current.iterrows():
                if (atom_in_cycle == row['atom_1'] and not row['atom_2'] in number_atom_in_cycle) or (
                        atom_in_cycle == row['atom_2'] and not row['atom_1'] in number_atom_in_cycle):
                    count_neighbors += 1

        # Элементы из которых он состоит
        elements_in_cycle = list(
            coordinates_current.loc[coordinates_current['number_atoms'].isin(number_atom_in_cycle)]['element'])

        # Если цикл чисто из углеродов - проверяем, сколько связей у всех элементов из цикла
        if all(x == 6 for x in elements_in_cycle):
            if len(number_atom_in_cycle) == count_neighbors:
                ar = aromatic_cycle(coor)
            else:
                ar = False
        else:
            ar = aromatic_cycle(coor)
        for c_i in c:
            lenght_cycle = len(cycle)
            cycle.loc[lenght_cycle, 'id'] = ID
            cycle.loc[lenght_cycle, 'number_cycle'] = n
            cycle.loc[lenght_cycle, 'atom_in_cycle'] = c_i
        lenght_arom = len(arom_cycle)
        arom_cycle.loc[lenght_arom, 'id'] = ID
        arom_cycle.loc[lenght_arom, 'number_cycle'] = n
        arom_cycle.loc[lenght_arom, 'arom'] = ar

        n += 1
    return cycle, arom_cycle

# In[ ]:




