
import re
from collections import Counter, defaultdict
from datetime import datetime
from sqlalchemy import text
from .functional_groups_function import *

def files_processing(file_content, user_db, password_db, name_database):
    # Переменная с ошибками в файлах
    ErrorFile = []
    error_type = None
    try:
        connection = psycopg2.connect(user='postgres',
                                      password='11111111',
                                      host="127.0.0.1",
                                      port="5432",
                                      database='GAUSSIAN')
        connection.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        print('Соединение с базой данных открыто')
        cursor = connection.cursor()

        engine = create_engine(f'postgresql://postgres:11111111@localhost:5432/GAUSSIAN', max_overflow=-1)

        # Возьмем все элементы, которые используются в молекулах
        cursor.execute('select * from elements')
        columns = []
        for idx, col in enumerate(cursor.description):
            columns.append(col.name)
        elements = pd.DataFrame(cursor.fetchall(), columns = columns)

        # Методы, по которым возможно добавить новые записи в базу
        cursor.execute('select * from methods')
        columns = []

        for idx, col in enumerate(cursor.description):
            columns.append(col.name)
        Methods = pd.DataFrame(cursor.fetchall(), columns = columns)

        # Проходимся по каждому файлу, в указанном каталоге
        itog = 0
        flag = True
        filename = file_content.name
        print(filename)
        try:
            # Берем все функциональные группы, которые есть в базе данных
            cursor.execute(f'select * from functional_groups')
            columns = []
            for idx, col in enumerate(cursor.description):
                columns.append(col.name)
            functional_groups = pd.DataFrame(cursor.fetchall(), columns=columns)

            # Берем все группы, которые есть в базе данных
            cursor.execute(f'select groups from contribution')
            columns = []
            for idx, col in enumerate(cursor.description):
                columns.append(col.name)
            contribution = pd.DataFrame(cursor.fetchall(), columns=columns)

            all_groups_in_contribution = contribution['groups'].to_list()

            # Последний используемый ID возьмем из базы
            cursor.execute('select count(id) from raschety')
            ID = cursor.fetchall()[0][0] + 1

            # Читаем файл из StringIO
            file_content.seek(0)  # Возвращаем курсор в начало файла
            text = file_content.read().decode('utf-8')

            # Берем только те файлы,с завершенным расчетом (в конце всегда есть дата выполнения)
            if not re.findall(r'((Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s*\d+\s*\d+:\d+:\d+\s*\d{4})',
                              text[-30:]):
                error_type = 1
            else:

                # Проверяем наличие CAS в файле:
                if not re.findall(r'\s*?CAS\s*?(\d+\-\d+\-\d+)', text):
                    error_type = 2
                    # Вот тут надо запрашивать кас, если его нет
                else:
                    # Вот тут надо спрашивать, верно ли найден кас

                    # Инициализируем переменные пустыми
                    duration = None
                    CAS = None
                    method = None
                    date = None
                    method_number = None

                    # Ищем в тексте CAS, метод, дату расчета и длительность
                    CAS = re.findall(r'\s*?CAS\s*?(\d+\-\d+\-\d+)', text)[-1]
                    method = find_method(text)
                    date = datetime.strptime(
                        re.findall('((Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\\s*\\d+\\s*\\d+:\\d+:\\d+\\s*\\d{4})',
                                   text)[-1][0], '%b %d %H:%M:%S %Y').strftime('%d-%m-%Y %H:%M:%S')
                    duration = find_duration(text)
                    print(CAS, method)

                    temp_and_press = text.split('Thermochemistry')[1].split('Atm.')[0]
                    if len(temp_and_press) > 0:
                        temperature = float(temp_and_press.split('Temperature')[1].split('Kelvin')[0].replace(' ', ''))
                        pressure = float(temp_and_press.split('Pressure')[1].split('Kelvin')[0].replace(' ', ''))
                    else:
                        temperature = 298
                        pressure = 1
                    # Если найдены CAS, метод и дата, то обрабатываем файл в соответствии с методом
                    if CAS and method and date:
                        if method == 'B3LYP':
                            method_number = int(Methods[Methods['name_method'] == method]['method'])
                            [coordinates_current, new_row_energy] = fun_B3LYP(text, method, ID)

                        elif method == 'G4':
                            method_number = int(Methods[Methods['name_method'] == method]['method'])
                            [coordinates_current, new_row_energy] = fun_G4(text, method, ID, elements)

                        elif method == 'G3MP2':
                            method_number = int(Methods[Methods['name_method'] == method]['method'].iloc[0])
                            [coordinates_current, new_row_energy] = fun_G3MP2(text, method, ID, elements)

                        elif method == 'M06':
                            method_number = int(Methods[Methods['name_method'] == method]['method'])
                            [coordinates_current, new_row_energy] = fun_M06(text, method, ID)

                        else:
                            error_type = 3
                            flag = False

                        if flag:
                            connection_current = new_connection(text, ID)
                            new_row_molecule = new_Molecule(text, method, CAS)
                            dict_group = grouping(coordinates_current, connection_current, ID, elements)
                            cycle, arom_cycle = find_cycle(coordinates_current, connection_current, ID)
                            connection_current = type_connections(connection_current, coordinates_current, cycle,
                                                                  arom_cycle)

                            if len(coordinates_current) > 0 and len(connection_current) > 0 and len(
                                    new_row_energy) == 1 and str(new_row_energy['energy_gibbs']) != 'nan':
                                new_row_raschety = pd.DataFrame([{'id': ID,
                                                                  'cas': CAS,
                                                                  'method': method_number,
                                                                  'date': date,
                                                                  'duration': duration,
                                                                  'temperature': temperature,
                                                                  'pressure': pressure,
                                                                  'date_added': datetime.now().strftime(
                                                                      '%d-%m-%Y %H:%M:%S')}])

                                cas = '\'' + CAS + '\''
                                cursor.execute(fr'select * from molecules where cas = {cas}')
                                columns = []
                                for idx, col in enumerate(cursor.description):
                                    columns.append(col.name)
                                molecules = pd.DataFrame(cursor.fetchall(), columns=columns)

                                cursor.execute(fr'select * from fg_molecules where cas = {cas}')
                                columns = []
                                for idx, col in enumerate(cursor.description):
                                    columns.append(col.name)
                                fg = pd.DataFrame(cursor.fetchall(), columns=columns)
                                if len(molecules) == 0:
                                    new_row_molecule.to_sql(con=engine, name='molecules', index=False,
                                                            if_exists='append')
                                    print('molecules добавлена')

                                if len(fg) == 0:
                                    new_row_fg = functional_groups_cas(coordinates_current, connection_current,
                                                                       dict_group, CAS, functional_groups, cycle,
                                                                       arom_cycle)
                                    new_row_fg.to_sql(con=engine, name='fg_molecules', index=False, if_exists='append')
                                    print('fg_molecules добавлена')

                                groups_in_current_raschety = set(dict_group['groups'].to_list())

                                for gr in set(groups_in_current_raschety):
                                    if not gr in set(all_groups_in_contribution):
                                        pd.DataFrame([{'groups': gr, 'vklad': np.nan}]).to_sql(con=engine,
                                                                                               name='contribution',
                                                                                               index=False,
                                                                                               if_exists='append')

                                try:
                                    method = method.lower()
                                    new_row_raschety.to_sql(con=engine, name='raschety', index=False,
                                                            if_exists='append')
                                    coordinates_current.to_sql(con=engine, name='coordinates', index=False,
                                                               if_exists='append')
                                    arom_cycle.to_sql(con=engine, name='arom_cycle', index=False, if_exists='append')
                                    cycle.to_sql(con=engine, name='cycle', index=False, if_exists='append')
                                    connection_current.to_sql(con=engine, name='connections', index=False,
                                                              if_exists='append')
                                    new_row_energy.to_sql(con=engine, name=f'{method}', index=False, if_exists='append')
                                    dict_group.to_sql(con=engine, name='groups', index=False, if_exists='append')
                                    cursor.execute(fr'SELECT calculation_new_dfH({ID})')

                                    print('Расчет добавлен в базу')
                                    # Подтверждение транзакции
                                    connection.commit()
                                    itog += 1
                                except (Exception, Error) as error:
                                    error_type = 4
                                    # Откат транзакции в случае ошибки
                                    connection.rollback()



        except Exception as error:
            error_type = 5
            ErrorFile.append([error, filename])
            print(error)
        if error_type is not None:
            print(filename)
            print('Код ошибки: ', error_type)
            return {'error_type': error_type, 'filename': filename}
        connection.commit()

    except (Exception, Error) as error:
        print("Ошибка при работе с PostgreSQL", error)

    finally:
        if connection:
            cursor.close()
            connection.close()
            print("Соединение с закрыто")
    print('Добавлено в базу данных: ', itog)
    return ErrorFile


# error


#     ------------------------------------------------------------------------------------------------------------------------

def fun_G4(text, method, ID, elements):
    enthalpy = None
    EGibbs = None
    new_row_energy = pd.DataFrame()
    coordinates_current = pd.DataFrame()

    if re.findall(rf"{method}\s*Energy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE):

        # В 16 версии в методе G4 идет смещение значений
        if bool(re.search("gaussian/g16", text).group()):
            enthalpy = re.findall(
                rf"{method}\s*Energy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
            )[-1]
            EGibbs = re.findall(rf"{method}\s*Enthalpy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE)[-1]

        # В других всё хорошо
        else:
            enthalpy = re.findall(
                rf"{method}\s*Enthalpy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
            )[-1]
            EGibbs = re.findall(
                rf"{method}\s*Free Energy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
            )[-1]

        text_coordinates = (
            text.split("Standard orientation:")[-1]
            .split("Rotational constants")[0]
            .split("-\n")[2]
            .split("-----")[0]
            .split("\n")[:-1]
        )
        coordinates_current = pd.DataFrame(
            columns=["number_atoms", "element", "x", "y", "z"]
        )

        for i, row in enumerate(text_coordinates):
            row = re.sub(r"\s+", " ", row).lstrip().split(" ")
            row.pop(2)
            coordinates_current.loc[len(coordinates_current)] = row
        coordinates_current["id"] = ID
        coordinates_current["number_atoms"] = coordinates_current[
            "number_atoms"
        ].astype(int)
        coordinates_current["element"] = coordinates_current["element"].astype(int)

        # Атомизация G4

        count_elements = pd.DataFrame(coordinates_current.groupby(['element'])['element'].count()).rename(
            columns={'element': 'count'}).reset_index().merge(elements, left_on='element',
                                                              right_on='element').set_index('symbol')['count'].to_dict()
        count_elements = defaultdict(lambda: 0, count_elements)
        at = count_elements['C'] * 711.185 + count_elements['H'] * 216.035 + count_elements['O'] * 246.79 + \
             count_elements['N'] * 470.82 + count_elements['S'] * 274.735 + count_elements['F'] * 77.284 + \
             count_elements['Cl'] * 119.621 + count_elements['Br'] * 117.917 - 2625.5 * (
                         count_elements['C'] * -37.83417 + count_elements['H'] * -0.50142 + count_elements[
                     'O'] * -75.0455 + count_elements['N'] * -54.57367 + count_elements['S'] * -397.98018 +
                         count_elements['F'] * -99.705 + count_elements['Cl'] * -460.015 + count_elements[
                             'Br'] * -2573.5854) + 2625.5 * float(enthalpy) - (
                         count_elements['C'] * 1.051 + count_elements['H'] * 4.234 + count_elements['O'] * 4.342 +
                         count_elements['N'] * 4.335 + count_elements['S'] * 4.412 + count_elements['F'] * 4.4125 +
                         count_elements['Cl'] * 4.5905 + count_elements['Br'] * 12.2545)
        entropy = 2625.5 * 1000 * (float(enthalpy) - float(EGibbs)) / 298
        new_row_energy = pd.DataFrame(
            [{"id": ID, "enthalpy": enthalpy, "energy_gibbs": EGibbs, 'at': at, 'entropy': entropy}]
        )

        print('fun_G4 отработала')
    return coordinates_current, new_row_energy


def fun_G3MP2(text, method, ID, elements):
    enthalpy = None
    EGibbs = None
    new_row_energy = pd.DataFrame()
    coordinates_current = pd.DataFrame()

    if re.findall(rf"{method}\s*Energy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE):
        enthalpy = re.findall(
            rf"{method}\s*Enthalpy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
        )[-1]
        EGibbs = re.findall(
            rf"{method}\s*Free Energy=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
        )[-1]

        text_coordinates = (
            text.split("Standard orientation:")[-1]
            .split("Rotational constants")[0]
            .split("-\n")[2]
            .split("-----")[0]
            .split("\n")[:-1]
        )
        coordinates_current = pd.DataFrame(
            columns=["number_atoms", "element", "x", "y", "z"]
        )
        for i, row in enumerate(text_coordinates):
            row = re.sub(r"\s+", " ", row).lstrip().split(" ")
            row.pop(2)
            coordinates_current.loc[len(coordinates_current)] = row
        coordinates_current["id"] = ID
        coordinates_current["number_atoms"] = coordinates_current[
            "number_atoms"
        ].astype(int)
        coordinates_current["element"] = coordinates_current["element"].astype(int)
        # Группировка для атомизации
        from collections import defaultdict

        count_elements = pd.DataFrame(coordinates_current.groupby(['element'])['element'].count()).rename(
            columns={'element': 'count'}).reset_index().merge(elements, left_on='element',
                                                              right_on='element').set_index('symbol')['count'].to_dict()
        count_elements = defaultdict(lambda: 0, count_elements)
        at = count_elements['C'] * 711.185 + count_elements['H'] * 216.035 + count_elements['O'] * 246.79 + \
             count_elements['N'] * 470.82 + count_elements['S'] * 274.735 + count_elements['F'] * 77.284 + \
             count_elements['Cl'] * 119.621 - 2625.5 * (
                         count_elements['C'] * -37.78934 + count_elements['H'] * -0.50184 + count_elements[
                     'O'] * -74.98977 + count_elements['N'] * -54.52519 + count_elements['S'] * -397.66376 +
                         count_elements['F'] * -99.6409 + count_elements['Cl'] * -459.687) + 2625.5 * float(
            enthalpy) - (count_elements['C'] * 1.051 + count_elements['H'] * 4.234 + count_elements['O'] * 4.342 +
                         count_elements['N'] * 4.335 + count_elements['S'] * 4.412 + count_elements['F'] * 4.4125 +
                         count_elements['Cl'] * 4.5905)
        entropy = 2625.5 * 1000 * (float(enthalpy) - float(EGibbs)) / 298
        new_row_energy = pd.DataFrame(
            [{"id": ID, "enthalpy": enthalpy, "energy_gibbs": EGibbs, 'at': at, 'entropy': entropy}]
        )
        print('fun_G3MP2 отработала')
    return coordinates_current, new_row_energy


def fun_B3LYP(text, method, ID):
    correction_enthalpy = None
    correction_EGibbs = None
    EGibbs = None
    new_row_energy = pd.DataFrame()
    coordinates_current = pd.DataFrame()

    if bool(
            re.findall(
                rf"SCF\s*Done:\s*E\(R{method}\)\s*=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
            )
    ) | bool(
        re.findall(rf"correction to Enthalpy=\s*(-*\d+\.\d+)", text, re.IGNORECASE)
    ):
        EGibbs = re.findall(
            rf"SCF\s*Done:\s*E\(R{method}\)\s*=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
        )[-1]
        correction_enthalpy = re.findall(
            rf"correction to Enthalpy=\s*(-*\d+\.\d+)", text, re.IGNORECASE
        )[-1]
        correction_EGibbs = re.findall(
            rf"correction to Gibbs Free Energy=\s*(-*\d+\.\d+)", text, re.IGNORECASE
        )[-1]
        new_row_energy = pd.DataFrame(
            [
                {
                    "id": ID,
                    "energy_gibbs": EGibbs,
                    "correction_enthalpy": correction_enthalpy,
                    "correction_energy": correction_EGibbs,
                }
            ]
        )

        text_coordinates = (
            text.split("Input orientation:")[-1]
            .split("Distance matrix")[0]
            .split("-\n")[2]
            .split("-----")[0]
            .split("\n")[:-1]
        )
        coordinates_current = pd.DataFrame(
            columns=["number_atoms", "element", "x", "y", "z"]
        )

        for i, row in enumerate(text_coordinates):
            row = re.sub(r"\s+", " ", row).lstrip().split(" ")
            row.pop(2)
            coordinates_current.loc[len(coordinates_current)] = row
        coordinates_current["id"] = ID
        coordinates_current["number_atoms"] = coordinates_current[
            "number_atoms"
        ].astype(int)
        coordinates_current["element"] = coordinates_current["element"].astype(int)
        print('fun_B3LYP отработала')
    return coordinates_current, new_row_energy


def fun_M06(text, method, ID):
    correction_enthalpy = None
    correction_EGibbs = None
    EGibbs = None
    new_row_energy = pd.DataFrame()
    coordinates_current = pd.DataFrame()

    if bool(
            re.findall(
                rf"SCF\s*Done:\s*E\(R{method}\)\s*=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
            )
    ) | bool(
        re.findall(rf"correction to Enthalpy=\s*(-*\d+\.\d+)", text, re.IGNORECASE)
    ):
        EGibbs = re.findall(
            rf"SCF\s*Done:\s*E\(R{method}\)\s*=\s*(-*\d+\.\d+)\s*", text, re.IGNORECASE
        )[-1]
        correction_enthalpy = re.findall(
            rf"correction to Enthalpy=\s*(-*\d+\.\d+)", text, re.IGNORECASE
        )[-1]
        correction_EGibbs = re.findall(
            rf"correction to Gibbs Free Energy=\s*(-*\d+\.\d+)", text, re.IGNORECASE
        )[-1]
        new_row_energy = pd.DataFrame(
            [
                {
                    "id": ID,
                    "energy_gibbs": EGibbs,
                    "correction_enthalpy": correction_enthalpy,
                    "correction_energy": correction_EGibbs,
                }
            ]
        )

        text_coordinates = (
            text.split("Input orientation:")[-1]
            .split("Distance matrix")[0]
            .split("-\n")[2]
            .split("-----")[0]
            .split("\n")[:-1]
        )
        coordinates_current = pd.DataFrame(
            columns=["number_atoms", "element", "x", "y", "z"]
        )

        for i, row in enumerate(text_coordinates):
            row = re.sub(r"\s+", " ", row).lstrip().split(" ")
            row.pop(2)
            coordinates_current.loc[len(coordinates_current)] = row
        coordinates_current["id"] = ID
        coordinates_current["number_atoms"] = coordinates_current[
            "number_atoms"
        ].astype(int)
        coordinates_current["element"] = coordinates_current["element"].astype(int)
        print('fun_M06 отработала')
    return coordinates_current, new_row_energy


def new_connection(text, ID):
    connection_current = pd.DataFrame()
    text_connection = (
        text.split("Optimized Parameters")[-1]
        .split("A1")[0]
        .split("-\n")[2]
        .split("\n")[:-1]
    )
    connection_current = pd.DataFrame(columns=["atom_1", "atom_2", "length"])
    for i, row in enumerate(text_connection):
        row = (
            re.sub("\\s+", ",", re.sub(r"!", "", row))
            .lstrip()
            .replace("R(", "")
            .replace(")", "")
            .replace("R", "")
            .split(",")
        )
        row = row[2:5]
        connection_current.loc[len(connection_current)] = row
    connection_current["id"] = ID
    connection_current[["atom_1", "atom_2"]] = connection_current[
        ["atom_1", "atom_2"]
    ].astype(int)
    print('new_connection отработала')
    return connection_current


def new_Molecule(text, method, CAS):
    molecule_name = (
        str(re.findall(rf"\n\s*Output=(.*)", text))
        .split(".")[0]
        .split("_")[0]
        .split("'")[1]
    )
    formula = re.findall(r"\n\s*Stoichiometry\s*(.*?)\s*\n", text)[-1]
    new_row_molecule = pd.DataFrame(
        [{"cas": CAS, "formula": formula, "name_molecules": molecule_name}]
    )
    print('new_Molecule отработала')
    return new_row_molecule


def find_method(text):
    method = []
    method_find = re.findall(
        r"#\s*(.*?)\s*geom|#.*?freq\s*(.*?)\/", text
    )
    for met in method_find:
        method.append("".join(list(filter(None, method_find[0]))))
    method = "".join(set(method)).split("_")[0].split(" ")[0].split("-")[0]

    #     if re.findall('b3lyp', text, re.IGNORECASE):
    #         method = 'B3LYP'
    #     elif re.findall('g4', text, re.IGNORECASE):
    #         method = 'G4'
    #     elif re.findall('g3mp2', text, re.IGNORECASE):
    #         method = 'G3MP2'
    #     elif re.findall('m06', text, re.IGNORECASE):
    #         method = 'M06'

    print('find_method отработала')
    return method.upper()


def find_duration(text):
    duration = re.findall("Elapsed time:(.*?)\n", text)
    duration = [item.strip() for item in duration]
    day = 0
    hours = 0
    minutes = 0
    seconds = 0
    for item in duration:
        day = day + int("".join(set(re.findall("(\\d+)\\s*days", item))))
        hours = hours + int("".join(set(re.findall("(\\d+)\\s*hours", item))))
        minutes = minutes + int("".join(set(re.findall("(\\d+)\\s*minutes", item))))
        seconds = seconds + float(
            "".join(set(re.findall("(\\d+\\.\\d*)\\s*seconds", item)))
        )
    seconds = day * (24 * 3600) + hours * 3600 + minutes * 60 + seconds

    day = seconds // (24 * 3600)
    sec = seconds % (24 * 3600)
    hour = sec // 3600
    sec %= 3600
    min = sec // 60
    sec %= 60
    return rf"{int(day)} day {int(hour)} hours {int(min)} minutes {round(float(sec))} seconds"


def grouping(coordinates, connections, id, elements):
    Atoms = pd.DataFrame(coordinates[(coordinates['id'] == id) & ((coordinates['element'] == 6) |
                                                                  (coordinates['element'] == 8) | (
                                                                              coordinates['element'] == 7) |
                                                                  (coordinates['element'] == 16) | (
                                                                              coordinates['element'] == 17) |
                                                                  (coordinates['element'] == 9) | (
                                                                              coordinates['element'] == 35) |
                                                                  (coordinates['element'] == 53))][
                             'number_atoms'].drop_duplicates())

    # Все связи между атомами в молекуле
    con = connections[connections['id'] == id][['atom_1', 'atom_2']].drop_duplicates()
    con2 = con.copy()
    con2.columns = ['atom_2', 'atom_1']
    connection_all = pd.concat([con, con2]).reset_index(drop=True)

    # Связи только с теми атомами, которые есть в Atoms
    connection_atoms = pd.merge(Atoms, connection_all, how='inner', right_on='atom_1', left_on='number_atoms')
    connection_atoms = connection_atoms[['atom_1', 'atom_2']]

    # Таблица с координатами и условным обозначением элемента
    coordinates_2 = pd.merge(coordinates[coordinates['id'] == id][['number_atoms', 'element']].drop_duplicates(),
                             elements,
                             left_on='element',
                             right_on='element',
                             how='left').set_index('number_atoms')

    # Список всех атомов, у которых мы будем смотреть группы
    a = list(set(connection_atoms['atom_1']))

    dict_group = pd.DataFrame()

    # Проходимся по каждому атому:
    for i in a:
        df = list(set(connection_atoms[connection_atoms['atom_1'] == i]['atom_2']))
        l = list()

        # Записываем все элементы, с которым соединен атом
        for j in df:
            if int(coordinates[coordinates['number_atoms'] == j]['element'].iloc[0]) == 8:
                con = connections[(connections['atom_1'] == j) | (connections['atom_2'] == j)]
                atom = set(list(con['atom_1']) + list(con['atom_2']))
                atom.remove(j)
                atom.remove(i)
                if (len(atom) == 1):
                    atom_n = list(atom)[0]
                    if int(coordinates.loc[coordinates['number_atoms'] == atom_n, 'element'].iloc[0]) == 1:
                        l.append('OH')
                    else:
                        l.append('O')
                else:
                    l.append('O')
            else:
                l.append(coordinates_2.loc[j, 'symbol'])

        d = dict(sorted(dict(Counter(l)).items()))  # Считаем кол-во каждых элементов и сортируем по названию элемента

        # Записываем все в одну строку:
        atom = coordinates_2.loc[i, 'symbol']
        res = f'{atom}-'
        for item in d:
            if str(d[item]) != '1':
                res += str(d[item]) + item + ','
            else:
                res += item + ','
        lenght = len(dict_group)
        dict_group.loc[lenght, 'id'] = id
        dict_group.loc[lenght, 'number_atom'] = i
        dict_group.loc[lenght, 'groups'] = res.rstrip(',')
    print('grouping отработала')
    return dict_group


def functional_groups_cas(coordinates_current, connection_current, groups_current, cas, functional_groups, cycle,
                          arom_cycle):
    new_row_fg = pd.DataFrame(columns=['cas', 'num_fg'])
    flag_COOH = R_COOH(coordinates_current, connection_current, groups_current)
    flag_CHO = R_CHO(coordinates_current, connection_current, groups_current)
    flag_ester = ester(coordinates_current, connection_current, groups_current)
    flag_R_O_R = R_O_R(coordinates_current, connection_current, groups_current)
    flag_OH = R_OH(coordinates_current, connection_current, groups_current)
    flag_N_ar, flag_N = R_N(coordinates_current, connection_current, groups_current, cycle, arom_cycle)
    flag_NO2 = R_NO2(coordinates_current, connection_current, groups_current)
    flag_R_Cl = R_Cl(coordinates_current, connection_current, groups_current)
    flag_R_Br = R_Br(coordinates_current, connection_current, groups_current)
    flag_R_F = R_F(coordinates_current, connection_current, groups_current)
    flag_R_I = R_I(coordinates_current, connection_current, groups_current)
    flag_R_SOOOH = R_SOOOH(coordinates_current, connection_current, groups_current)
    flag_R_SH = R_SH(coordinates_current, connection_current, groups_current)
    flag_R_CO_R = R_CO_R(coordinates_current, connection_current, groups_current)
    g = arom_cycle.groupby('arom').size()
    flag_cycle = False
    flag_cycle_arom = False
    if True in g.keys():
        flag_cycle_arom = True
    if False in g.keys():
        flag_cycle = True
    list_flag = {'r_cooh': flag_COOH,
                 'r_cho': flag_CHO,
                 'r_co': flag_R_CO_R,
                 'r_coh': flag_OH,
                 'r_no2': flag_NO2,
                 'r_sh': flag_R_SH,
                 'r_o_r': flag_R_O_R,
                 'r_coo_r': flag_ester,
                 'r_so3h': flag_R_SOOOH,
                 'r_cl': flag_R_Cl,
                 'r_br': flag_R_Br,
                 'r_f': flag_R_F,
                 'r_i': flag_R_I,
                 'r_n_arom': flag_N_ar,
                 'r_n': flag_N,
                 'cycle': flag_cycle,
                 'cycle_arom': flag_cycle_arom
                 }
    for item in list_flag.items():
        if item[1]:
            name_fg = item[0]
            num_fg = functional_groups[functional_groups['symbol_fg'] == name_fg]['num_fg']
            new_row_fg.loc[len(new_row_fg)] = [cas, int(num_fg.iloc[0])]
    print('functional_groups_cas отработала')
    return new_row_fg


def new_graph(coordinates_current, connections_current):
    n = max(coordinates_current['number_atoms'])
    graph = np.zeros((n, n))
    for i in range(len(connections_current)):
        atom_1 = connections_current.loc[i, 'atom_1']
        atom_2 = connections_current.loc[i, 'atom_2']
        type_i = connections_current.loc[i, 'type']
        flag = False
        for index, row in enumerate(graph):
            if flag:
                break
            else:
                if index == atom_1 - 1:
                    for column, value in enumerate(row):
                        if column == atom_2 - 1:
                            if np.isnan(type_i):
                                graph[index][column] = 1
                                graph[column][index] = 1
                                flag = True
                            else:
                                graph[index][column] = type_i
                                graph[column][index] = type_i
                                flag = True
                            break
    return graph


def custom_sort(key_value):
    if key_value[1] == 8:
        return (0, key_value[1])  # Помещаем значения 8 в начало
    elif key_value[1] == 7:
        return (1, key_value[1])  # Затем значения 7
    elif key_value[1] == 16:
        return (2, key_value[1])  # Затем значения 16
    else:
        return (3, key_value[1])  # Все остальные значения в конце


def valence(coordinates_current):
    valence_atom = {}
    for i, row in coordinates_current.iterrows():
        if row['element'] == 6:
            valence_atom[row['number_atoms']] = 4
        elif row['element'] == 8:
            valence_atom[row['number_atoms']] = 2
        elif row['element'] == 7:
            valence_atom[row['number_atoms']] = 3
        elif row['element'] == 16:
            valence_atom[row['number_atoms']] = 2
        else:
            valence_atom[row['number_atoms']] = 1
    return valence_atom


def type_connections(connections_current, coordinates_current, cycle, arom_cycle):
    con_element = connections_current.merge(coordinates_current, left_on='atom_1', right_on='number_atoms', how='left')[
        ['atom_1', 'element', 'atom_2', 'length']]
    con_element = con_element.merge(coordinates_current, left_on='atom_2', right_on='number_atoms', how='left')[
        ['atom_1', 'element_x', 'atom_2', 'element_y', 'length']]

    for i, row in con_element.iterrows():
        length = float(row['length'])
        flag_type = False
        flag_type2 = False
        # Если в молекуле есть ароматические циклы, проверим, есть ли такая связь внутри цикла:
        for i_cycle in arom_cycle[arom_cycle['arom'] == True]['number_cycle']:
            number_atom_in_cycle = list(cycle[cycle['number_cycle'] == i_cycle]['atom_in_cycle'])
            number_atom_in_cycle = list(map(int, number_atom_in_cycle))
            if (int(row['atom_1']) in number_atom_in_cycle) and (int(row['atom_2']) in number_atom_in_cycle):
                flag_type = True
                break
        if not flag_type:
            # Углерод-углерод
            if (row['element_x'] == 6) and (row['element_y'] == 6):
                if length > 1.4:
                    connections_current.loc[i, 'type'] = 1
                elif (length > 1.28) and (length < 1.4):
                    connections_current.loc[i, 'type'] = 2
                else:
                    connections_current.loc[i, 'type'] = 3
            # Углерод-кислород
            elif ((row['element_x'] == 6) and (row['element_y'] == 8)) or (
                    (row['element_y'] == 6) and (row['element_x'] == 8)):
                if length > 1.28:
                    connections_current.loc[i, 'type'] = 1
                else:
                    connections_current.loc[i, 'type'] = 2
            # Углерод-азот
            elif ((row['element_x'] == 7) and (row['element_y'] == 6)) or (
                    (row['element_y'] == 7) and (row['element_x'] == 6)):
                if length >= 1.335:
                    connections_current.loc[i, 'type'] = 1
                elif (length > 1.19) and (length < 1.4):
                    connections_current.loc[i, 'type'] = 2
                else:
                    connections_current.loc[i, 'type'] = 3
                    # Азот-кислород
            elif ((row['element_x'] == 7) and (row['element_y'] == 8)) or (
                    (row['element_y'] == 7) and (row['element_x'] == 8)):
                if length > 1.33:
                    connections_current.loc[i, 'type'] = 1
                else:
                    connections_current.loc[i, 'type'] = 2
            # Углерод-сера
            elif ((row['element_x'] == 6) and (row['element_y'] == 16)) or (
                    (row['element_y'] == 6) and (row['element_x'] == 16)):
                if length > 1.7:
                    connections_current.loc[i, 'type'] = 1
                else:
                    connections_current.loc[i, 'type'] = 2
                    # Азот-азот
            elif ((row['element_x'] == 7) and (row['element_y'] == 7)):
                if length > 1.4:
                    connections_current.loc[i, 'type'] = 1
                elif length > 1.15 and length < 1.4:
                    connections_current.loc[i, 'type'] = 2
                else:
                    connections_current.loc[i, 'type'] = 2
            # Азот-кислород
            elif ((row['element_x'] == 7) and (row['element_y'] == 8)) or (
                    (row['element_y'] == 7) and (row['element_x'] == 8)):
                if length > 1.4:
                    connections_current.loc[i, 'type'] = 1
                else:
                    connections_current.loc[i, 'type'] = 2
            # Кислород-кислород
            elif ((row['element_x'] == 8) and (row['element_y'] == 8)):
                if length > 1.35:
                    connections_current.loc[i, 'type'] = 1
                else:
                    connections_current.loc[i, 'type'] = 2
            # Сера-сера
            elif ((row['element_x'] == 16) and (row['element_y'] == 16)):
                if length > 2:
                    connections_current.loc[i, 'type'] = 1
                else:
                    connections_current.loc[i, 'type'] = 2
                    # Остальные связи - одинарные
            else:
                connections_current.loc[i, 'type'] = 1

    valence_atom = valence(coordinates_current)

    # Создаем матрицу связей
    graph = new_graph(coordinates_current, connections_current)
    # Проверяем все ли связи у кислородов:
    num_O = list(coordinates_current[coordinates_current['element'] == 8]['number_atoms'])
    for O in num_O:
        if sum(graph[O - 1]) < 2:
            column = np.where(graph[O - 1] == 1)[0][0]
            graph[O - 1][column] = 2
            graph[column][O - 1] = 2

    order_cycle = {}
    for i_cycle in arom_cycle[arom_cycle['arom'] == True]['number_cycle']:
        number_atom_in_cycle = list(cycle[cycle['number_cycle'] == i_cycle]['atom_in_cycle'])
        number_atom_in_cycle = list(map(int, number_atom_in_cycle))
        # Элементы из которых он состоит
        elements_in_cycle = list(
            coordinates_current.loc[coordinates_current['number_atoms'].isin(number_atom_in_cycle)]['element'])
        # Если цикл чисто из углеродов - ставим в type 0
        if all(x == 6 for x in elements_in_cycle):
            order_cycle[i_cycle] = True
        else:
            order_cycle[i_cycle] = False
    sorted_dict = dict(sorted(order_cycle.items(), key=lambda item: item[1], reverse=True))

    # Если в молекуле есть ароматические циклы, заполним связи внутри них:
    for i_cycle, all_C in sorted_dict.items():
        # Найдем номера атомов, из которых состоит цикл
        number_atom_in_cycle = list(cycle[cycle['number_cycle'] == i_cycle]['atom_in_cycle'])
        number_atom_in_cycle = list(map(int, number_atom_in_cycle))
        # Элементы из которых он состоит
        elements_in_cycle = list(
            coordinates_current.loc[coordinates_current['number_atoms'].isin(number_atom_in_cycle)]['element'])

        # Создадим массив с элементами и номерами атомов так, чтобы вначале шли кислород, азот и сера, а после углерод
        element_and_number = {}
        n = 0
        for i in number_atom_in_cycle:
            element_and_number[i] = elements_in_cycle[n]
            n += 1
        # Массив, в котором элементы идут в порядке, как в цикле
        path_in_cycle = []

        if all_C:
            sum_one = 0
            for num_atom in number_atom_in_cycle:
                if sum_one < sum(graph[num_atom - 1]):
                    start = num_atom
                    sum_one = sum(graph[num_atom - 1])
            path_in_cycle.append(start)
        else:
            atom_not_C = [key for key, value in element_and_number.items() if value != 6]
            path_in_cycle.append(atom_not_C[0])

        while len(path_in_cycle) != len(number_atom_in_cycle):
            for i in np.where(graph[path_in_cycle[-1] - 1] != 0.0)[0]:
                if (i + 1 in number_atom_in_cycle) and (i + 1 not in path_in_cycle):
                    path_in_cycle.append(i + 1)
                    break
        for atom in path_in_cycle:
            # Если у этого атома кол-во связей меньше валентности
            while sum(graph[atom - 1]) < valence_atom[atom]:
                # Находим все индексы элементов, где в матрице у него стоят 1
                index_1 = np.where(graph[atom - 1] == 1.0)[0]
                # Проходимся по каждой связи
                for i_1 in index_1:
                    # Если мы поменяем эту связь на двойную, проверим, пройдет ли тот атом, с которым у него связь проверку на валентность
                    if np.all((sum(graph[i_1]) + 1 <= valence_atom[i_1 + 1]) & (
                            sum(graph[atom - 1]) < valence_atom[atom]) & (i_1 + 1 in number_atom_in_cycle)):
                        graph[atom - 1][i_1] = 2
                        graph[i_1][atom - 1] = 2
                        index_for_type = connections_current[
                            ((connections_current['atom_1'] == atom) & (connections_current['atom_2'] == i_1 + 1)) | (
                                        (connections_current['atom_2'] == atom) & (
                                            connections_current['atom_1'] == i_1 + 1))].index[0]
                        connections_current.iloc[index_for_type, 4] = 2
                    else:
                        index_for_type = connections_current[
                            ((connections_current['atom_1'] == atom) & (connections_current['atom_2'] == i_1 + 1)) | (
                                        (connections_current['atom_2'] == atom) & (
                                            connections_current['atom_1'] == i_1 + 1))].index[0]
                        connections_current.iloc[index_for_type, 4] = 1

        connections_current['type'] = connections_current['type'].fillna(1)
    return connections_current

