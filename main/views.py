import mimetypes
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.shortcuts import render, redirect
from django.db import connections
from django.views import View
from main.forms import UserCreationForm
from django import forms
from main.forms import RegisterForm
from django.http import JsonResponse
from django.contrib.auth.views import PasswordResetView, PasswordResetDoneView, PasswordResetCompleteView, \
    PasswordChangeDoneView
from django.conf import settings
from django.shortcuts import render
from .script import create_molecule_image, get_unique_elements, get_groups_data, get_atom_color
from django.http import HttpResponse
import os
from django.http import HttpResponseNotFound, HttpResponseServerError
from django.db import connection
def proba(request):
    return render(request, 'main/about.html')

def start(request):
    return render(request, 'main/start.html')

@login_required
def profile_view(request):
    user_email = request.user.email
    return render(request, 'web/profile.html', {'user_email': user_email})


class CustomPasswordResetView(PasswordResetView):
    template_name = 'registration/custom_password_reset_form.html'

class CustomPasswordResetDoneView(PasswordResetDoneView):
    template_name = "registration/custom_password_reset_done.html"


class CustomPasswordResetCompleteView(PasswordResetCompleteView):
    template_name = "registration/custom_password_reset_complete.html"


class CustomPasswordChangeDoneView(PasswordChangeDoneView):
    template_name = "registration/custom_password_change_done.html"


class CustomUserCreationForm(UserCreationForm):
    email = forms.EmailField()

class RegisterView(View):
    template_name = 'registration/register.html'

    def get(self, request):
        context = {
            'form': RegisterForm()
        }
        return render(request, self.template_name, context)

    def post(self, request):
        form = RegisterForm(request.POST)

        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password1')
            user = authenticate(username=username, password=password)
            login(request, user)
            return redirect('profile')

        context = {
            'form': form
        }
        return render(request, self.template_name, context)

@login_required()
def download_specific_file(request, raschety_id):
    try:
        # Выполняем SQL-запрос для получения пути к файлу и его имени по raschety_id
        with connection.cursor() as cursor:
            cursor.execute("""
                SELECT i.file, i.original_filename
                FROM raschety r
                INNER JOIN import_uploadedfile i ON r.id = i.id
                WHERE r.id = %s
            """, [raschety_id])
            row = cursor.fetchone()

        if row:
            file_relative_path = row[0]  # Теперь только относительный путь к файлу
            original_filename = row[1]  # Получаем оригинальное имя файла

            # Полный путь к файлу в папке uploaded_files
            file_path = os.path.join(settings.MEDIA_ROOT, file_relative_path)

            if os.path.exists(file_path):
                # Определяем тип содержимого
                content_type, _ = mimetypes.guess_type(file_path)
                if content_type is None:
                    content_type = 'application/octet-stream'

                # Создаем HTTP-ответ с содержимым файла
                with open(file_path, 'rb') as file:
                    response = HttpResponse(file.read(), content_type=content_type)
                # Устанавливаем заголовок Content-Disposition с оригинальным именем файла
                response['Content-Disposition'] = f'attachment; filename="{original_filename}"'

                return response
            else:
                return HttpResponseNotFound('Файл не найден')
        else:
            return HttpResponseNotFound('Расчёт не найден')

    except Exception as e:
        return HttpResponseServerError(f'Произошла ошибка: {str(e)}')


def search(request):
    if request.method == 'GET':
        userId_param = request.GET.get('userId')
        userId = 0 if userId_param is None or userId_param == 'NaN' else int(userId_param)
        search_query = request.GET.get('search_query')

        # Разбиваем строку поиска на отдельные слова
        keywords = search_query.split()

        # Создаем строки SQL-запросов
        favorites = f"""
            select favorites_cas from favorites
            where favorites_cas = %s and id_user = {userId}
        """

        molecule_id = """
            select id from raschety
            where cas = %s
        """

        sql_query_molecules = """
            SELECT *
            FROM molecules
            WHERE cas = %s
        """

        sql_query_g4 = """
            SELECT raschety.id, date, temperature, pressure, energy_gibbs, enthalpy, dfh,  at, entropy
            FROM raschety, molecules, methods, g4
            WHERE raschety.cas = molecules.cas AND raschety.method = methods.method AND raschety.id = g4.id AND raschety.cas = %s
            order by date desc
        """

        sql_query_g3mp2 = """
            SELECT raschety.id, date, temperature, pressure, energy_gibbs, enthalpy, dfh,  at, entropy
            FROM raschety, molecules, methods, g3mp2
            WHERE raschety.cas = molecules.cas AND raschety.method = methods.method AND raschety.id = g3mp2.id AND raschety.cas = %s
            order by date desc

        """

        sql_query_b3lyp = """
            SELECT raschety.id, date, temperature, pressure, energy_gibbs, correction_enthalpy, correction_energy, dfh
            FROM raschety, molecules, methods, b3lyp
            WHERE raschety.cas = molecules.cas AND raschety.method = methods.method AND raschety.id = b3lyp.id AND raschety.cas = %s
            order by date desc

        """

        sql_query_m06 = """
            SELECT raschety.id, date, temperature, pressure, energy_gibbs, correction_enthalpy, correction_energy, dfh
            FROM raschety, molecules, methods, m06
            WHERE raschety.cas = molecules.cas AND raschety.method = methods.method AND raschety.id = m06.id AND raschety.cas = %s
            order by date desc

        """

        # Выполняем запросы к базе данных
        with connections['default'].cursor() as cursor:

            cursor.execute(favorites, [search_query])
            favorites = cursor.fetchall()

            cursor.execute(molecule_id, [search_query])
            molecule_id = cursor.fetchall()

            cursor.execute(sql_query_molecules, [search_query])
            results_molecules = cursor.fetchall()

            cursor.execute(sql_query_g4, [search_query])
            results_g4 = cursor.fetchall()

            cursor.execute(sql_query_g3mp2, [search_query])
            results_g3mp2 = cursor.fetchall()

            cursor.execute(sql_query_b3lyp, [search_query])
            results_b3lyp = cursor.fetchall()

            cursor.execute(sql_query_m06, [search_query])
            results_m06 = cursor.fetchall()

        # Подготавливаем контекст для передачи в шаблон
        data = {
            'molecule_id': molecule_id,
            'results_molecules': results_molecules,
            'results_g4': results_g4,
            'results_g3mp2': results_g3mp2,
            'results_b3lyp': results_b3lyp,
            'results_m06': results_m06,
            'search_query': search_query,
            'favorites': favorites,
        }
        return JsonResponse(data)


def homepage(request):
    if request.method == 'POST':
        userId_param = request.POST.get('userId')
        userId = 0 if userId_param is None or userId_param == 'NaN' else int(userId_param)
        checkbox_values = []
        for i in range(0, 18):
            checkbox_value = request.POST.get(str(i), 'off')
            if checkbox_value == 'on':
                checkbox_values.append(str(i))

        # Проверяем, содержится ли цифра '0' в списке
        if '0' in checkbox_values:
            flagfavorites = 1
            checkbox_values.remove('0')
            #print("Цифра '0' содержится в списке.")
        else:
            flagfavorites = 0
            #print("Цифра '0' не содержится в списке.")

        F = int(request.POST.get('F', 0) or 0)
        Cl = int(request.POST.get('Cl', 0) or 0)
        O = int(request.POST.get('O', 0) or 0)
        C = int(request.POST.get('C', 0) or 0)
        Br = int(request.POST.get('Br', 0) or 0)
        H = int(request.POST.get('H', 0) or 0)
        N = int(request.POST.get('N', 0) or 0)
        S = int(request.POST.get('S', 0) or 0)
        I = int(request.POST.get('I', 0) or 0)

        # Формируем строку fg из значений чекбоксов
        fg = ','.join(checkbox_values)
        list_count = {'F': F, 'Cl': Cl, 'O': O, 'C': C, 'Br': Br, 'H': H, 'N': N, 'S': S, 'I': I}

        # Формируем строку для запроса
        string_to_sql = ''

        for key, value in list_count.items():
            if value == 1:
                string_to_sql += f'''formula ~ '(^|[^A-Z]){key}([^0-9l]|$)' AND '''
            elif value > 1:
                string_to_sql += f'''formula ~ '(^|[^A-Z]){key}{value}([^0-9]|$)' AND '''

        # Удалим последний and из строки
        last_and_index = string_to_sql.rfind(" AND ")
        if last_and_index != -1:
            string_to_sql = string_to_sql[:last_and_index]

        # Печать строки для отладки
        print(string_to_sql)

        # Проверяем, есть ли вообще необходимость формирования запроса. Если все нули - выводим все CAS
        all_zero = all(value == 0 for value in list_count.values())



        if not fg and all_zero and not flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы не выбраны и молекулы не выбраны
            sql_query = f'''
                SELECT cas
                FROM molecules
                ORDER BY cas ASC
                LIMIT 28 OFFSET {value}
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(sql_query)
                cas = cursor.fetchall()
            cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
            if not cas_values:
                cas_values = "'None_Values'"
            print("Чекбоксы не выбраны, молекулы не выбраны, не избранное", cas_values)

        if not fg and all_zero and flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы не выбраны и молекулы не выбраны
            sql_query = f'''
                SELECT m.cas
                FROM molecules AS m
                JOIN favorites AS f ON m.cas = f.favorites_cas
                WHERE f.id_user = {userId}
                ORDER BY m.cas ASC
                LIMIT 28 OFFSET {value}
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(sql_query)
                cas = cursor.fetchall()
            cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
            if not cas_values:
                cas_values = "'None_Values'"
            #print("Чекбоксы не выбраны, молекулы не выбраны, избранное", cas_values)


        elif fg and all_zero and not flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы выбраны и молекулы не выбраны
            sql_query = f'''
                SELECT b.cas
                FROM (
                    SELECT cas, count(num_fg) AS count_fg
                    FROM fg_molecules
                    WHERE num_fg IN ({fg})
                    GROUP BY cas
                ) b,
                (
                    SELECT LENGTH('{fg}') - LENGTH(REPLACE('{fg}', ',', '')) + 1 AS count_query
                ) a
                WHERE b.count_fg = a.count_query
                GROUP BY cas
                ORDER BY cas ASC
                LIMIT 28 OFFSET {value}
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(sql_query)
                cas = cursor.fetchall()
            cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
            if not cas_values:
                cas_values = "'None_Values'"
            #print("Чекбоксы выбраны, молекулы не выбраны, не избранное",cas_values)



        elif fg and all_zero and flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы выбраны и молекулы не выбраны
            sql_query = f'''
                SELECT b.cas
                FROM (
                    SELECT cas, count(num_fg) AS count_fg
                    FROM fg_molecules
                    WHERE num_fg IN ({fg})
                    GROUP BY cas
                ) b
                JOIN (
                    SELECT LENGTH('{fg}') - LENGTH(REPLACE('{fg}', ',', '')) + 1 AS count_query
                ) a ON b.count_fg = a.count_query
                JOIN favorites AS f ON b.cas = f.favorites_cas
                WHERE f.id_user = {userId} 
                GROUP BY b.cas
                ORDER BY b.cas ASC
                LIMIT 28 OFFSET {value} 
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(sql_query)
                cas = cursor.fetchall()
            cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
            if not cas_values:
                cas_values = "'None_Values'"
            #print("Чекбоксы выбраны, молекулы не выбраны, избранное",cas_values)



        elif not fg and not all_zero and not flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы не выбраны и молекулы выбраны
            resultat = f'''
                select DISTINCT molecules.cas 
                from molecules 
                where (molecules.cas IN (SELECT cas FROM molecules ORDER BY cas ASC)) 
                AND ({string_to_sql}) 
                ORDER BY cas ASC
                LIMIT 28 OFFSET {value}
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(resultat)
                cas = cursor.fetchall()
                cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
                if not cas_values:
                    cas_values = "'None_Values'"
                #print("Чекбоксы не выбраны, молекулы выбраны, не избранное", cas_values)



        elif not fg and not all_zero and flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы не выбраны и молекулы выбраны
            resultat = f'''
                SELECT DISTINCT m.cas
                FROM molecules m
                JOIN favorites f ON m.cas = f.favorites_cas
                WHERE f.id_user = {userId}
                  AND ({string_to_sql})
                ORDER BY m.cas ASC
                LIMIT 28 OFFSET {value}
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(resultat)
                cas = cursor.fetchall()
                cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
                if not cas_values:
                    cas_values = "'None_Values'"
                #print("Чекбоксы не выбраны, молекулы выбраны, избранное", cas_values)


        elif fg and not all_zero and not flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы выбраны и молекулы выбраны
            resultat = f'''
                select DISTINCT molecules.cas 
                from molecules 
                where (molecules.cas IN (SELECT b.cas
                FROM (
                    SELECT cas, count(num_fg) AS count_fg
                    FROM fg_molecules
                    WHERE num_fg IN ({fg})
                    GROUP BY cas
                ) b,
                (
                    SELECT LENGTH('{fg}') - LENGTH(REPLACE('{fg}', ',', '')) + 1 AS count_query
                ) a
                WHERE b.count_fg = a.count_query
                GROUP BY cas
                ORDER BY cas ASC)) 
                AND ({string_to_sql}) ORDER BY cas ASC
                LIMIT 28 OFFSET {value}
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(resultat)
                cas = cursor.fetchall()
                cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
                if not cas_values:
                    cas_values = "'None_Values'"
                #print("Чекбоксы выбраны, молекулы выбраны, не избранное", cas_values)


        elif fg and not all_zero and flagfavorites:
            value = int(request.POST.get('value', 0))
            # Чекбоксы выбраны и молекулы выбраны
            resultat = f'''
                SELECT DISTINCT m.cas
                FROM molecules m
                JOIN favorites f ON m.cas = f.favorites_cas
                WHERE f.id_user = {userId}
                  AND m.cas IN (
                    SELECT b.cas
                    FROM (
                        SELECT cas, count(num_fg) AS count_fg
                        FROM fg_molecules
                        WHERE num_fg IN ({fg})
                        GROUP BY cas
                    ) b,
                    (
                        SELECT LENGTH('{fg}') - LENGTH(REPLACE('{fg}', ',', '')) + 1 AS count_query
                    ) a
                    WHERE b.count_fg = a.count_query
                    GROUP BY cas
                    ORDER BY cas ASC
                  )
                  AND ({string_to_sql})
                ORDER BY m.cas ASC
                LIMIT 28 OFFSET {value}
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(resultat)
                cas = cursor.fetchall()
                cas_values = ', '.join([f"'{cas[0]}'" for cas in cas])
                if not cas_values:
                    cas_values = "'None_Values'"
                #print("Чекбоксы выбраны, молекулы выбраны, избранное", cas_values)



        favorites = f'''

            WITH cas_values_table AS (
                SELECT unnest(ARRAY[{cas_values}]) AS cas_value,
                       generate_subscripts(ARRAY[{cas_values}], 2) AS array_index
            )
            SELECT f.favorites_cas
            FROM cas_values_table
            LEFT JOIN favorites AS f 
                ON cas_values_table.cas_value = f.favorites_cas 
                AND f.id_user = {userId}
            ORDER BY cas_values_table.array_index;
        '''
        name = f'''
            SELECT name_molecules
            FROM molecules
            WHERE cas IN ({cas_values})
            ORDER BY cas ASC
        '''
        formula = f'''
            SELECT formula
            FROM molecules
            WHERE cas IN ({cas_values})
            ORDER BY cas ASC
        '''
        metody = f'''
            SELECT 
                COALESCE(
                    STRING_AGG(DISTINCT methods.name_method, ', '), 
                    'Расчёты не найдены'
                ) AS methods_list
            FROM methods
            INNER JOIN raschety ON raschety.method = methods.method
            RIGHT JOIN molecules ON molecules.cas = raschety.cas
            WHERE molecules.cas IN ({cas_values})
            GROUP BY molecules.cas
            ORDER BY molecules.cas ASC
        '''


        with connections['default'].cursor() as cursor:
            cursor.execute(favorites)
            favorites = cursor.fetchall()
        with connections['default'].cursor() as cursor:
            cursor.execute(name)
            name = cursor.fetchall()
        with connections['default'].cursor() as cursor:
            cursor.execute(formula)
            formula = cursor.fetchall()
        with connections['default'].cursor() as cursor:
            cursor.execute(metody)
            metody = cursor.fetchall()

        # Возвращаем результаты в формате JSON
        return JsonResponse({
            'results': [
                {'type': 'favorites', 'value': favorites},
                {'type': 'name', 'value': name},
                {'type': 'formula', 'value': formula},
                {'type': 'cas', 'value': cas},
                {'type': 'metody', 'value': metody},
            ],
        })
    # Если GET-запрос, просто отображаем страницу
    return render(request, 'main/homepage.html')

def check(request):
    if request.method == 'POST':
        user_id = request.POST.get('userId', None)
        checkbox_name = request.POST.get('name', None)
        checkbox_state = request.POST.get('state', None)
        if checkbox_state == "1":
            INSERT = f'''
            INSERT INTO favorites (id_user, favorites_cas) 
            VALUES ({user_id}, '{checkbox_name}');
            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(INSERT)
            #print("Произошла вставка")
        elif checkbox_state == "2":
            DELETE = f'''
            DELETE FROM favorites
            WHERE id_user = {user_id} 
            AND favorites_cas = '{checkbox_name}';

            '''
            with connections['default'].cursor() as cursor:
                cursor.execute(DELETE)
            #print("Произошло удаление")
        else:
            print("ОШИБКА ЧЕКБОКСОВ")
        return JsonResponse({})


def molecule_viewer(request, value):
    # Получаем конфигурацию базы данных из настроек Django
    connection_params = {
        'user': settings.DATABASES['default']['USER'],
        'password': settings.DATABASES['default']['PASSWORD'],
        'host': settings.DATABASES['default']['HOST'],
        'port': settings.DATABASES['default']['PORT'],
        'database': settings.DATABASES['default']['NAME'],
    }

    molecule_data = create_molecule_image(connection_params, value)
    unique_elements = get_unique_elements(connection_params, value)
    groups_data = get_groups_data(connection_params, value)
    return render(request, 'main/molecule_viewer.html', {
        'molecule_data': molecule_data,
        'unique_elements': unique_elements,
        'groups_data': groups_data,

    })
