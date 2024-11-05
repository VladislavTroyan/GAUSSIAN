import os
import hashlib
import psycopg2
from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.core.exceptions import ValidationError
from django.views.decorators.csrf import csrf_protect
from django.db import connection
from .models import UploadedFile, Error
from .file_function import files_processing
from GAUSSIAN2 import settings
from django.core.paginator import Paginator
from rdkit import Chem
from rdkit.Chem import Draw


@csrf_protect
@login_required
def profile(request):
    user = request.user

    # Получаем историю загрузок и ошибок пользователя
    upload_error_history = []

    # Получаем историю загрузок пользователя из UploadedFile
    for uploaded_file in UploadedFile.objects.filter(user=user).order_by('-uploaded_at'):
        raschety_id = uploaded_file.id  # Получаем id из UploadedFile
        cas = None
        name_method = None
        processed_successfully = uploaded_file.processed_successfully

        # Получаем значение cas и name_method из таблиц raschety и methods соответственно, используя JOIN
        with connection.cursor() as cursor:
            cursor.execute('''
                SELECT r.cas, m.name_method
                FROM raschety r
                INNER JOIN methods m ON r.method = m.method
                WHERE r.id = %s;
            ''', [raschety_id])
            raschety_data = cursor.fetchone()

            # Если есть данные о raschety, извлекаем cas и name_method
            if raschety_data:
                cas, name_method = raschety_data

        # Добавляем информацию о загрузке в список upload_error_history
        upload_error_history.append({
            'type': 'upload',
            'original_filename': uploaded_file.original_filename,
            'uploaded_at': uploaded_file.uploaded_at,
            'cas': cas,
            'name_method': name_method,
            'processed_successfully': processed_successfully
        })

    # Получаем историю ошибок пользователя из Error
    for error in Error.objects.filter(user=user).order_by('-created_at'):
        upload_error_history.append({
            'type': 'error',
            'filename': error.filename,
            'error_type': error.error_type,
            'error_message': error.error_message,
            'created_at': error.created_at,
            'processed_successfully': error.processed_successfully
        })

    # Сортируем объединенный список по дате загрузки или создания, чтобы последние события были вверху
    upload_error_history.sort(key=lambda x: x.get('uploaded_at') or x.get('created_at'), reverse=True)

    # Инициализируем Paginator
    paginator = Paginator(upload_error_history, 10)  # Разбиваем на страницы по 10 элементов

    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)

    return render(request, 'web/profile.html', {'user': user, 'page_obj': page_obj})


@csrf_protect
@login_required
def upload_files(request):
    if request.method == 'POST':
        user = request.user
        uploaded_files = request.FILES.getlist('fileInput')
        errors = []
        success_files = []

        for uploaded_file in uploaded_files:
            try:
                if not uploaded_file.name.lower().endswith('.log'):
                    raise ValidationError('Недопустимый формат файла. Разрешены только файлы с расширением .log')

                print(f"Processing file: {uploaded_file.name}")

                error_info = files_processing(uploaded_file, user_db=user.username, password_db=user.password,
                                              name_database='postgres')
                if error_info:
                    error_type = error_info.get('error_type')
                    filename = error_info.get('filename')

                    # Определяем сообщение об ошибке на основе error_type
                    error_messages = {
                        1: 'Незавершённый расчёт.',
                        2: 'Нет CAS.',
                        3: 'Неизвестный метод.',
                        4: 'Расчёт уже есть в базе.',
                        5: 'Неизвестная ошибка. Пожалуйста, повторите попытку позже или обратитесь к администратору.'
                    }
                    error_message = error_messages.get(error_type, 'Неизвестная ошибка')

                    # Сохраняем информацию об ошибке в модель Error
                    Error.objects.create(
                        user=user,
                        filename=filename,
                        error_type=error_type,
                        error_message=error_message
                    )
                    errors.append({'file_name': filename, 'error_type': error_type})
                else:
                    # Если нет ошибки, добавляем файл в список успешно загруженных
                    success_files.append({'file_name': uploaded_file.name})

                    # Вычисляем хэш файла
                    hash_object = hashlib.sha1()
                    for chunk in uploaded_file.chunks():
                        hash_object.update(chunk)
                    hash_str = hash_object.hexdigest()

                    # Генерируем уникальное имя файла на основе хэша
                    filename = f"{uploaded_file.name}_{hash_str}"
                    file_path = os.path.join(settings.MEDIA_ROOT, 'uploaded_files', filename)

                    with open(file_path, 'wb') as file:
                        for chunk in uploaded_file.chunks():
                            file.write(chunk)

                    with connection.cursor() as cursor:
                        cursor.execute('SELECT id FROM raschety ORDER BY id DESC LIMIT 1')
                        raschety_id = cursor.fetchone()[0]

                        cursor.execute('''
                            SELECT raschety.id, molecules.name_molecules, methods.name_method
                            FROM molecules
                            INNER JOIN raschety ON molecules.cas = raschety.cas
                            INNER JOIN methods ON raschety.method = methods.method
                            WHERE raschety.id = %s;
                        ''', [raschety_id])
                        result = cursor.fetchone()
                        if result:
                            raschety_id, name_molecules, name_method = result
                            original_filename = f"{raschety_id}_{name_molecules}_{name_method}.log"
                        else:
                            raise Exception("Values not found for the specified raschety ID.")

                    UploadedFile.objects.create(
                        id=raschety_id,
                        user=user,
                        file=file_path,
                        original_filename=original_filename,
                        processed_successfully=True
                    )

                    # Генерация изображения молекулы и сохранение его в media/photos
                    print(f"Generating molecule image for raschety_id: {raschety_id}")
                    generate_molecule_image(raschety_id)

            except ValidationError as ve:
                errors.append({'file_name': uploaded_file.name, 'error_message': ve.message})
                print(f"Error processing file {uploaded_file.name}: {ve.message}")

            except Exception as e:
                errors.append({'file_name': uploaded_file.name, 'error_message': str(e)})
                print(f"Error processing file {uploaded_file.name}: {e}")
                print(f"Error type: {type(e)}")

        # Передаем ошибки и успешно загруженные файлы в шаблон профиля пользователя для их отображения
        return JsonResponse({'errors': errors, 'success_files': success_files})

    return JsonResponse({'message': 'Invalid request method.'}, status=400)


def generate_molecule_image(raschety_id):
    # Подключение к базе данных
    connection_params = {
        'user': "postgres",
        'password': "11111111",
        'host': "127.0.0.1",
        'port': "5432",
        'database': "GAUSSIAN"
    }

    conn = psycopg2.connect(**connection_params)
    cur = conn.cursor()

    # Запрос к таблице elements для получения символов элементов
    cur.execute("SELECT element, symbol FROM elements")
    elements_data = cur.fetchall()

    # Получение cas для данного raschety_id
    cur.execute("SELECT cas FROM raschety WHERE id = %s", (raschety_id,))
    cas_data = cur.fetchone()
    if not cas_data:
        print(f"Error: No CAS found for raschety ID {raschety_id}.")
        return

    cas = cas_data[0]

    # Получение координат атомов
    cur.execute("""
        SELECT c.x, c.y, c.z, c.element 
        FROM coordinates c 
        INNER JOIN raschety r ON c.id = r.id 
        WHERE r.id = %s
    """, (raschety_id,))
    coordinates = cur.fetchall()

    # Создание молекулы с атомами
    mol = Chem.RWMol()
    atom_indices = {}  # Словарь для отображения id атома в индекс в молекуле
    for idx, (_, _, _, element) in enumerate(coordinates):
        # Находим символ элемента по имени в полученных данных
        symbol = next((symbol for _, symbol in elements_data if _ == element), None)
        if symbol:
            atom = Chem.Atom(symbol)
            atom_idx = mol.AddAtom(atom)
            atom_indices[idx + 1] = atom_idx  # Сохраняем соответствие id атома и его индекса в молекуле
        else:
            print(f"Warning: Element symbol not found for atom {element}.")

    # Получение связей для данного id
    cur.execute("SELECT atom_1, atom_2, length, type FROM connections WHERE id = %s", (raschety_id,))
    connections = cur.fetchall()

    # Добавление связей между атомами
    for atom1, atom2, _, bond_type in connections:
        bond = Chem.BondType.SINGLE if bond_type == 1 else Chem.BondType.DOUBLE  # Пример: использование одиночной или двойной связи
        mol.AddBond(atom_indices[atom1], atom_indices[atom2], bond)

    # Генерация 2D координат для молекулы
    Chem.rdDepictor.Compute2DCoords(mol)

    # Создание изображения молекулы
    img = Draw.MolToImage(mol, size=(1800, 1200), kekulize=True, wedgeBonds=True)

    # Сохранение изображения в файл
    img_path = os.path.join(settings.MEDIA_ROOT, 'photos', f'{cas}.png')
    img.save(img_path)

    print(f"Image saved for CAS {cas}")

    # Закрытие курсора и соединения с базой данных
    cur.close()
    conn.close()



