# import.views - это ваш модуль с представлениями. Пожалуйста, убедитесь, что путь к нему верный.

from django.urls import path
from . import views
from .views import upload_files, profile

urlpatterns = [
    path('upload-files/', upload_files, name='upload_files'),
    path('profile/', profile, name='profile'),
]
