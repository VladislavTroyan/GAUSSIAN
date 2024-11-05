from django.db import models
from django.contrib.auth.models import User
from django.utils import timezone

class UploadedFile(models.Model):
    # Поле id, которое будет связываться с id в таблице Raschety
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    file = models.FileField(upload_to='uploaded_files/', max_length=255)
    original_filename = models.CharField(max_length=999, unique=False)  # Уникальность имени файла
    uploaded_at = models.DateTimeField(default=timezone.now)
    processed_successfully = models.BooleanField(default=False)
    error_message = models.TextField(null=True, blank=True)


    def __str__(self):
        return self.original_filename

class Error(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    filename = models.CharField(max_length=255)
    error_type = models.CharField(max_length=255)
    error_message = models.TextField()
    created_at = models.DateTimeField(default=timezone.now)
    processed_successfully = models.BooleanField(default=False)

    def __str__(self):
        return f"{self.filename} - {self.error_type}"
