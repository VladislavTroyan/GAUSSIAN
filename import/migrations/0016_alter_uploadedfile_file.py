# Generated by Django 5.0.2 on 2024-04-15 06:29

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('import', '0015_alter_uploadedfile_original_filename'),
    ]

    operations = [
        migrations.AlterField(
            model_name='uploadedfile',
            name='file',
            field=models.FileField(max_length=255, upload_to='uploaded_files/'),
        ),
    ]
