# Generated by Django 5.0.2 on 2024-03-17 10:52

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('import', '0007_uploadedfile_cas_uploadedfile_method_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='uploadedfile',
            name='original_filename',
            field=models.CharField(max_length=255),
        ),
    ]