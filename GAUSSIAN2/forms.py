from django import forms
from .models import Get

class CASForm(forms.ModelForm):

    class Meta:
        model = Get
        fields = ('text')