from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth import get_user_model

Buyer = get_user_model()

class RegisterForm(UserCreationForm):
    email = forms.EmailField(label=("Email"), required=True)
    first_name = forms.CharField(label=("First Name"), max_length=30, required=True)
    last_name = forms.CharField(label=("Last Name"), max_length=30, required=True)

    class Meta(UserCreationForm.Meta):
        model = Buyer
        fields = ['username', 'email', 'first_name', 'last_name', 'password1', 'password2']
        labels = {
            'username': ('Username'),
            'email': ('Email'),
            'first_name': ('First Name'),
            'last_name': ('Last Name'),
            'password1': ('Password'),
            'password2': ('Confirm Password'),
        }

    def clean_email(self):
        email = self.cleaned_data.get('email')
        if Buyer.objects.filter(email__iexact=email).exists():
            raise forms.ValidationError('Пользователь с таким email уже существует.')
        return email

    def save(self, commit=True):
        user = super().save(commit=False)
        user.email = self.cleaned_data['email']
        user.first_name = self.cleaned_data['first_name']
        user.last_name = self.cleaned_data['last_name']
        if commit:
            user.save()
        return user
