from django.contrib.auth.views import PasswordChangeView, PasswordChangeDoneView, PasswordResetConfirmView
from django.urls import path
from .import views
from .views import proba, start, homepage, search, check, profile_view, RegisterView, CustomPasswordResetView, \
    CustomPasswordResetDoneView, CustomPasswordResetCompleteView, molecule_viewer
from django.contrib.auth import views as auth_views
from django.conf.urls.static import static
from django.conf import settings
from .views import download_specific_file


urlpatterns = [
    path('', views.proba),
    path('start/', views.start),
    path('homepage/', views.homepage, name='homepage'),
    path('check/', views.check, name='check'),
    path('search_results/', search, name='search_results'),
    path('profile/', profile_view, name='profile'),
    path('register/', RegisterView.as_view(), name='register'),
    path('password_reset/', CustomPasswordResetView.as_view(success_url='/password_reset/done/'), name='custom_password_reset'),
    path('password_change/', PasswordChangeView.as_view(template_name="registration/custom_password_change_form.html",success_url='/password_change/done/'), name='custom_password_change'),
    path('password_change/done/', PasswordChangeDoneView.as_view(template_name='registration/custom_password_change_done.html'), name='custom_password_change_done'),
    path('password_reset/done/', CustomPasswordResetDoneView.as_view(template_name = "registration/custom_password_reset_done.html"), name='password_reset_done'),
    path('accounts/reset/<uidb64>/<token>/', PasswordResetConfirmView.as_view(template_name='registration/custom_set_password.html'), name='password_reset_confirm'),
    path('accounts/reset/done/', CustomPasswordResetCompleteView.as_view(template_name = "registration/custom_password_reset_complete.html"), name='password_done'),
    path('logout/', auth_views.LogoutView.as_view(), name='logout'),
    path('molecule-viewer/<int:value>/', molecule_viewer, name='molecule_viewer'),
    path('download-specific-file/<int:raschety_id>/', download_specific_file, name='download_specific_file'),
]
if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)