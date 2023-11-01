from django.urls import path
from . import views


app_name = "protein_profiler"

urlpatterns = [
    path('test_view/', views.test_view, name = "test_view"),
    path('prot_char/', views.prot_char, name = "prot_char"),
    path('submitted_prot_char/', views.submitted_prot_char, name = "submitted_prot_char"),
    path('download_csv/', views.download_csv, name = "download_csv"),
    path('view_protein/<str:id>/', views.view_protein, name='view_protein'),
    path('plot/<str:choice>/', views.plot, name='plot')
    
]

