from django.urls import path
from . import views


app_name = "protein_profiler"

urlpatterns = [
    path('test_view/', views.test_view, name = "test_view"),
    path('prot_char/', views.prot_char, name = "prot_char"),
    path('submitted_prot_char/', views.submitted_prot_char, name = "submitted_prot_char")
    
]

