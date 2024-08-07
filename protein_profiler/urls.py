from django.urls import path
from . import views


app_name = "protein_profiler"

urlpatterns = [
    path('test_view/', views.test_view, name = "test_view"),
    path('prot_char/', views.prot_char, name = "prot_char"),
    path('submitted_prot_char/', views.submitted_prot_char, name = "submitted_prot_char"),
    path('download_csv/', views.download_csv, name = "download_csv"),
    path('view_protein/<str:id>/', views.view_protein, name='view_protein'),
    path('plot/<str:choice>/', views.plot, name='plot'),
    path('how_to_use/', views.how_to_use, name = 'how_to_use'),
    path('download/', views.download, name = 'download'),
    path('download_zip_file/', views.download_zip_file, name = 'download_zip_file'),
    path('publication/', views.publication, name = 'publication'),
    path('download_plots/', views.download_plots, name = 'download_plots'),
    path('contact/', views.contact, name = 'contact')
]

