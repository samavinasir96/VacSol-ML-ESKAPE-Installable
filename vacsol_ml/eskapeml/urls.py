from django.urls import path
from . import views

app_name = "eskapeml"

urlpatterns = [
    path('', views.home, name='home'),
    path('upload/', views.upload_sequence, name='upload_sequence'),
    path('processing/', views.processing_view, name='processing'),
    path('faqs/', views.faqs, name='faqs'),
    path('glossary/', views.glossary, name='glossary'),
    path('results/', views.get_results, name='results'),
    path("get_progress/", views.get_progress, name="get_progress"),
    path("download-results/", views.eskapeml_download_csv, name="eskapeml_download_csv"),
    path("download/biological/", views.download_biological_csv, name="download_biological"),
    path("download/physiochemical/", views.download_physiochemical_csv, name="download_physiochemical"),
    path("track_citation/", views.track_citation, name="track_citation"),
    path('clear_temp/', views.clear_temp, name='clear_temp'),
]
