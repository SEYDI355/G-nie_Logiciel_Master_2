from django.urls import path

from . import views

app_name = "notre_application"

urlpatterns = [
    path("", views.accueil, name="accueil"),
    path("accueil.html", views.accueil, name="accueil_html"),
    path("datasets/", views.home, name="home"),
    path("datasets/<int:pk>/", views.dataset_detail, name="detail"),
    path("datasets/<int:pk>/network.json", views.dataset_network_json, name="network_json"),
]
