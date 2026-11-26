from django.shortcuts import get_object_or_404, redirect, render

from .forms import DatasetForm
from .models import Dataset


def accueil(request):
    """Page d'accueil du site."""
    return render(request, "accueil.html")


def home(request):
    """Liste des jeux de données + formulaire d'import."""
    if request.method == "POST":
        form = DatasetForm(request.POST, request.FILES)
        if form.is_valid():
            ds = form.save()
            return redirect("notre_application:detail", pk=ds.pk)
    else:
        form = DatasetForm()
    return render(request, "notre_application/home.html", {"form": form, "datasets": Dataset.objects.all()})


def dataset_detail(request, pk):
    """Page d'analyse d'un jeu de données (coque vide : les graphes sont à intégrer)."""
    ds = get_object_or_404(Dataset, pk=pk)
    return render(request, "notre_application/detail.html", {"dataset": ds})
