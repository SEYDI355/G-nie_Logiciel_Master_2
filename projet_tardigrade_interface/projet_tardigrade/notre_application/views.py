from django.contrib import messages
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, redirect, render

from .forms import DatasetForm
from .models import Dataset
from .services.kegg_network import build_dataset_network_payload
from .services.metabolomics import analyze_and_store_dataset


def accueil(request):
    """Page d'accueil du site."""
    return render(request, "accueil.html")


def home(request):
    """Liste des jeux de données + formulaire d'import."""
    if request.method == "POST":
        form = DatasetForm(request.POST, request.FILES)
        if form.is_valid():
            ds = form.save()
            try:
                measurements = analyze_and_store_dataset(ds)
                mapped = sum(1 for m in measurements if m.compound_id and ";" not in m.compound_id)
                messages.success(
                    request,
                    f"Jeu de données importé : {len(measurements)} comparaisons analysées "
                    f"({mapped} associées à un identifiant KEGG).",
                )
            except ValueError as exc:
                ds.delete()
                messages.error(request, f"Format de fichier invalide : {exc}")
                return render(
                    request,
                    "notre_application/home.html",
                    {"form": form, "datasets": Dataset.objects.all()},
                )
            return redirect("notre_application:detail", pk=ds.pk)
    else:
        form = DatasetForm()
    return render(request, "notre_application/home.html", {"form": form, "datasets": Dataset.objects.all()})


def dataset_detail(request, pk):
    """Page d'analyse d'un jeu de données."""
    ds = get_object_or_404(Dataset, pk=pk)
    return render(request, "notre_application/detail.html", {"dataset": ds})


def dataset_network_json(request, pk):
    """Retourne les mesures métabolomiques mappées sur KEGG pour Cytoscape.js."""
    ds = get_object_or_404(Dataset, pk=pk)
    condition = request.GET.get("condition")
    try:
        payload = build_dataset_network_payload(ds, condition=condition)
    except ValueError as exc:
        return JsonResponse({"error": str(exc)}, status=400)
    return JsonResponse(payload)
