from django.db import models

# Create your models here.
class Dataset(models.Model):
    """Un jeu de données de métabolomique."""

    name = models.CharField("Nom", max_length=120)
    description = models.TextField("Description", blank=True)
    matrix_file = models.FileField("Matrice métabolites × échantillons", upload_to="datasets/")
    metadata_file = models.FileField("Métadonnées des échantillons", upload_to="datasets/")
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return self.name