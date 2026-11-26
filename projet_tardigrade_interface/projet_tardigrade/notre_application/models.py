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


class MetaboliteMeasurement(models.Model):
    """Mesure statistique et mapping KEGG d'un métabolite pour un jeu de données."""

    dataset = models.ForeignKey(
        Dataset,
        related_name="measurements",
        on_delete=models.CASCADE,
    )
    metabolite_name = models.CharField("Nom du métabolite", max_length=200)
    pathway_label = models.CharField("Voie métabolique (CSV)", max_length=160, blank=True)
    compound_id = models.CharField("Identifiant KEGG", max_length=32, blank=True, db_index=True)
    mapping_status = models.CharField("Statut du mapping", max_length=32, default="non trouvé")
    condition = models.CharField("Condition comparée", max_length=80)
    control_condition = models.CharField("Condition contrôle", max_length=80, blank=True)
    value = models.FloatField("Abondance moyenne", null=True, blank=True)
    control_value = models.FloatField("Abondance contrôle", null=True, blank=True)
    log2fc = models.FloatField("log2 Fold Change", null=True, blank=True)
    pvalue = models.FloatField("p-value (Welch)", null=True, blank=True)
    padj_bh = models.FloatField("p-value ajustée (BH)", null=True, blank=True)
    significance = models.CharField("Significativité", max_length=16, default="ns")

    class Meta:
        ordering = ["metabolite_name", "condition"]
        indexes = [
            models.Index(fields=["dataset", "compound_id"]),
        ]

    def __str__(self):
        cid = self.compound_id or "non mappé"
        return f"{self.metabolite_name} ({cid}) [{self.condition}]"