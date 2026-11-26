import tempfile
from pathlib import Path

from django.conf import settings
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings
from django.urls import reverse

from .models import Dataset, MetaboliteMeasurement


class NetworkIntegrationTests(TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp_dir.cleanup)

    def _upload_synthetic(self):
        data_dir = Path(settings.BASE_DIR) / "data" / "synthetic"
        matrix_bytes = (data_dir / "synthetic_matrix.csv").read_bytes()
        meta_bytes = (data_dir / "synthetic_metadata.csv").read_bytes()
        return self.client.post(
            reverse("notre_application:home"),
            {
                "name": "Jeu synthétique tardigrade",
                "description": "Hydraté vs Desséché",
                "matrix_file": SimpleUploadedFile("synthetic_matrix.csv", matrix_bytes, content_type="text/csv"),
                "metadata_file": SimpleUploadedFile("synthetic_metadata.csv", meta_bytes, content_type="text/csv"),
            },
        )

    def test_upload_synthetic_and_network_json(self):
        with override_settings(MEDIA_ROOT=self.tmp_dir.name):
            resp = self._upload_synthetic()
            self.assertEqual(resp.status_code, 302)
            ds = Dataset.objects.get()

            trehalose = MetaboliteMeasurement.objects.get(dataset=ds, metabolite_name="Trehalose")
            self.assertEqual(trehalose.compound_id, "C01083")
            self.assertEqual(trehalose.control_condition, "hydraté")
            self.assertEqual(trehalose.condition, "desséché")
            self.assertGreater(trehalose.log2fc, 1.0)
            self.assertEqual(trehalose.significance, "hausse")

            unknown = MetaboliteMeasurement.objects.get(dataset=ds, metabolite_name="Unknown_SUC_01")
            self.assertEqual(unknown.compound_id, "")
            self.assertEqual(unknown.mapping_status, "non trouvé")

            json_resp = self.client.get(reverse("notre_application:network_json", args=[ds.pk]))
            self.assertEqual(json_resp.status_code, 200)
            payload = json_resp.json()

            self.assertEqual(payload["summary"]["selected_condition"], "desséché")
            self.assertEqual(payload["summary"]["control_condition"], "hydraté")
            self.assertEqual(payload["summary"]["mapped_count"], 53)
            self.assertEqual(payload["summary"]["in_graph_count"], 51)
            self.assertIn("C01083", payload["measurements"])
            self.assertTrue(payload["measurements"]["C01083"]["in_kegg_graph"])

            detail_resp = self.client.get(reverse("notre_application:detail", args=[ds.pk]))
            self.assertEqual(detail_resp.status_code, 200)
            self.assertContains(detail_resp, "Réseau métabolique")
            self.assertContains(detail_resp, "cytoscape_network.js")
            self.assertContains(detail_resp, reverse("notre_application:network_json", args=[ds.pk]))

    def test_direct_kegg_id_matrix(self):
        matrix = (
            "compound_id,C1,C2,C3,T1,T2,T3\n"
            "C01083,10,12,11,100,110,105\n"
            "cpd:C00031,100,95,105,10,12,11\n"
        ).encode("utf-8")
        meta = (
            "sample,condition\n"
            "C1,Control\nC2,Control\nC3,Control\n"
            "T1,Treated\nT2,Treated\nT3,Treated\n"
        ).encode("utf-8")

        with override_settings(MEDIA_ROOT=self.tmp_dir.name):
            resp = self.client.post(
                reverse("notre_application:home"),
                {
                    "name": "Direct KEGG",
                    "description": "",
                    "matrix_file": SimpleUploadedFile("m.csv", matrix, content_type="text/csv"),
                    "metadata_file": SimpleUploadedFile("meta.csv", meta, content_type="text/csv"),
                },
            )
            self.assertEqual(resp.status_code, 302)
            ds = Dataset.objects.get()
            payload = self.client.get(reverse("notre_application:network_json", args=[ds.pk])).json()
            self.assertEqual(payload["summary"]["mapped_count"], 2)
            self.assertEqual(payload["measurements"]["C01083"]["mapping_status"], "direct")
            self.assertEqual(payload["measurements"]["C00031"]["significance"], "baisse")

    def test_invalid_metadata_rejects_dataset(self):
        matrix = b"metabolite,S1\nTrehalose,10\n"
        bad_meta = b"bad_col,other\nS1,Control\n"
        with override_settings(MEDIA_ROOT=self.tmp_dir.name):
            resp = self.client.post(
                reverse("notre_application:home"),
                {
                    "name": "Invalide",
                    "description": "",
                    "matrix_file": SimpleUploadedFile("m.csv", matrix, content_type="text/csv"),
                    "metadata_file": SimpleUploadedFile("meta.csv", bad_meta, content_type="text/csv"),
                },
            )
            self.assertEqual(resp.status_code, 200)
            self.assertEqual(Dataset.objects.count(), 0)
            self.assertContains(resp, "Format de fichier invalide")
