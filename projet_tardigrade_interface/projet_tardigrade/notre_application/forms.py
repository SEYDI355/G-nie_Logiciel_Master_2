from django import forms

from .models import Dataset


class DatasetForm(forms.ModelForm):
    class Meta:
        model = Dataset
        fields = ["name", "description", "matrix_file", "metadata_file"]
        widgets = {"description": forms.Textarea(attrs={"rows": 2})}
