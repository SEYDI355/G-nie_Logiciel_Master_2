from django.contrib import admin

from .models import Dataset, MetaboliteMeasurement

admin.site.register(Dataset)
admin.site.register(MetaboliteMeasurement)
