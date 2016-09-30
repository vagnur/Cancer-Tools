from django.contrib import admin

from .models import Cancer
from .models import Bibliography
from .models import Parameters

admin.site.register(Cancer)
admin.site.register(Bibliography)
admin.site.register(Parameters)