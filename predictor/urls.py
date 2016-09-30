from django.conf.urls import include, url
from . import views

urlpatterns = [
	url(r'^index/', views.index,name='index'),
	url(r'^contact/', views.contact,name='contact'),
	url(r'^tool/', views.tool,name='tool'),
	url(r'^bibliography/', views.bibliography,name='bibliography'),
	url(r'^results/',views.results,name='results'),
]