from __future__ import unicode_literals

from django.db import models

class Cancer(models.Model):
	nombre = models.CharField(max_length=200)
	descripcion = models.CharField(max_length=2048)

	def __str__(self):
		return self.nombre

class Bibliography(models.Model):
	nombre = models.CharField(max_length=200)
	link = models.CharField(max_length=2048)
	cancer = models.ForeignKey(Cancer, on_delete=models.CASCADE)

	def __str__(self):
		return self.nombre

class Parameters(models.Model):
	tc = models.FloatField()
	delta = models.FloatField()
	u1=models.FloatField()
	u2=models.FloatField()
	C0=models.FloatField()
	eta=models.FloatField()
	psi=models.FloatField()
	gamma0=models.FloatField()
	gamma1=models.FloatField()
	Cs=models.FloatField()
	Mth=models.FloatField()
	xini=models.FloatField()
	alfa=models.FloatField()
	nu2=models.FloatField()
	lmin=models.FloatField()
	epxi=models.FloatField()
	Tmax = models.FloatField()
	T1ini=models.FloatField()
	T1end=models.FloatField()
	T2ini=models.FloatField()
	T2end=models.FloatField()
	W=models.FloatField()
	L=models.FloatField()
	nx=models.IntegerField()
	ny=models.IntegerField()
	rx=models.FloatField()
	ry=models.FloatField()
	e = models.FloatField()
	eth = models.FloatField()
	Sini = models.FloatField()
	qini = models.FloatField()
	flag = models.FloatField()
	cfl = models.FloatField()
	bt = models.FloatField()