# -*- coding: utf-8 -*-
# Generated by Django 1.11.dev20160928140452 on 2016-09-28 19:07
from __future__ import unicode_literals

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Bibliography',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('nombre', models.CharField(max_length=200)),
                ('link', models.CharField(max_length=2048)),
            ],
        ),
        migrations.CreateModel(
            name='Cancer',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('nombre', models.CharField(max_length=200)),
                ('descripcion', models.CharField(max_length=2048)),
            ],
        ),
        migrations.CreateModel(
            name='Parameters',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tc', models.FloatField()),
                ('delta', models.FloatField()),
                ('u1', models.FloatField()),
                ('u2', models.FloatField()),
                ('C0', models.FloatField()),
                ('eta', models.FloatField()),
                ('psi', models.FloatField()),
                ('gamma0', models.FloatField()),
                ('gamma1', models.FloatField()),
                ('Cs', models.FloatField()),
                ('Mth', models.FloatField()),
                ('xini', models.FloatField()),
                ('alfa', models.FloatField()),
                ('nu2', models.FloatField()),
                ('lmin', models.FloatField()),
                ('epxi', models.FloatField()),
                ('Tmax', models.FloatField()),
                ('T1ini', models.FloatField()),
                ('T1end', models.FloatField()),
                ('T2ini', models.FloatField()),
                ('T2end', models.FloatField()),
                ('W', models.FloatField()),
                ('L', models.FloatField()),
                ('nx', models.IntegerField()),
                ('ny', models.IntegerField()),
                ('rx', models.FloatField()),
                ('ry', models.FloatField()),
                ('e', models.FloatField()),
                ('eth', models.FloatField()),
                ('Sini', models.FloatField()),
                ('qini', models.FloatField()),
                ('flag', models.FloatField()),
                ('cfl', models.FloatField()),
                ('bt', models.FloatField()),
            ],
        ),
        migrations.AddField(
            model_name='bibliography',
            name='cancer',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='predictor.Cancer'),
        ),
    ]
