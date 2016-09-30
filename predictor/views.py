from django.shortcuts import render
from django.http import HttpResponse
from django.template import loader

from .models import Cancer,Bibliography,Parameters
from . import analysis_func

def index(request):
	#template = loader.get_template('predictor/index.html')
	#return HttpResponse(template.render(request))
	return render(request,'predictor/index.html')

def contact(request):
	#template = loader.get_template('predictor/contact.html')
	#return HttpResponse(template.render(request))
	return render(request,'predictor/contact.html')

def tool(request):
	#template = loader.get_template('predictor/tool.html')
	#return HttpResponse(template.render(request))
	parameters_list = Parameters.objects.get(id=1)

	return render(request,'predictor/tool.html',{'parameters_list':parameters_list})

def results(request,chartID = 'chart_ID', chart_type = 'line', chart_height = 500):
	Tc = float(str(request.POST.get('Tc')))
	delta = float(str(request.POST.get('delta')))
	u1 = float(str(request.POST.get('u1')))
	u2 = float(str(request.POST.get('u2')))
	C0 = float(str(request.POST.get('C0')))
	eta = float(str(request.POST.get('eta')))
	psi = float(str(request.POST.get('psi')))
	gamma0 = float(str(request.POST.get('gamma0')))
	gamma1 = float(str(request.POST.get('gamma1')))
	Cs = float(str(request.POST.get('Cs')))
	Mth = float(str(request.POST.get('Mth')))
	xini = float(str(request.POST.get('xini')))
	alfa = float(str(request.POST.get('alfa')))
	nu2 = float(str(request.POST.get('nu2')))
	l = float(str(request.POST.get('l')))
	epxi = float(str(request.POST.get('epxi')))
	Tmax = float(str(request.POST.get('Tmax')))
	T1ini = float(str(request.POST.get('T1ini')))
	T1end = float(str(request.POST.get('T1end')))
	T2ini = float(str(request.POST.get('T2ini')))
	T2end = float(str(request.POST.get('T2end')))
	W = float(str(request.POST.get('W')))
	L = float(str(request.POST.get('L')))
	nx = int(str(request.POST.get('nx')))
	ny = int(str(request.POST.get('ny')))
	rx = float(str(request.POST.get('rx')))
	ry = float(str(request.POST.get('ry')))
	e = float(str(request.POST.get('e')))
	eth = float(str(request.POST.get('eth')))
	Sini = float(str(request.POST.get('Sini')))
	qini = float(str(request.POST.get('qini')))
	flag = float(str(request.POST.get('flag')))
	cfl = float(str(request.POST.get('cfl')))
	bt = float(str(request.POST.get('bt')))

	dater =  analysis_func.ChartData()
	data = dater.check_valve_data(Tc,delta,u1,u2,C0,eta,psi,gamma0,gamma1,Cs,Mth,xini,alfa,nu2,l,epxi,Tmax,T1ini,T1end,T2ini,T2end,W,L,nx,ny,rx,ry,e,eth,Sini,qini,flag,cfl,bt)

	chart = {"renderTo": chartID, "type": chart_type, "height": chart_height,}  
	title = {"text": 'GIST Tumor Growth'}
	xAxis = {"title": {"text": 'Days'}, "categories": data['days']}
	yAxis = {"title": {"text": 'Area'}}
	series = [
		{"name": 'Patient', "data": data['area']}
	]

	parameters_list = Parameters.objects.get(id=1)
	parameters_list.tc = Tc
	parameters_list.delta = delta
	parameters_list.u1 = u1
	parameters_list.u2 = u2
	parameters_list.C0 = C0
	parameters_list.eta = eta
	parameters_list.psi = psi
	parameters_list.gamma0 = gamma0
	parameters_list.gamma1 = gamma1
	parameters_list.Cs = Cs
	parameters_list.Mth = Mth
	parameters_list.xini = xini
	parameters_list.alfa = alfa
	parameters_list.nu2 = nu2
	parameters_list.lmin = l
	parameters_list.epxi = epxi
	parameters_list.Tmax = Tmax
	parameters_list.T1ini = T1ini
	parameters_list.T1end = T1end
	parameters_list.T2ini = T2ini
	parameters_list.T2end = T2end
	parameters_list.W = W
	parameters_list.L = L
	parameters_list.nx = nx
	parameters_list.ny = ny
	parameters_list.rx = rx
	parameters_list.ry = ry
	parameters_list.e = e
	parameters_list.eth = eth
	parameters_list.Sini = Sini
	parameters_list.qini = qini
	parameters_list.flag = flag
	parameters_list.cfl = cfl
	parameters_list.bt = bt
	parameters_list.save()

	return render(request, 'predictor/results.html', {'chartID': chartID, 'chart': chart,
	                                            'series': series, 'title': title, 
	                                            'xAxis': xAxis, 'yAxis': yAxis})

def bibliography(request):
	template = loader.get_template('predictor/biblio.html')
	cancer_list = Cancer.objects.order_by('nombre')
	biblio_list = Bibliography.objects.order_by('cancer')
	numbers = []
	i = 1
	for element in cancer_list:
		numbers.append(i)
		numbers.append(i)
		i = i + 1
	context = {
		'cancer_list' : cancer_list,
		'numbers' : numbers,
		'biblio_list' : biblio_list,
	}
	return render(request,'predictor/biblio.html',{'cancer_list' : cancer_list,
													'numbers' : numbers,
													'biblio_list' : biblio_list,})