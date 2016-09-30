from . import predictor_func

class ChartData():    
	def check_valve_data(self,Tc,delta,u1,u2,C0,eta,psi,gamma0,gamma1,Cs,Mth,xini,alfa,nu2,l,epxi,Tmax,T1ini,T1end,T2ini,T2end,W,L,nx,ny,rx,ry,e,eth,Sini,qini,flag,cfl,bt):
		data = {'days': [], 'area': []}

		predictor = predictor_func.Calculator()

		areas = predictor.initAndSolver(Tc,delta,u1,u2,C0,eta,psi,gamma0,gamma1,Cs,Mth,xini,alfa,nu2,l,epxi,Tmax,T1ini,T1end,T2ini,T2end,W,L,nx,ny,rx,ry,e,eth,Sini,qini,flag,cfl,bt)

		i = 1

		for area in areas:
			data['area'].append(area)
			data['days'].append(i)
			i = i + 1

		return data