###Imports###
import numpy as np
import scipy as sc
from scipy.sparse import linalg
from scipy import sparse
from rogues.matrices import poisson

class Calculator():

	###Functions###

	def bord(self,F):
		### take a function and it remove the boundaries ###
		F = np.delete(F,-1,0) # Ultima fila
		F = np.delete(F,0,0) # Primera fila
		F = np.delete(F,-1,1) # Ultima columna
		F = np.delete(F,0,1) # Primera columna
		return F

	def Fxi(self,P1,P2,P3,xi,g,GammaPP,alfa,gamma0,epxi,l,hx,hy):
		### El termino sum( sum(... representa una aprox. de la integral ###
		F4 = alfa*hx*hy*np.sum(np.sum((1.0+epxi-GammaPP/gamma0)*(g*(P1+P2)+P3)))-l*xi
		return F4

	def antiangio(self,t,Cg,Tgini,Tgend):
		g=Cg*((t>=Tgini) and (t<Tgend))
		g=1.0-g
		return g

	def div(self,GammaPP,P1,P2,P3,delta,N,M):
		return GammaPP*( P1 + P2 + P3 )-delta*N*(M+1.0)

	def grad_vf(self,X,hx,hy):
		m,n=np.shape(X)
		### Extend S by 1 outside of the domain ###
		X=np.vstack((np.ones((1,n),float),X,np.ones((1,n),float)))
		X=np.column_stack((((np.ones((m+2,2)),X,np.ones((m+2,2))))))
		### Now we compute the gradient (Fx,Fy) of -X: ###
		### First we compute in the middle grid points ###
		### (i+1/2,j) for Fx and (i,j+1/2) for Fy ###
		I=np.arange(0,m+1)
		I = I.reshape(121,1)
		Fx = -(X[I+1,np.arange(1,n+1)]-X[I,np.arange(1,n+1)])/hx # (m+1) x n
		J=np.arange(0,n)
		J = J.reshape(120,1)
		Fy = -(X[np.arange(1,m+2),J+1]-X[np.arange(1,m+2),J])/hy # m x (n+1)
		### Second we compute Fx and Fy at center grid points ###
		Fx = (Fx[np.arange(0,m),:] + Fx[np.arange(1,m+1),:])/2.0
		Fy = (Fy[:,np.arange(0,n)] + Fy[:,np.arange(1,n+1)])/2.0
		return Fx,Fy

	def solver_vf(self,A,F,hx,hy):
		m,n=np.shape(F)
		### Computing the pressure ###
		F=np.reshape(F,(m*n,1),order='F')
		F = sc.sparse.csr_matrix(F)
		P=-sc.sparse.linalg.spsolve(A,F)
		i = 0
		P=np.reshape(P,(m,n),order='F')
		### Homogeneous Dirichlet conditions for the pressure ###
		P=np.vstack((np.zeros((1,n),float),P,np.zeros((1,n),float)))
		P=np.column_stack((np.zeros((m+2,1),float),P,np.zeros((m+2,1),float)))
		### Now we compute the velocity field of expansion/decay of the tissue ###
		### First we compute in the middle grid points ###
		### (i+1/2,j) for U and (i,j+1/2) for V ###
		I=np.arange(0,m+1)
		I = I.reshape(121,1)
		U = -(P[I+1,np.arange(1,n+1)]-P[I,np.arange(1,n+1)])/hx # (m+1) x n
		J=np.arange(0,n)
		J = J.reshape(120,1)
		V = -(P[np.arange(1,m+2),J+1]-P[np.arange(1,m+2),J] )/hy # m x (n+1)
		### Second we compute U and V at center grid points ###
		U = (U[np.arange(0,m),:] + U[np.arange(1,m+1),:])/2.0
		V = (V[:,np.arange(0,n)] + V[:,np.arange(1,n+1)])/2.0
		return U,V

	def condInicialMod(self,theta,e,rx,ry,m0,xini,eth,S,q,L,W,X,Y):
		c=(2.0*np.pi)/(2.0*np.pi-np.arccos(1.0-2.0*eth))
		Dx=X-L/2.0
		Dy=Y-W/2.0
		Dxr=Dx*np.cos(theta)+Dy*np.sin(theta)
		Dyr=-Dx*np.sin(theta)+Dy*np.cos(theta)
		D = (np.sqrt((Dxr/(e*Dxr+c*rx))**2+(Dyr/(e*Dyr+c*ry))**2))
		Y=1.0*(D<=0.5)+0.5*(1-np.cos(2*np.pi*D))*(np.logical_and(D>0.5,D<1))
		P1=(1.0-S)*Y
		P2=(S/(1.0+q))*Y
		P3=q*P2
		N=np.zeros(np.shape(X))
		M=m0*np.ones(np.shape(X))
		xi=xini
		S=1.0-P1-P2-P3-N
		return P1,P2,P3,N,S,M,xi

	def evolutiontime_weno5_M_vf(self,vitx,vity,U,U0,lambdax,lambday,flag,m0,A,LD,dt,coef,hx,hy):
		U0xp,U0xm=self.LUx_weno5(U0,flag)
		U0yp,U0ym=self.LUy_weno5(U0,flag)
		I=U0xm.shape[0]
		J=U0ym.shape[1]
		i=np.arange(1,I)
		j=np.arange(1,J)
		Htx = vitx*((U0xp[i,:]-U0xp[i-1,:])*(vitx<0)+(U0xm[i,:]-U0xm[i-1,:])*(vitx>0))
		Hty = vity*((U0yp[:,j]-U0yp[:,j-1])*(vity<0)+(U0ym[:,j]-U0ym[:,j-1])*(vity>0))
		m,n=np.shape(U)
		U = np.vstack((m0*np.ones((1,n)),U,m0*np.ones((1,n))))
		U = np.column_stack((m0*np.ones((m+2,1)),U,m0*np.ones((m+2,1))))
		i = np.arange(1,m+1)
		j = np.arange(1,n+1)
		LapU = ((U[i-1,j]-2*U[i,j]+U[i+1,j])/hx**2)+((U[i,j-1]-2*U[i,j]+U[i,j+1])/hy**2) 
		U=self.bord(U)
		F=-lambdax*Htx-lambday*Hty+0.5*dt*coef*LD+U+0.5*dt*coef*LapU
		F= np.reshape(F,(m*n,1),order='F')
		value = sparse.eye(m*n,m*n)-0.5*dt*coef*A
		Unp1=sc.sparse.linalg.spsolve(value,F)
		Unp1=np.reshape(Unp1,(m,n))
		return Unp1

	def evolutiontime_weno5(self,vitx,vity,U,lambdax,lambday,flag,bd):
		Uxp,Uxm=self.LUx_weno5(U,flag)
		Uyp,Uym=self.LUy_weno5(U,flag)
		I=Uxm.shape[0]
		J=Uym.shape[1]
		i=np.arange(1,I)
		j=np.arange(1,J)
		Vx = vitx*((Uxp[i,:]-Uxp[i-1,:])*(vitx<0)+(Uxm[i,:]-Uxm[i-1,:])*(vitx>0))
		Vy = vity*((Uyp[:,j]-Uyp[:,j-1])*(vity<0)+(Uym[:,j]-Uym[:,j-1])*(vity>0))
		Unp1=U-lambdax*Vx-lambday*Vy
		i= np.absolute(Unp1)<np.spacing(1)
		Unp1[i]=0
		Unp1[1,:]=bd*(vitx[1,:]>0)+Unp1[1,:]*(vitx[1,:]<=0)
		Unp1[-1,:]=bd*(vitx[-1,:]<0)+Unp1[-1,:]*(vitx[-1,:]>=0)
		Unp1[:,1]=bd*(vity[:,1]>0)+Unp1[:,1]*(vity[:,1]<= 0)
		Unp1[:,-1]=bd*(vity[:,-1]<0)+Unp1[:,-1]*(vity[:,-1]>=0)
		return Unp1

	def LUx_weno5(self,U,flag):
		Up=U
		Um=U
		### Natural extrapolation of flag-order (flag>=0; flag=0 means periodic BC) ###
		if flag>0:
			Up=self.extension_1d(Up,flag,0)
			Up=self.extension_1d(Up,flag,0)
			Up=self.extension_weno5(Up,flag,0,1)
			Um=self.extension_1d(Um,flag,0)
			Um=self.extension_1d(Um,flag,0)
			Um=self.extension_weno5(Um,flag,0,0)
		else:
			I=Up.shape[0]
			Up=np.vstack(Up[I-2,:], Up[I-1,:], Up, Up[2,:], Up[3,:], Up[4,:])
			Um=np.vstack(Um[I-3,:], Um[I-2,:], Um[I-1,:], Um, Um[2,:], Um[3,:])
		I=Up.shape[0]
		i=np.arange(2,I-2)
		Up=self.weno5(Up[i+2,:],Up[i+1,:],Up[i,:],Up[i-1,:],Up[i-2,:])
		i=np.arange(3,I-1)
		Um=self.weno5(Um[i-1,:],Um[i,:],Um[i-1,:],Um[i,:],Um[i+1,:])
		return Up,Um

	def LUy_weno5(self,V,flag):
		Vp=V
		Vm=V
		### Natural extrapolation of flag-order (flag varying between 1 and 5) ###
		if flag>0:
			Vp=self.extension_1d(Vp,flag,1)
			Vp=self.extension_1d(Vp,flag,1)
			Vp=self.extension_weno5(Vp,flag,1,1)
			Vm=self.extension_1d(Vm,flag,1)
			Vm=self.extension_1d(Vm,flag,1)
			Vm=self.extension_weno5(Vm,flag,1,0)
		else:
			J=Vp.shape[1]
			Vp=np.concatenate(Vp[:,J-2],Vp[:,J-1],Vp,Vp[:,2],Vp[:,3],Vp[:,4])
			Vm=np.concatenate(Vm[:,J-3],Vm[:,J-2],Vm[:,J-1],Vm,Vm[:,2],Vm[:,3])
		J=Vp.shape[1]
		j=np.arange(2,J-2)
		Vp=self.weno5(Vp[:,j+2],Vp[:,j-1],Vp[:,j],Vp[:,j-1],Vp[:,j-2])
		j=np.arange(3,J-1)
		Vm=self.weno5(Vm[:,j-3],Vm[:,j-2],Vm[:,j-1],Vm[:,j],Vm[:,j+1])
		return Vp,Vm

	def extension_1d(self,X,s,fgd):
		# fgd=0: ext. en la dir. x
		# fgd=1: ext. en la dir. y
		if fgd==1:
			X=np.transpose(X)
		m=X.shape[0]
		# s-th order extrapolation
		if s>=1:
			k=np.arange(1,int(s)+1)
			a = ((-1)**(k+1)*sc.misc.factorial(s)/(sc.misc.factorial(k)*sc.misc.factorial(s-k)))
			b = X[k-1,:]
			a = np.reshape(a,(1,5))
			XexI = a.dot(b)
			a = ((-1)**(k+1)*sc.misc.factorial(s)/(sc.misc.factorial(k)*sc.misc.factorial(s-k)))
			b = X[m-k-1,:]
			XexD = a.dot(b)
		else:
		# Ext. periodica
			XexI=X[m-1,:]
			XexD=X[2,:]
		X=np.vstack((XexI, X, XexD))
		if fgd==1:
			X=np.transpose(X)
		return X

	### s-th order extrapolation to the left(fg_dir=0) or to the right(fg_dir=1) ###
	def extension_weno5(self,X,s,fg_dim,fg_dir):
		if s>=1:
			k=np.arange(1,int(s)+1)
			if fg_dim==0:
				m=X.shape[0]
				if fg_dir==1:
					#XexD=( (-1).^(k+1).*factorial(s)./(factorial(k).*factorial(s-k)) )*X(m+1-k,:);
					a = ((-1)**(k+1)*sc.misc.factorial(s)/(sc.misc.factorial(k)*sc.misc.factorial(s-k)))
					b = X[m-k-1,:]
					XexD = a.dot(b)
					X=np.vstack((X,XexD))
				else:
					a = ((-1)**(k+1)*sc.misc.factorial(s)/(sc.misc.factorial(k)*sc.misc.factorial(s-k)))
					b = X[k-1,:]
					XexI = a.dot(b)
					X=np.vstack((XexI,X))
			else:
				n=X.shape[1]
				if fg_dir==1:
					#XexAr=X(:,n+1-k)*( (-1).^(k+1).*factorial(s)./(factorial(k).*factorial(s-k)) )';
					a = X[:,n-k-1]
					b = ((-1)**(k+1)*sc.misc.factorial(s)/(sc.misc.factorial(k)*sc.misc.factorial(s-k)))
					XexAr = a.dot(b)
					XexAr = np.transpose(XexAr)
					X=np.column_stack((X,XexAr))
				else:
					#XexAb=X(:,k)*( (-1).^(k+1).*factorial(s)./(factorial(k).*factorial(s-k)) )';
					a = X[:,k-1]
					b = ((-1)**(k+1)*sc.misc.factorial(s)/(sc.misc.factorial(k)*sc.misc.factorial(s-k)))
					XexAb = a.dot(b)
					XexAb = np.transpose(XexAb)
					X=np.column_stack((XexAb,X))
		else:
		#Ext. per.
			if fg_dim==0:
				if fg_dir==0:
					m=X.shape[0]
					XexI=X[m-1,:]
					X=np.vstack((XexI, X))
				else:
					XexD=X[2,:]
					X=np.vstakc((X, XexD))
			else:
				if fg_dir==0:
					n=X.shape[1]
					XexAb=X[:,n-1]
					X=np.column_stack((XexAb,X))
				else:
					XexAr=X[:,2]
					X=np.column_stack((X,XexAr))
		return X

	def weno5(self,a,b,c,d,e):
		### Polynomial reconstructions in the three stencils ###
		q1=a/3.0-7.0*b/6.0+11.0*c/6.0
		q2=-b/6.0+5.0*c/6.0+d/3.0
		q3=c/3.0+5.0*d/6.0-e/6.0
		### Smoothness indicators ###
		S1=( 13.0*(a-2.0*b+c)**2+3.0*(a-4.0*b+3.0*c)**2 )/12.0
		S2=( 13.0*(b-2.0*c+d)**2+3.0*(d-b)**2 )/12.0
		S3=( 13.0*(c-2.0*d+e)**2+3.0*(3.0*c-4.0*d+e)**2 )/12.0
		ep=1.0e-6
		a1=1.0/( 10.0*(ep+S1)**2 )
		a2=6.0/( 10.0*(ep+S2)**2 )
		a3=3.0/( 10.0*(ep+S3)**2 )
		z=(a1*q1+a2*q2+a3*q3)/(a1+a2+a3)
		return z

	def solver(self,Tc,delta,u1,u2,C0,eta,psi,gamma0,gamma1,Cs,Mth,m0,xini,alfa,nu2,l,epxi,T1ini,T1end,T2ini,T2end,W,L,nx,ny,hx,hy,h,R,X,Y,A,LD,rx,ry,e,eth,Sini,qini,P10,P20,P30,N0,S0,M0,xi0,t,Tmax,A1,A2,A3,An,At,M1,M2,M3,Mn,Mt,Xi,T,flag,dta,cfl,bt):
		areas = []
		areas.append(At)
		#while(t<Tmax):
		while(t<5.0):
			### *** First step in the Splitting Method *** ###
			GammaPP0=0.5*gamma0*(1.0+np.tanh(R*(M0-Mth)))
			GammaPD0=0.5*gamma1*(1.0-np.tanh(R*(M0-Mth)))
			GammaSD0=Cs*gamma1*np.maximum(0.0,-np.tanh(R*(M0-Mth)))
			### Calcula la divergencia de la velocidad ###
			F0=self.div(GammaPP0,P10,P20,P30,delta,N0,M0)
			### Resuelve la ecuacion de la presion y encuentra la velocidad ###
			U,V=self.solver_vf(A,F0,hx,hy)
			##dt_cfl=min( hx/max( max( abs(U) ) ), hy/max( max( abs(V) ) ) )
			dt_cfl=np.minimum(hx/np.absolute(U).max(),hy/np.absolute(V).max())
			### Velocidad para la ecuacion de M ###
			Um,Vm=self.grad_vf(S0,hx,hy)
			### Vel. normalizadas ###
			vel=np.sqrt(Um**2.0+Vm**2.0)+np.spacing(1)
			Um=(xi0*Um)/vel
			Vm=(xi0*Vm)/vel
			dt_cflxi=h/xi0
			dt_cfl=np.minimum(dt_cfl,dt_cflxi)
			### Celulas proliferativas p1, p2 y p3 ###
			G1 = GammaPP0-GammaPD0-(u1*((t>=T1ini) and (t<T1end)) + u2*((t>=T2ini) and (t<T2end)))*(M0+1.0)
			G2 = GammaPP0-GammaPD0-(u2*((t>=T2ini) and (t<T2end)))*(M0+1.0)
			G3 = GammaPP0-GammaPD0
			Gn = delta*(M0+1.0)
			dtW=np.minimum(1.0/np.absolute(G1-F0).max(),1.0/np.absolute(G2-F0).max())
			dtW=np.minimum(dtW,1.0/np.absolute(G3-F0).max())
			dtW=np.minimum(dtW,1.0/np.absolute(Gn-F0).max())
			dtW=np.minimum(dtW,1.0/np.absolute(GammaSD0+F0).max())
			dt_cfl=np.minimum(dt_cfl,dtW)
			dt=cfl*np.minimum(dta,dt_cfl)
			P1=P10+(dt/2.0)*(G1-F0)*P10
			P2=P20+(dt/2.0)*(G2-F0)*P20
			P3=P30+(dt/2.0)*(G3-F0)*P30
			### Celulas necroticas ### 
			N = N0+(dt/2.0)*(-(F0+Gn)*N0+GammaPD0*(P10+P20+P30)+GammaSD0*(1.0-P10-P20-P30-N0)+(u1*((t>=T1ini) and (t<T1end))*P10+u2*((t>=T2ini) and (t<T2end))*(P10+P20))*(M0+1.0))
			### Nutrientes y flujo sang. ###
			M = M0+(dt/2.0)*(C0*S0*(1.0-M0/(2.0*Mth))-eta*(P10+P20+P30)*M0)
			### Resolucion de la ec. de xi ###
			g=self.antiangio(t,nu2,T2ini,T2end)
			Fchi=self.Fxi(P10,P20,P30,xi0,g,GammaPP0,alfa,gamma0,epxi,l,hx,hy)
			xi=xi0+dt*Fchi
			### Cheking the max and min of variables ###
			I=np.logical_or((np.logical_and(P1<0,P1>-1e-5)),(np.absolute(P1)<np.spacing(1)))
			P1[I]=0
			I=np.logical_or((np.logical_and(P2<0,P2>-1e-5)),(np.absolute(P2)<np.spacing(1)))
			P2[I]=0
			I=np.logical_or((np.logical_and(P3<0,P3>-1e-5)),(np.absolute(P3)<np.spacing(1)))
			P3[I]=0
			I=np.logical_or((np.logical_and(N<0,N>-1e-5)),(np.absolute(N)<np.spacing(1)))
			N[I]=0
			i = 0
			#************************************************************************#
			### *** Second step in the Splitting Method *** ###
			lambdax=dt/hx
			lambday=dt/hy
			P1=self.evolutiontime_weno5(U,V,P1,lambdax,lambday,flag,0.0)
			P2=self.evolutiontime_weno5(U,V,P2,lambdax,lambday,flag,0.0)
			P3=self.evolutiontime_weno5(U,V,P3,lambdax,lambday,flag,0.0)
			N=self.evolutiontime_weno5(U,V,N,lambdax,lambday,flag,0.0)
			M=self.evolutiontime_weno5_M_vf(Um,Vm,M,M0,lambdax,lambday,flag,m0,A,LD,dt,psi,hx,hy)
			### Cheking the max and min of variables ###
			I=np.logical_or((np.logical_and(P1<0,P1>-1e-5)),(np.absolute(P1)<np.spacing(1)))
			P1[I]=0
			I=np.logical_or((np.logical_and(P2<0,P2>-1e-5)),(np.absolute(P2)<np.spacing(1)))
			P2[I]=0
			I=np.logical_or((np.logical_and(P3<0,P3>-1e-5)),(np.absolute(P3)<np.spacing(1)))
			P3[I]=0
			I=np.logical_or((np.logical_and(N<0,N>-1e-5)),(np.absolute(N)<np.spacing(1)))
			N[I]=0
			### ************************************************************************ ###
			### *** Third step in the Splitting Method *** ###
			P1 = P1+(dt/2.0)*( G1-F0 )*P10
			P2 = P2+(dt/2.0)*( G2-F0 )*P20
			P3 = P3+(dt/2.0)*( G3-F0 )*P30
			N = N+(dt/2.0)*(-(F0+Gn)*N+GammaPD0*(P10+P20+P30)+GammaSD0*(1.0-P10-P20-P30-N0)+(u1*((t>=T1ini) and (t<T1end))*P10+u2*((t>=T2ini) and (t<T2end))*(P10+P20))*(M0+1.0))
			### Defining S
			S=1.0-(P1+P2+P3+N)
			I=S>1.0
			S[I]=1.0
			P1[I]=0.0
			P2[I]=0.0
			P3[I]=0.0
			N[I]=0.0
			### Nutrientes y flujo sang. (MODIFICADO) ###
			M = M + (dt/2.0)*(C0*S0*(1.0-M0/(2.0*Mth))-eta*( P10+P20+P30 )*M0)
			### Cheking the max and min of variables ###
			I=np.logical_or((np.logical_and(P1<0,P1>-1e-5)),(np.absolute(P1)<np.spacing(1)))
			P1[I]=0
			I=np.logical_or((np.logical_and(P2<0,P2>-1e-5)),(np.absolute(P2)<np.spacing(1)))
			P2[I]=0
			I=np.logical_or((np.logical_and(P3<0,P3>-1e-5)),(np.absolute(P3)<np.spacing(1)))
			P3[I]=0
			I=np.logical_or((np.logical_and(N<0,N>-1e-5)),(np.absolute(N)<np.spacing(1)))
			N[I]=0
			I=np.logical_or((np.logical_and(S<0,S>-1e-5)),(np.absolute(S)<np.spacing(1)))
			S[I]=0.0
			### Updating of variables for the next iteration ###
			t=t+dt
			P10=P1
			P20=P2
			P30=P3
			S0=S
			N0=N
			xi0=xi
			M0=M	
			A1=100.0*hx*hy*np.sum( np.sum( 1.0*(P10>eth) ) )
			A2=100.0*hx*hy*np.sum( np.sum( 1.0*(P20>eth) ) )
			A3=100.0*hx*hy*np.sum( np.sum( 1.0*(P30>eth) ) )
			An=100.0*hx*hy*np.sum( 1.0*(N0>eth) ) 
			At=100.0*hx*hy*np.sum( np.sum( 1.0*( P10+P20+P30+N0>eth ) ) )
			M1=100.0*hx*hy*np.sum( np.sum( P10 ) )
			M2=100.0*hx*hy*np.sum( np.sum( P20 ) )
			M3=100.0*hx*hy*np.sum( np.sum( P30 ) )
			Mn=100.0*hx*hy*np.sum( N0 )
			Mt=100.0*hx*hy*np.sum( np.sum( (P10+P20+P30+N0) ) )
			Xi=xini
			T=Tc*t
			areas.append(At)
		return areas

	###Main###
	def initAndSolver(self,Tc,delta,u1,u2,C0,eta,psi,gamma0,gamma1,Cs,Mth,xini,alfa,nu2,l,epxi,Tmax,T1ini,T1end,T2ini,T2end,W,L,nx,ny,rx,ry,e,eth,Sini,qini,flag,cfl,bt):
		###Data initialization###
		### Characteristic time (one month) ###
		#Tc=30.0
		### Elimination rate of N ###
		#delta = 0.4
		### Cell death rate due to the drugs ###
		#u1=0.21
		#u2=0.142
		### Angiogenic capacity of healthy tissue ###
		#C0=1.0
		### Rate of comsuption for M due to the tumor ###
		#eta=2.0
		### Diffusion rate of the oxygen ###
		#psi=0.4
		### Rate of proliferation ###
		#gamma0=0.638
		### Dead parameter ###
		#gamma1=0.24
		### Healthy tissue apoptosis relative to gamma1 ###
		#Cs=10.0
		### Parameters of blood vessels ###
		#Mth=2.0
		### Parameter for controlling increasing of M ###
		m0=2.0*Mth
		### Parameters for xi ###
		### Initial condition ###
		#xini=0.1
		### Angiogenic excitability ###
		#alfa=1.0
		### Decay due to antiangiogenic drug ###
		#nu2=0.8
		### Decay (mean life) ###
		#l=0.6
		### Residual production of growth factor ###
		#epxi=0.1
		### * Time of delivery * ###
		### For antiproliferative drug ###
		T1ini=T1ini/Tc
		T1end=T1end/Tc
		### For antiangiogenic drug ###
		T2ini=T2ini/Tc
		T2end=T2end/Tc
		### Dimensions of the computational domain ###
		#W=6.0
		#L=6.0
		### nx: numbers of interior grid nodes in x-direction ###
		### ny: idem in the y-direction ### 
		#nx=119
		#ny=119
		hx=L/(nx+1.0)
		hy=W/(ny+1.0)
		h=np.minimum(hx,hy)
		R=5.0
		X,Y = np.meshgrid(np.arange(hx/2,L-hx/2+hx,hx), np.arange(hy/2,W-hy/2+hy,hy))
		X = np.transpose(X)
		Y = np.transpose(Y)
		### Matrix for solving diffusion equations (for pressure and M) ###
		A=-poisson(nx+1)/h**2.0
		### Boundary conditions for diffusion part in eq. of M  ###
		LD=np.zeros((nx+1,ny+1))
		LD[np.arange(1,nx),0]=m0/(hy**2.0)
		LD[np.arange(1,nx),ny]=m0/(hy**2.0)
		LD[0,np.arange(1,ny)]=m0/(hx**2.0)
		LD[nx,np.arange(1,ny)]=m0/(hx**2.0)
		### Corner grid points ###
		LD[0,0]=m0*((1.0/hx**2)+(1.0/hy**2))
		LD[nx,0]=m0*((1.0/hx**2)+(1.0/hy**2))
		LD[0,ny]=m0*((1.0/hx**2)+(1.0/hy**2))
		LD[nx,ny]=m0*((1.0/hx**2)+(1.0/hy**2))
		### * Initial conditions * ###
		### Radii ###
		#rx=0.47 
		#ry=0.36
		### Excentricity ###
		#e=0.35
		### Minimal threshold for the numerical location of the tumor ###
		#eth=1.0e-2
		#Sini=3.0e-6
		#qini=7.5e-3
		### Initial condition ####
		P10,P20,P30,N0,S0,M0,xi0=self.condInicialMod((25.0*np.pi)/32.0,e,rx,ry,m0,xini,eth,Sini,qini,L,W,X,Y)
		#### * Initialization of time, iteration, and variables to stock tumor areas * ###
		t=0.0
		Tmax = Tmax/Tc
		### Computing the areas of tumor cells (in mm^2) ###
		A1=100.0*hx*hy*np.sum( np.sum( 1.0*(P10>eth) ) )
		A2=100.0*hx*hy*np.sum( np.sum( 1.0*(P20>eth) ) )
		A3=100.0*hx*hy*np.sum( np.sum( 1.0*(P30>eth) ) )
		An=100.0*hx*hy*np.sum( 1.0*(N0>eth) ) 
		At=100.0*hx*hy*np.sum( np.sum( 1.0*( P10+P20+P30+N0>eth ) ) )
		M1=100.0*hx*hy*np.sum( np.sum( P10 ) )
		M2=100.0*hx*hy*np.sum( np.sum( P20 ) )
		M3=100.0*hx*hy*np.sum( np.sum( P30 ) )
		Mn=100.0*hx*hy*np.sum( N0 )
		Mt=100.0*hx*hy*np.sum( np.sum( (P10+P20+P30+N0) ) )
		Xi=xini
		T=Tc*t
		### We stock time in order to plot the areas ###
		### Order for extrapolation ###
		#flag=5.0
		### Computing initial time step size ###
		### For angiogenesis signal equation ###
		dta=np.minimum(1.0/eta,1.0/l)
		### CFL number ###
		#cfl=0.4
		dta=np.minimum(Tc*h,dta);
		### Parameter for Twin-Weno5 method ###
		#bt=0.3
		areas = self.solver(Tc,delta,u1,u2,C0,eta,psi,gamma0,gamma1,Cs,Mth,m0,xini,alfa,nu2,l,epxi,T1ini,T1end,T2ini,T2end,W,L,nx,ny,hx,hy,h,R,X,Y,A,LD,rx,ry,e,eth,Sini,qini,P10,P20,P30,N0,S0,M0,xi0,t,Tmax,A1,A2,A3,An,At,M1,M2,M3,Mn,Mt,Xi,T,flag,dta,cfl,bt)
		return areas
