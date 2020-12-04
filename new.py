                       #Project D.A.R.T.
#Coding challenge
#Simillation of RAFT polymerizaion
import numpy as np
from numpy import random
#from pyculib import rand as random
from multiprocessing.dummy import Pool as ThreadPool
pool = ThreadPool()
#from pyculib import
I = float( input ( "I_concentration:"))
M = 5
T = float(input("T_concentrarion:"))
#Frametime = float(input(“frametime in millisecond”))
#For future analysis
frametime = 1
Rx = np.zeros(3)
dTRx = np.zeros(3)
TRx = np.zeros(3)
RmTRx = np.zeros((3,3), dtype = float)
Px = np.zeros(3)
Frametime = frametime /1000
#Reaction Engine In Loop
Loops = 0
while(Loops<3):
	Loops = Loops + 1
	RxA = np.arange(Rx.size)
	dTRxA = np.arange(dTRx.size)
	TRxA = np.arange(TRx.size)
	#Initiation
	R11=0.0000001*2*I
	#propagetion
	R21n=10*Rx*M
	R21 = np.sum(R21n)
	#Pre equilibrium
	#3.1
	R31n=10000*Rx*T
	R31 = np.sum(R31n)
	#3.2/3.3
	R32n=500*dTRx
	R32 = np.sum(R32n)
	#core equilibrium
	#4.1
	R41n = np.zeros((TRx.size,Rx.size))
	R53n = np.zeros((TRx.size,Rx.size,Rx.size))
	R51n = np.zeros((Rx.size,Rx.size))
	for a in RxA:
		for b in TRxA:
			R41n[b,a] = 10000 * Rx[a] * TRx[b]
			for c in RxA:
				R53n[b,a,c] =  100000 * RmTRx[b,a] * Rx[c]
				R51n[a,c] = 100000 * Rx[a] *Rx[c]
	R41 = np.sum (R41n)
	R53 = np.sum (R53n)
	R51 = np.sum (R51n)
	#4.2/4.3
	R42n = 500 * RmTRx
	R42 = np.sum (R42n)
	if Rx[-1] != 0:
		Rx = np.append(Rx,0)
	if np.size(dTRx)<np.size(Rx):
		dTRx = np.append(dTRx, np.zeros(np.size(Rx)-np.size(dTRx)))
	if np.size(TRx)<np.size(dTRx):
		TRx = np.append(TRx, np.zeros(np.size(dTRx)-np.size(TRx)))
	if np.size(RmTRx,axis=1)<np.size(Rx):
		Array = np.zeros(((np.size(RmTRx,axis = 0),1)),dtype = float)
		RmTRx = np.c_[RmTRx,Array]
	if np.size(RmTRx,axis=0) < np.size(TRx):
		Array = np.zeros(((1,np.size(RmTRx,axis = 1))),dtype = float)
		RmTRx = np.r_[RmTRx,Array]
	if np.size(Px) < np.size(RmTRx)+np.size(Rx):
		Px =  np.append(Px,0)
	#Px
	Rand = random.rand()
	Rtotal = R11 + R21 + R31 + R32*2 + R41 + R42*2 + R51*2 + R53*2
	P11  = R11/Rtotal  - Rand
	P21n = R21n/Rtotal - Rand
	P31n = R31n/Rtotal - Rand
	P32n = R32n/Rtotal - Rand
	P41n = R41n/Rtotal - Rand
	P42n = R42n/Rtotal - Rand
	P51n = R51n/Rtotal - Rand
	P53n = R53n/Rtotal - Rand
	#1.1
	if P11 >= 0 :
		 I = I - R11*Frametime
		 Rx[0] = Rx[0] + R11*Frametime*2
	for a in RxA:
		 if  P21n[a] > 0 :
			 Rx[a] -= R21n[a] * Frametime
			 M = M - R21n[a]*Frametime
			 Rx[a+1] = Rx[a+1] +R21n[a]*Frametime
		 if P31n[a] > 0 and Rx[a]>0:
			 Rx[a] = Rx[a] - R31n[a]*Frametime
			 T = T - R31n[a]*Frametime
			 dTRx[a] = dTRx[a] + R31n[a]*Frametime
		 for b in TRxA:
				if P41n[b,a] > 0 :
					Rx[a] = Rx[a] - R41n [b,a] *Frametime
					TRx[b] = TRx[b] - R41n [b,a] *Frametime
					RmTRx[b,a] = RmTRx[b,a] + R41n [b,a] *Frametime
				if P42n[b,a] > 0 :
					Rx[a] = Rx[a] + R41n [b,a] *Frametime
					TRx[b] = TRx[b] + R41n [b,a] *Frametime
					RmTRx[b,a] = RmTRx[b,a] - R41n [b,a] *Frametime
					RmTRx[b,a] = RmTRx[b,a] - R41n [b,a] *Frametime
					Rx[b] = Rx[b] + R41n [b,a] *Frametime
					TRx[a] = TRx[a] + R41n [b,a] *Frametime
				for c in RxA:
					 if P53n[b,a,c] > 0 :
						 RmTRx[b,a] -= R53n[b,a,c]*Frametime
						 Px[a+b+c] += R53n[b,a,c]*Frametime
					 if P51n[a,c] > 0:
						 Rx[a] = Rx[a] - R51n[a,c]*Frametime
						 Rx[c] = Rx[c] - R51n[a,c]*Frametime
						 Rx[a] = Rx[a] - R51n[a,c]*Frametime
						 Rx[c] = Rx[c] - R51n[a,c]*Frametime
						 Px[a+c] = Px[a+c] - R51n[a,c]*Frametime
						 Px[c] = Px[c] - R51n[a,c]*Frametime
						 Px[a] = Px[a] - R51n[a,c]*Frametime
	#3.2/3.3
	for a in dTRxA:
		if P32n[a] > 0 :
			Rx[a] = Rx[a] + R32n[a]*Frametime
			T = T + R32n[a]*Frametime
			dTRx[a] = dTRx[a] - R32n[a]*Frametime
			Rx[0] = Rx[0] + R32n[a]*Frametime
			T = T + R32n[a]*Frametime
			TRx[a] = TRx[a]+ R32n[a]*Frametime
