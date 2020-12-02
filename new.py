                       #Project D.A.R.T.
#Coding challenge 
#Simillation of RAFT polymerizaion
import numpy as np
from numpy import random
#from pyculib import rand as random
import _thread
#from pyculib import 
I = float( input ( "I_concentration:"))
M = 5
T = float(input("T_concentrarion:"))
Frametime = 1
Rx = np.zeros(1)
dTRx = np.zeros(1)
TRx = np.zeros(1)
RmTRx = np.zeros((1,2), dtype = float)
Px = np.zeros(1)
#Frametime = float(input(“frametime in millisecond”))
#For future analysis
#Reaction Engine In Loop
Loops = 0
while(Loops<5):
	Loops = Loops+1
	#Initiation
	R11=0.0000001*2*I
	#propagetion
	RxA = np.arange(Rx.size)
	dTRxA = np.arange(dTRx.size)
	R21n = np.zeros( Rx.size )
	for N in RxA:
		R21n[N]=10*Rx[N]*M
	R21 = np.sum(R21n)
	#Pre equilibrium
	#3.1
	R31n = np.zeros(Rx.size)
	for a in RxA :
		R31n[a]=10000*Rx[a]*T
	R31 = np.sum(R31n)
	#3.2
	R32n = np.zeros(dTRx.size)
	for a in dTRxA :
		R32n[a]=500*dTRx[a]
	R32 = np.sum(R32n)
	#core equilibrium
	#Termination
	#Px
	Rand= random.rand()
	Rtotal = R11 + R21 + R31 + R32
	P11=R11/Rtotal- Rand
	P21n = np.zeros(R21n.size)
	P31n = np.zeros(R31n.size)
	P32n = np.zeros(R32n.size)
	for a in RxA :
		P21n[a] = R21n[a]/Rtotal
		P31n[a] = R31n[a]/Rtotal
	for a in dTRxA :
		P32n[a] = R32n[a]/Rtotal
	#Array update
	if Rx[-1] != 0:
	 Rx = np.append(Rx,0)
	if dTRx[-1] != 0:
	 dTRx = np.append(dTRx,0)
	#1.1
	if P11 >= 0 :
	 I = I - R11*Frametime
	 Rx[0] = Rx[0] + R11*Frametime*2
	#2.1
	for k in RxA:
	 if  P21n[k]- Rand > 0 :
	  Rx[k] = Rx[k] - R21[k]*Frametime
	  M = M - R21[k]*Frametime
	  Rx[k+1] = Rx[k+1] +R21[k]*Frametime
	#3.1
	for a in RxA:
	 if P31n[a]-Rand > 0 :
	  Rx[a] = Rx[a] - R31n[a]*Frametime
	  T = T - R31n[a]*Frametime
	  dTRx[a] = dTRx[a] + R31n[a]*Frametime
			
	#3.2/3.3
	for a in dTRxA:
	 if P32n[a]-Rand > 0 :
	  Rx[a] = Rx[a] + R32n[a]*Frametime
	  T = T + R32n[a]*Frametime
	  dTRx[a] = dTRx[a] - R32n[a]*Frametime
	  Rx[0] = Rx[0] + R32n[a]*Frametime
	  T = T + R32n[a]*Frametime
	  TRx[a] = TRx[a]+ R32n[a]*Frametime
print (I)
			
	
			
			


	
	
	   
    
     








