import numpy as np
from util import *
import sys
import copy



## DMODriv
def DMODriv(A,B,Freq,NAtm3,NVib):
    Tresh = 1e-50
    # CH3+H2 -> 1e-10 -> 1e-17
    # HCOOH  -> 1e-10
    # PH2H2  -> 1e-7
    #

#    print A
#    print "Before RotDg1"

    A,B,TMat = RotDg1(0,NAtm3,NVib,A,B)
#    print "A", A
#    print "TMat",TMat

#    print "After RotDg1"
    Freq,NSp,NDim,Lst1 = AssSp1(NVib,Freq,Tresh)     #BUG !
#    print "Freq",Freq
#    print "NSp",NSp
#    print "NDim",NDim
#    print "Lst1",Lst1

    NSp,TMat,NDim,Lst1,Lst2 = AssSp2(NVib,NSp,TMat,NDim,Lst1)
#    print "NSp",NSp
#    print "TMat",TMat,
#    print "NDim",NDim
#    print "Lst1",Lst1
#    print "Lst2",Lst2

#    print "Before RotDeg"

    B,NDim = RotDeg(NVib,NAtm3,NSp,A,B,NDim,Lst1,Lst2,TMat)
#    print "B",B
#    print "NDim",NDim
    B,Freq = Order1(NAtm3,NVib,NSp,NDim,Lst1,Lst2,B,Freq)

    return (B,Freq)


## Order1
def Order1(NParm,NVib,NSp,NDim,Lst1,Lst2,B,Freq):
    Scr=[]
    for i in range(NVib):
	Scr.append([])
    for i in range(NVib):
	for j in range(NParm):
            Scr[i].append(0)
	
    for Mu in range(NVib):
	Scr[Mu][0] = Freq[Mu]   # Corrected
    for ISp in range(NSp):
        NDim1 = NDim[ISp]   
        for K in range(NDim1): # Corrected
            I1 = Lst1[K][ISp]
	    I2 = Lst2[K][ISp]
            Freq[I2-1] = Scr[I1-1][0] # Corrected
    Scr = copy.deepcopy(B)
    for ISp in range(NSp):
	NDim1 = NDim[ISp]
        for K in range(NDim1): # Corrected
            I1 = Lst1[K][ISp]
	    I2 = Lst2[K][ISp]
            for i in range(NParm):
		B[I2-1][i] = Scr[I1-1][i]
		#Scr[i][I1] = B[i][I2]
    for ISp in range(NSp):
	NDim1 = NDim[ISp]
        for K in range(NDim1): # Corrected
	    I1 = Lst1[K][ISp]
	    I2 = Lst2[K][ISp]

    return (B,Freq)

########################################

## RotDeg
def RotDeg(NVib,NParm,NSp,A,B,NDim,Lst1,Lst2,T):
   # A1= matrix( [1,2,3, ...,NParm],
   #             [2,...
   #             [3,...
   #             [4,...
   #             [.
   #             [NVib

   # B1= matrix( [1,2,3, ...,NParm],
   #             [2,...
   #             [3,...
   #             [4,...
   #             [.
   #             [NVib

   # NParm = 3*NAtom

    A1=[]
    B1=[]
    for i in range(NVib):
	A1.append([])
	B1.append([])
    for i in range(NVib):
	for j in range(NParm):
            A1[i].append(0)
	    B1[i].append(0)


    for ISp in range(NSp):
	NDim1 = NDim[ISp]
	for K in range(NDim1): # Corrected
            for i in range(NParm):
                A1[K][i] = A[Lst2[K][ISp]-1][i]
                B1[K][i] = B[Lst1[K][ISp]-1][i]  # Corrected
                #A1[i][K] = A[i][Lst2[K][ISp]]
                #B1[i][K] = B[i][Lst2[K][ISp]]
        #print "NDim1",NDim1
	#print B1
	A1 = A1[0:(NDim1)]
	B1 = B1[0:(NDim1)]
        A1,B1,T = RotDg1(1,NParm,NDim1,A1,B1)
             #def RotDg1(RotB,NParm,NVib,A,B)
	B1=np.array(B1)
	for K in range(NDim1): # Corrected 
            for i in range(NParm):
                B[Lst1[K][ISp]-1][i] = B1[K][i]
                #B[i][Lst1[K][ISp]] = B1[i][K]

    return (B,NDim)
###############################################

## AssSp2
def AssSp2(NVib,NSp,T,NDim,Lst1):
    #print "NVib=",NVib,"NSp=",NSp
    T = np.array(T) 
    Ampl=[]
    for i in range(NSp):
	Ampl.append([])
    for i in range(NSp):
	for j in range(NVib):
	    Ampl[i].append(0) # Corrected
            #Ampl[i][j] = 0
    Scr=[]
    for i in range(NVib):
	Scr.append([])
    for i in range(NVib):
	for j in range(NVib):
	    Scr[i].append(0)

    Lst2=[]
    for i in range(NVib):
        Lst2.append([])
    for i in range(NVib):
	for j in range(NVib):
            Lst2[i].append(0)
 

    for ISp in range(NSp):
	NDim1 = NDim[ISp] #- 1 # correction for python
        for Mu in range(NVib):
            X = 0.0
	    for K in range(NDim1): # Corrected 
		Y = T[Mu][Lst1[K][ISp]-1]
		#print type(T)
                X = X + Y**2
            Ampl[ISp][Mu] = X

    for ISp in range(NSp):
        X = 0.0
        for Mu in range(NVib):
            X = X + Ampl[ISp][Mu]
        for Mu in range(NVib):
            Ampl[ISp][Mu] = Ampl[ISp][Mu] / X

    #Scr = copy.deepcopy(Ampl)
    for i in range(NSp):
	for j in range(NVib):
            Scr[i][j] = Ampl[i][j]

    for ISp in range(NSp):
        NDim1 = NDim[ISp]
        SumAmp = 0.0
        for K in range(NDim1): # Corrected
	    XMax = max( Scr[ISp] )  #Scr is a list 
            IMax = Scr[ISp].index(XMax)  # Corrected 
            Lst2[K][ISp] = IMax + 1  # Corrected
	    for Mu in range(NVib):
		Scr[Mu][IMax] = 0.0
		#Scr[IMax][Mu] = 0.0
            SumAmp = SumAmp + XMax
    
    IErr = ChkLs2(NVib,NSp,NDim,Lst2)
    if IErr != 0:
       print ("AssSp2: Inconsistency in space assignment! \nTry smaller NSTEP.")
       sys.exit()

    return (NSp,T,NDim,Lst1,Lst2)

##################################################

## ChkLs2
def ChkLs2(NVib,NSp,NDim,Lst2):
    IErr = 0
    for ISp1 in range(NSp):
	NDim1 = NDim[ISp1]
	for K1 in range(NDim1): # Corrected
            L1 = Lst2[K1][ISp1]
            for ISp2 in range(NSp):
                NDim2 = NDim[ISp2]
		for K2 in range(NDim2): # Corrected
                    L2 = Lst2[K2][ISp2]
                    Same =( (ISp1==ISp2) and (K1==K2) )
                    if (not Same):
		       if L1 == L2:
			  IErr =1

    return IErr 
#####################################################

## AssSp1
def AssSp1(NVib,Freq,Tresh):
    #Tresh = 1.0e-3 # For mass reation, this value is set to small number

    NDim=[]
    Ind=[]
    Freq1=[]
    Freq2=[]
    for i in range(NVib):
	NDim.append(0)
	Ind.append(0)
	Freq1.append(0)
	Freq2.append(0)


    Lst1=[]
    for i in range(NVib):
        Lst1.append([])
    for i in range(NVib):
	for j in range(NVib):
            Lst1[i].append(0)


    ISp = 0   # correction for python
    LOld = NVib
    for Mu in range(NVib): # start from 0
	Ind[Mu] = Mu + 1
	Freq1[Mu] = Freq[Mu] 
    goto5 = 1
    #print "Ind",Ind
    #print "Freq1",Freq1
    while goto5 == 1:
	  ISp = ISp + 1
	  FrqMin = min(Freq1[:LOld])
          #print "ISp,FrqMin",ISp,FrqMin 
          Freq2[ISp-1] = FrqMin     ###########
          K = 0  # correction for python
	  LNew = 0  # correction for python
          for Mu in range(LOld):
              #print "RUN"
              Delta = abs(FrqMin - Freq1[Mu])
              if Delta <= Tresh :             ######### Tresh
		 K = K + 1
		 #print "RUN"
		 #print "Ind",Ind
		 #print "K,ISp",K,ISp
		 Lst1[K-1][ISp-1] = Ind[Mu] ############
              else:
		 LNew = LNew + 1
                 Ind[LNew-1] = Ind[Mu]
		 Freq1[LNew-1] = Freq1[Mu]
          NDim[ISp-1] = K   ## ??
          LOld = LNew 
	  if LNew <= 0:   ## ??
	     goto5 = 0
	  else:
             goto5 = 1 # repeat the while loop
    NSp = ISp

    return (Freq,NSp,NDim,Lst1) 
#################################################

## RotDg1 
def RotDg1(RotB,NParm,NVib,A,B):
    "Interfacing RotDg1 into python"
    
   # A    [NVib][NParm]  
   # B    [NVib][NParm]
   # Sab  [Nvib][Nvib ]

   # A = matrix( [1,2,3, ...,NParm],
   #             [2,...
   #             [3,...
   #             [4,...
   #             [.
   #             [NVib

   # B = matrix( [1,2,3, ...,NParm],
   #             [2,...
   #             [3,...
   #             [4,...
   #             [.
   #             [NVib

   # Sab=matrix( [1,2,3, ...,NVib2],
   #             [2,...
   #             [3,...
   #             [4,...
   #             [.
   #             [NVib1

   # Scr=matrix( [1,2,3, ...,NVib1],
   #             [2,...
   #             [3,...
   #             [4,...
   #             [.
   #             [NVib1


   # T = matrix( [1,2,3, ...,NV],
   #             [2,...
   #             [3,...
   #             [4,...
   #             [.
   #             [NV


    #print A
    #print "A stop"

    Sab = MatMp1(A,B)
    
    #print Sab

    Scr = MatMp1(Sab,Sab) # symmetric matrix
  
    #print Scr

    Scr1= SqrtMp(-1,Scr)

    T   = MPACMF(Scr1,Sab,3) # ? 

    Scr = MMpyMF(Sab, T)

    if RotB == 1:
       Scr = MMpyMF(B,T)
       B = copy.deepcopy(Scr)  

    return (A,B,T)
	

###############################################

def ABscIn(m,n,A,B):
    # C = (B-A)* (m/n) + A
    m = float(m)
    n = float(n)
    if m > n:
       print ("Error: Illegal m value. m should not exceed n in ABscIn!")
       sys.exit()
    la = len(A)
    lb = len(B)
    if la != lb:
       print ("Different length for A and B!")
       sys.exit()
#    C = []   
#    for i in range(la):
#	C.append(0)
#    for i in range(la):
#	C[i] = (B[i] - A[i]) * (m/n) + A[i]
    A = np.array(A)
    B = np.array(B)
    C = (B - A)*(m/n) + A
    C = C.tolist()



    return C



