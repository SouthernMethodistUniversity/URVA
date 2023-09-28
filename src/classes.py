import numpy as np


###############################
class xPathpoint:
      'Trajectory Point of XYZ file'
      pCount=0
 
      def __int__(self):
	  xPathpoint.pCount += 1
      NAtom = 0
      Vmass=[]
      Vcoor=[]
      VAnZ = []
      energy = 0.0


##############################


###############################
class nPathpoint:
     'Path Point of New Browsing File'
     pCount=0

     def __int__(self):
         nPathpoint.pCount += 1

     NAtom=0
     Vmass=[]
     Vcoor=[]
     Veta=[]
     Vkappa=[]
     s=0.0
     energy=0.0
#################################

###############################
class oPathpoint:
     'Path Point of Old Browsing File'
     pCount=0

     def __int__(self):
         oPathpoint.pCount += 1

     NAtom=0
     Vmass=[]
     Vcoor=[]
     Vgrad=[]
     Vhess=[]
     Veta=[]
     Vkappa=[]
     VAnZ=[]

     Vnm_0=[]

     s=0.0
     energy=0.0
#################################


#################################
class userInputs:
     'The class collecting all user input parameters.'
     pCount = 0

     def __int__(self):
	 userInputs.pCount += 1

     Datfiletype = -1               # Required
     Datfilename = "InputDataFile"  # Required
     Basename = "InputBaseName"     # Required 


     TFparm = 0

     TFvib    = 0

     TFcurvcorts = 0
     NcurvcorLn = 25
     NcurvcorRn = 25

     TFautosmth  = 0
     NautosmthLn = 3
     NautosmthRn = 3
     Nautosmthstpsize = -1         # Required
     Nautosmthd2y = 0.5

     TFremovespk = 0
     Nrmspkcut   = 20.0 
     Nrmspkperc  = 0.85
     Nrmspkgratio= 1.2


     TFdircurv = 0
     TFavam = 0
     TFcurvcpl = 0
     TFcoriolis = 0

     TFenergy = 1
     TFadiabfc = 0




#################################


def nBscalcurv(nBlist):
    results = []

    for i in range(len(nBlist)):
	sval = nBlist[i].s
        kappa0 = [item for sublist in nBlist[i].Vkappa for item in sublist]
	kappaval = np.linalg.norm( kappa0  )
        results.append( [sval, kappaval] )

    return results

#kappa0 = [item for sublist in Brlist[i].Vkappa for item in sublist] 




def oBscalcurv(oBlist):
    results=[]

    for i in range(len(oBlist)):
        sval = oBlist[i].s
        kappaval = np.linalg.norm( oBlist[i].Vkappa )
	#print oBlist[i].Vkappa
        results.append( [ sval, kappaval ] )

    return results




