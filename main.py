###########################################################
#   To Do List:
#   1. DMO
#     1.1 Add overlap check after each DMO step -> basically done
#     1.2 Fully test it -> basically done
#     1.3 Tested on CH3+H2 and PH2+H2 -> pass
#     1.4 Adjust User input 
#     1.5 Add skip point(s) option -> done  
#
#   2. Parameter filter ( not useful )
#   3. Ring Coordinates (interfacing done) 
#   4. Dihedral component result -> problematic ? -> check passed
#   5. Torsion of reaction path 
#   6. New general local mode force constant  -> 
#      6.a) new formula without projection -> implemented 
#      7.a) with projection -> lead to discontinuity
#
#   *. Recalculate the reaction parameter s value -> useless 
#   7. Process pure trajectory files (.xyz) 
#   8. Make the program read the user input  --> improving
#   9. Print parameter-related properties as csv files ... q_n, direction, curvature, -> done 
#  10. Put in all missing functionalities e.g. Curvature coupling coefficients, ...
#     #a) Coriolis coupling coefficient -> done 
#      b) Decomposition of normal modes into local modes -> Local mode program 
#      c) Force decomposition -> difficult to interpret   
#     #d) Curvature coupling coefficient -> done
#     #e) Adiabatic mode coupling coefficient -> done  
#
#  11. curvcor direct output option (for some small systems)
#  12. New browsing file  
#  13. Base name -> test not finished yet
#  14. Automatic correction of k^a along the reaction path   
#  15. Energy second derivative correction problem   
#
#  16. Spline fitting for all points on IRC path   
#
#
############################################################
#


import numpy as np
import math
import sys
import os
import copy
import time
import scipy

from sympy import *
from calc import * 
from classes import *
from fileio import *
from cologne import *

from interf_curvcor_main import *
from interfaces import *



logo()


########## User input ##########
#inpf="inp.inp"
inpf = sys.argv[1]
################################


cmdgui = userInputs()

usrinps = rdcmd(inpf,cmdgui)


#stop("pause after reading user input...")
#print usrinps


# Assign parameter from usrinps

oldnew = usrinps[0]


print("\nData file: "+usrinps[1])
print("Base folder: "+usrinps[22]+"\n")


numPn=0
## Initial scan of browsing file
print ("Info: Initial scanning of data file...")
if oldnew == 1:
   nBrowsf = usrinps[1] 
   numPn=iscNB(nBrowsf,numPn)

elif oldnew == 0:
   oBrowsf = usrinps[1] 
   numPn=iscOB(oBrowsf,numPn)

elif oldnew == 2:
   xyzf = usrinps[1] 
   numPn=iscXYZ(xyzf,numPn)

else:
   stop("Error: Please specify the data file type correctly!")


print ("Info: Total number of points in the data file: "+str(numPn))

## Create list
nBlist=[]
oBlist=[]
xBlist=[]


#stop("stopppppp")

print ("Info: Reading in data file...")
for i in range(numPn):
   if oldnew == 1:
      nBlist.append( nPathpoint() ) # append object
   elif oldnew == 0:
      oBlist.append( oPathpoint() ) # 
   else: # == 2
      xBlist.append( xPathpoint() )




## Read in browsing file
if oldnew == 1:
   rdNB(nBrowsf,nBlist)
   nBlist = checkduply( nBlist )
   numPn = len(nBlist)
elif oldnew == 0:
   rdOB(oBrowsf,oBlist)
   oBlist = checkduply( oBlist )
   numPn = len(oBlist)
else:
   rdXYZ(xyzf,xBlist)


# Check duplicate points on reaction path
## Delete one point if two points are identical in coordinates

#checkduply( Brlist )

#stop("check duplicate done...")


print ("Info: Data file loaded.")


###################################
# Fitting for the IRC points
####################################
if oldnew == 0 and 1 == 2:
   para,IRCpoints = GetPoints(oBlist)
  
   #print IRCpoints
   npts = len(para)
   ds = 0.02
   xs = np.arange( para[0],para[-1],ds/10 )

   yslib = []
   plib = []
   p1dlib = []
   p2dlib = []

#   print "XS:"
   
#   for i in xs:
#        print i


   for i in range(len(IRCpoints)): # 3*natom
       z = np.polyfit( para, IRCpoints[i] , 28 ) # Up to 16th-order polynomial fitting
       p = np.poly1d( z )
       #print "z",z
       #print "p",p
       p1 = np.polyder(p)
       #print "p1",p1
       p2 = np.polyder(p,m=2)
       #print "p2",p2

       plib.append( p )
       p1dlib.append( p1 )
       p2dlib.append( p2 )
       #print len(plib),len(p1dlib),len(p2dlib)



       ys = []
       for j in range(len( xs )):
           ys.append( p(xs[j]) )

#       print "YS:"
#       for j in ys:
#	   print j

       yslib.append( ys )

   #
   for i in range(npts):
       s = para[i]
       natom = len(IRCpoints)/3


       # Calculate r'(s) and evaluate it as well as its magnitude
       r1s = []
       for j in range(3*natom):
	   r1s.append( p1dlib[j](s) )
       #print r1s
       r1snorm = np.linalg.norm( r1s )

       # calculate T'(s) and its magnitude

       t1s = []
       for j in range(3*natom):
	   t1s.append( p2dlib[j](s) / r1snorm )

       t1snorm = np.linalg.norm( t1s )

       k = t1snorm / r1snorm 

       print k


       pass	   

   #




   

#





# Print Out Energy for 1) old browsing file and 2) new browsing file
## At the same time, first derivatives and second derivatives are also needed.
## Do we need to set up a switch for energy output?
### Current defaut: print all 


if oldnew == 1 and usrinps[20] == 1:
   PrtEnergy(nBlist, "energy.csv", "energy_1_d.csv", "energy_2_d.csv")	
   pass
elif oldnew == 0 and usrinps[20] == 1:
   PrtEnergy(oBlist, "energy.csv", "energy_1_d.csv", "energy_2_d.csv")
   pass



## Get the atom ID for parameters
if usrinps[2] >= 1:
   parID=[]
   parID=rdPar(inpf,parID)
   #print parID
   #stop("pause")

   ## Calculate the value of parameters
   print ("Info: Calculating q_n values...")
   parVal=[]
   if oldnew == 1:
      parVal=calparVal(numPn,nBlist,parID,usrinps[22])
      printqn2file( parVal, nBlist, "q_n.csv" )
   elif oldnew == 0:
      parVal=calparVal(numPn,oBlist,parID,usrinps[22])
      printqn2file( parVal, oBlist, "q_n.csv" )
   else:
      parVal=calparVal(numPn,xBlist,parID,usrinps[22])
      printqn2file( parVal, xBlist, "q_n.csv" )

   #print parVal
   #printqn2file has been checked and passed 


#stop("pause after print q_n ...   ")





## Calculate the B matrix
if usrinps[2] == 2:
   print ("Info: Calculating B matrices of q_n...")
   parB=[]
   if oldnew == 1:
      parB=calparB(numPn,nBlist,parID,usrinps[22])
   elif oldnew == 0:
      parB=calparB(numPn,oBlist,parID,usrinps[22])
   else:
      print("Warning: B matrices are not going to be calculated for XYZ data file. Skipping...")

   #print parB
   print ("Info: B matrix calculation completed.")

#stop("pause after calculating B matrices...")



#------------------------------------------
## Calculate curvature for NEW browsing file
if oldnew == 1:
   if usrinps[16] == 1:
      pass
      origcorprt2(nBlist,"originalkappa.dat")
      print("Output: Scalar curvature written to "+str("originalkappa.dat"))




## Calculate curvature and path direction for OLD browsing file
if oldnew == 0 :# and usrinps[16] == 1:
   pass
   dircur(oBlist) #

#	find atoms with largest contribution in curvature vector -> debug use only  
#	findcontri(oBlist,1.6,5, "CurvAtoms.dat" )


   if usrinps[16] == 1:
      origcorprt(oBlist,"originalkappa.dat")
      print("Output: Scalar curvature written to "+str("originalkappa.dat"))  # just for old browsing file


   # CURVCOR module 
   if usrinps[4] == 1:
      runcorvcor(usrinps[5],usrinps[6],oBlist,"correctedkappa.dat","originalkappa.dat")

      # AUTOSMTH module
      if usrinps[7] == 1:
        runautosmth("originalkappa.dat","correctedkappa.dat","corrsmthkappa.dat","merged.dat",usrinps[10],usrinps[8],usrinps[9],usrinps[11],usrinps[22])
        print("Output: Scalar curvature written to "+str("merged.dat"))

        # RMSPK module 
        if usrinps[12] == 1:
           runrmspk("merged.dat","merged-nospk.dat",usrinps[13],usrinps[14],usrinps[15],usrinps[22])
	   print("Output: Scalar curvature written to "+str("merged-nospk.dat"))
        # By using indent, impose conditions implicitly


## Calculate the components of Path direction and curvature
#  Input information:
#   1. oldnew
#   2. Veta
#   3. Vkappa
#   4. numPn
#   5. parB 
#   6. s  
#   7. parID
#
#
#  Output file:
#   1. eta-q_n.csv
#   2. kappa-q_n.csv
#
#

# caldircurcomp(oldnew, Brlist, numPn, parB, parID, "eta-q_n.csv", "kappa-q_n.csv" ) 


if usrinps[2] == 2 and usrinps[16] == 1:
#
   if oldnew == 0:
      qnetadata,qnkappadata = caldircurcomp( oldnew, oBlist, numPn, parB )
   if oldnew == 1: 
      qnetadata,qnkappadata = caldircurcomp( oldnew, nBlist, numPn, parB )

   printqncomp(qnetadata,parID,   "eta-q_n.csv",1)
   printqncomp(qnkappadata,parID, "kappa-q_n.csv",2)


#stop("pause after decomposition...")



# 




## Calculate Adiabatic force constant -(New formula)
# for old browsing file onlu
# Input: 1. B matrix -> parB 
#        2. Hessian matrix -> oBlist.Vhess 
#

if oldnew == 0 and usrinps[2] == 2 and usrinps[21] == 1:
   pass
   #print parB
   print ("Info: calculating adiabatic force constant...")
   
   calKaNew( oBlist, parB,parID ,"adiabfc-ka.csv")

   #calKa2p( oBlist,parB ) 

#   stop("pause before calculating k^a...")



## Calculate Adiabatic mode coupling coefficient (AvAM)
# for old browsing file only 
# Input: 1. Curvature vector -> oBlist
#        2. Mass -> oBlist     
#        3. Generalized adiabatic modes 
#           3.1 Hessian -> oBlist
#           3.2 B matrix -> parB
#           3.3 Cartesian coordiantes -> oBlist
#           3.4 Mass -> oBlist
#
# This part relies on reaction path direction and curvature


if oldnew == 0 and usrinps[17] == 1 and usrinps[2] == 2: 
   pass
   calAvAM( oBlist, parB,parID,"adiab-coupling.csv" )  # Output file not ready for AvAM yet

#stop("Pause after AvAM...")


#-----------------------------------------------------------------------------------------------------------
## Calculate vibrational frequencies if requested
#
# CAUTION: Abortion in IRC might lead to failure of DMO, so this situation should be considered.      
#          The abortion in IRC could be identified via overlap checking of path direction between two points.      
#         
#          Solutions to cure this problem -> requires further testing...
#                1. Check frequency value (find the closest partner)
#                2. If (1.) fails, ask the user for an internal coordinates parameter, then use the product between      
#                   the B matrix of this parameter and vector(s) of normal modes to decide the closest partner   
#                3. ...
#
#  Notes on DMO: An overlap check after each DMO step is necessary. This is could be used as an indicator of local difficulty  
#                Smaller step should be employed for such situation. (e.g. PH2 + H2) 
#
#
#  Steps of calculation:
#   A. Store frequencies and normal modes for every point on the path
#
#   B. Just calculate 2 consecutive points and do the DMO right away  
#      1. Take less memory, but more computer time
#      2. Find error in DMO early
#      3. Report crossings in time 
#     *4. Save ordered normal modes for later use -> not necessary  
#      5. Easier to analyze just part of data file
#
#  Further considerations for normal mode -dependent properties   
#      1. Whether or not to store the normal mode vectors -> current strategy is NO. Calculate them on the fly.
#      2.   
#

if oldnew == 0 and usrinps[3] == 1:
   pass
   print("Info: Calculating vibrational frequencies...")
   
   ####################################GenTRQ(oBlist[0].Vcoor,oBlist[0].Vmass)

   # Get the index for the TS point
   for i in range(len(oBlist)):
      NoDMO = 0

      if NoDMO == 1:
	  e,nm,cartdisp=solWils(oBlist[i].Vcoor,oBlist[i].Vhess,oBlist[i].Vmass)
	  print("s: "+str(oBlist[i].s))
          prtVib(e)
          if i > 3:
	     stop(str(i+1)+" vibs.")

      if oBlist[i].s == 0.0: 
          k = i

          #k = 563
          kpool = [k -2, k -1, k, k +1, k +2]
          
	  for pp in kpool:
              e,nm,cartdisp=solWils(oBlist[pp].Vcoor,oBlist[pp].Vhess,oBlist[pp].Vmass)
              oBlist[pp].Vnm_0 = nm[0]

     #     e,nm,cartdisp=solWils(oBlist[k].Vcoor,oBlist[k].Vhess,oBlist[k].Vmass)
#         print "With Eckart conditions"
#         prtVib(e)
     #     oBlist[k].Vnm_0 = nm[0]




   
          #stop("pause after vibration analysis...")

#         print ("With projection of nm[0]... ")
          e2,nm2,cartdisp2=solWilsProj( oBlist[k].Vcoor,oBlist[k].Vhess,oBlist[k].Vmass,nm[0] )

          #prtVib(e2)


   # Check whether the path is continuous or not.
   # store the information of the location of discontinuities 

   if oldnew == 1:
      disctlist = chkcont(nBlist,1)

   elif oldnew == 0:
      disctlist = chkcont(oBlist,0)

   #print disctlist

   # enter DMO engine and get reordered frequencies as well as normal modes 
   calfrq_dmo(oBlist,disctlist,usrinps[18],usrinps[19],inpf,"freq_dmo.csv")







stop("\nNormal termination at "+time.strftime("%c"))













   #print oBlist[k].Vnm_0




   #
   #e3,nm3,cartdisp3=solWilsProj( oBlist[1].Vcoor,oBlist[1].Vhess,oBlist[1].Vmass,oBlist[1].Veta )
   #prtVib(e3,cartdisp3)
   #print " "
   #


   # Print out the gradient norm 
   # 
   #for j in range(numPn):
   #    print (  '%s %7.4f %s %10.7f %s '   %  ("s:",oBlist[j].s, "force:", np.linalg.norm(oBlist[j].Vgrad), " a.u."  ) )
   ## print oBlist[j].Vkappa
   ## print (len(oBlist))



######### TEST DMO  ############

#e1,nm1,cartdisp1=solWilsProj( oBlist[0].Vcoor,oBlist[0].Vhess,oBlist[0].Vmass,oBlist[0].Veta )
#e2,nm2,cartdisp2=solWilsProj( oBlist[1].Vcoor,oBlist[1].Vhess,oBlist[1].Vmass,oBlist[1].Veta )
#e3,nm3,cartdisp3=solWilsProj( oBlist[2].Vcoor,oBlist[2].Vhess,oBlist[2].Vmass,oBlist[2].Veta )
#e4,nm4,cartdisp4=solWilsProj( oBlist[3].Vcoor,oBlist[3].Vhess,oBlist[3].Vmass,oBlist[3].Veta )
#e5,nm5,cartdisp5=solWilsProj( oBlist[4].Vcoor,oBlist[4].Vhess,oBlist[4].Vmass,oBlist[4].Veta )
#e6,nm6,cartdisp6=solWilsProj( oBlist[5].Vcoor,oBlist[5].Vhess,oBlist[5].Vmass,oBlist[5].Veta )
#e7,nm7,cartdisp7=solWilsProj( oBlist[6].Vcoor,oBlist[6].Vhess,oBlist[6].Vmass,oBlist[6].Veta )
#e8,nm8,cartdisp8=solWilsProj( oBlist[7].Vcoor,oBlist[7].Vhess,oBlist[7].Vmass,oBlist[7].Veta )
#e9,nm9,cartdisp9=solWilsProj( oBlist[8].Vcoor,oBlist[8].Vhess,oBlist[8].Vmass,oBlist[8].Veta )
#
#e,nm,cartdisp=solWils(oBlist[9].Vcoor,oBlist[9].Vhess,oBlist[9].Vmass)
#e10,nm10,cartdisp10=solWilsProj( oBlist[9].Vcoor,oBlist[9].Vhess,oBlist[9].Vmass,oBlist[9].Veta )

#e11,nm11,cartdisp11=solWilsProj( oBlist[10].Vcoor,oBlist[10].Vhess,oBlist[10].Vmass,oBlist[10].Veta )

#print nm2
#print e2

#nm2,e2 = DMODriv(nm1,nm2,e2,18,11)
#nm3,e3 = DMODriv(nm2,nm3,e3,18,11)
#nm4,e4 = DMODriv(nm3,nm4,e4,18,11)
#nm5,e5 = DMODriv(nm4,nm5,e5,18,11)
#nm6,e6 = DMODriv(nm5,nm6,e6,18,11)
#nm7,e7 = DMODriv(nm6,nm7,e7,18,11)
#nm8,e8 = DMODriv(nm7,nm8,e8,18,11)
#nm9,e9 = DMODriv(nm8,nm9,e9,18,11)
#nm10,e10=DMODriv(nm9,nm10,e10,18,11)
#nm11,e11=DMODriv(nm10,nm11,e11,18,11)

#prtVib(e1,nm1)
#prtVib(e2,nm2)
#prtVib(e3,nm3)
#prtVib(e4,nm4)
#prtVib(e5,nm5)
#prtVib(e6,nm6)
#prtVib(e7,nm7)
#prtVib(e8,nm8)
#prtVib(e9,nm9)
#prtVib(e10,nm10)
#prtVib(e11,nm11)
#
#for i in range(11):
    #for j in range(11):
#	 prd = np.dot( nm4[i], nm5[i] )
#	 print i+1,i+1,prd




# Solve local difficulty in PH2 + H2 between 
# -0.509996648415 ->  9
# -0.504996917815 -> 10
print ("Solving local difficulty for PH2 + H2...")
e9,nm9,cartdisp9=solWilsProj( oBlist[9].Vcoor,oBlist[9].Vhess,oBlist[9].Vmass,oBlist[9].Veta )
e10,nm10,cartdisp10=solWilsProj( oBlist[10].Vcoor,oBlist[10].Vhess,oBlist[10].Vmass,oBlist[10].Veta )

nstep = 0

#prtVib(e10)

aE = e9
aNM = nm9

for i in range(nstep):
    j = i + 1
    bE,bNM,bCAR=solWilsProj( \
		             ABscIn(j,nstep,oBlist[9].Vcoor,oBlist[10].Vcoor),\
                             ABscIn(j,nstep,oBlist[9].Vhess,oBlist[10].Vhess),\
                             oBlist[9].Vmass,\
                             ABscIn(j,nstep,   oBlist[9].Veta,oBlist[10].Veta) \
			     )
    nvib = len(bNM)
    natom3 = len(bNM[0])

    bNM,bE = DMODriv(aNM,bNM,bE,natom3,nvib)
    tmpL = []
    for k in range(nvib):
        tmpL.append( round(np.dot(aNM[k],bNM[k]),4) )

    #print i,nstep,tmpL
    print ABscIn(j,nstep,[oBlist[9].s],[oBlist[10].s]),prtViblist(bE)

    aE = bE
    aNM = bNM







### Solve the local difficulty between e4 and e5
#e4,nm4,cartdisp4=solWilsProj( oBlist[3].Vcoor,oBlist[3].Vhess,oBlist[3].Vmass,oBlist[3].Veta )
#e5,nm5,cartdisp5=solWilsProj( oBlist[4].Vcoor,oBlist[4].Vhess,oBlist[4].Vmass,oBlist[4].Veta )
#print type(oBlist[3].Vcoor)


#K1=[0,0]
#K2=[1,2]
#print ABscIn(2,2,K1,K2)
#nstep = 10

#aE=e4
#aNM=nm4
#for i in range(1,nstep+1):
      #print i
#      bE,bNM,bCAR=solWilsProj( \
#		               ABscIn(i,nstep,oBlist[3].Vcoor,oBlist[4].Vcoor),\
#		               ABscIn(i,nstep,oBlist[3].Vhess,oBlist[4].Vhess),\
#		               oBlist[3].Vmass,\
#		               ABscIn(i,nstep,   oBlist[3].Veta,oBlist[4].Veta) \
#			       )

 #     bNM,bE = DMODriv(aNM,bNM,bE,18,11)
 #     tmpL=[]
 #    for j in range(11):
#	  tmpL.append( round(np.dot(aNM[j],bNM[j]),4) )
  #    print i,nstep,tmpL
#      aE = bE
#      aNM = bNM
	  

if 2 == 1:
 print ("testing basic DMO...")
 a=[[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0]]
 b=[[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[1,0,0,0,0,0]]

 #a=[[1,0,0,0,0,0],[0,1,0,0,0,0]]
 #b=[[0,1,0,0,0,0],[1,0,0,0,0,0]]

 frq=[2.5,3.5,4.5,1.5]
 b,frq = DMODriv(a,b,frq,6,4)
 print frq
 print b


