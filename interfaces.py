import sys
import os
from interf_curvcor_main import *
from fileio import expList2File
from classes import oBscalcurv
from classes import nBscalcurv
from fileio import smergel


   ## This module migrates the curvature correction, smoothening and spike removal
   ## interfaces from main.py into this place.

   ## Correct the curvature near the TS point
    # 1. calculate the corrected curvature value first
    # 1.1 modify current script and make it work on objects 
    # 1.2 Export these corrected curvature result as an external file
    #
    # NOTICE: It is recommended to check the overlap between two directions   
    # WARNING: If the IRC calculation aborts very close the TS, the curvature could hardly
    #          be corrected using this approach
    #
    # 2. smoothen the curvature plot   
    # 2.1 Read in two data text file
    # 2.2 Run the script
    # 2.3 No need to modify this script
    #



def origcorprt2(Brlist,outfile):
    origkappa = nBscalcurv(Brlist)
    expList2File( origkappa, outfile)
    # For NEW browsing file



def origcorprt(Brlist,outfile):
    origkappa = oBscalcurv(Brlist)
    expList2File( origkappa, outfile )

def runcorvcor(Ln,Rn,Brlist,outfile1,outfile2):
    
    corrkappa = corcurv(Brlist,Ln,Rn)   
    expList2File( corrkappa, outfile1 )

    #origkappa = oBscalcurv(Brlist)
    #expList2File( origkappa, outfile2 )




#   # CURVCOR module 
#   if usrinps[4] == 1:
#
#    Ln = usrinps[5]   #25
#    Rn = usrinps[6]   #25
#    corrkappa = corcurv(oBlist,Ln,Rn)
#    expList2File( corrkappa, "correctedkappa.dat" )
#
#   #stop("debug corrected kappa...")
#
#   # collect original curvature
#    origkappa = oBscalcurv(oBlist)
#   # print origkappa     
#    expList2File( origkappa, "originalkappa.dat" )
#


def runautosmth(infile1,infile2,outfile,outfile2,stepsize,ln,rn,d2ythresh,basepath):
      infile1 = infile1 
      infile2 = infile2 
      outfile = outfile 
      smthrunstr = "python "+basepath+"inter_autosmooth_main.py "+ infile1 +" " + infile2 + " "+ \
		    str(stepsize)+" "+str(ln)+" "+str(rn)+" "+\
		    str(d2ythresh)+" "+outfile+" "

      os.system( smthrunstr )

      smergel(outfile,infile1,outfile2)

   # Start the automatic smoothening
   ## parameters required
   ## 1. original curvature file
   ## 2. corrected curvature file
   ## 3. step size
   ## 4. Number of points for cubic spline fitting (Left, Right)    
   ## 5. minimum *curvature* threshold value      
   ## 6. output file 
   #
   #  Mechanism of working: os.system("...")

#    smthorigpath = "originalkappa.dat"+" "
#    smthcorrpath = "correctedkappa.dat"+" "
#    smthoutfpath = "corrsmthkappa.dat"+" "
#    smthstepsize =  usrinps[10] #0.03 
#    smthleftn = usrinps[8] #4
#    smthrighn = usrinps[9] #4
#    smthcurthre = usrinps[11]  #0.5
#
#    smthrunstr = "python inter_autosmooth_main.py "+smthorigpath+smthcorrpath+\
#                str(smthstepsize)+" "+str(smthleftn)+" "+str(smthrighn)+" "+\
#		str(smthcurthre)+" "+smthoutfpath
#
#   os.system( smthrunstr ) 
#
#
##    # Merge the corrected curvature curve around TS into original curvature curve
#    smergel("corrsmthkappa.dat","originalkappa.dat","merged.dat")




def runrmspk(infile,outfile,cuthigh,percent,gradratio,basepath):
     rmspkrunstr = "python "+basepath+"inter_removespikes_main.py " + infile + " " + str(cuthigh)+" "+\
                   str(percent) + " " + str(gradratio) +" "+outfile
      
     os.system( rmspkrunstr )


#    rmspkinpath = "merged.dat" + " "
#    rmspkotpath = "merged-nospk.dat"+ " " 
#    rmspkcuthigh= usrinps[13]  #0.8 # normally 20
#    rmspkpercent= usrinps[14]  #0.90 # range from 0.85 to 0.98
#    rmspkgratio = usrinps[15]  #1.2
#  
#    rmspkrunstr = "python inter_removespikes_main.py "+rmspkinpath + str(rmspkcuthigh)+" "+\
#		  str(rmspkpercent)+" "+str(rmspkgratio)+" "+rmspkotpath
#   
#    os.system( rmspkrunstr ) 







