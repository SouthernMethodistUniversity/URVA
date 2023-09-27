# This is the library for 
#

import sys
from scipy import optimize
import math
from scipy.interpolate import interp1d
import numpy as np



def diffe(x,Snew,Vnew,Vold,dstep): # Here x is the variable and others are all parameters
    # calculate the difference between the scaled value and original value(calculated by old formula)
    if abs(Vnew) > abs(Vold):
       sign0 = -1 # scale down
    else:
       sign0 = 1  # scale up
    #print x
    k = math.exp(sign0*(abs(Snew)/dstep)/x)   
    err = abs(k*Vnew - Vold) 

    return err
#

 
def graddiff(x,ps,js,qs,pv,jv,jo,qv,dstep):
    # This function returns the direction difference (with regard to gradient)
    # ps -> s value for the point before that in question -> kp, pv
    # js -> s value for the point in question -> kj, jv
    #    jo -> value of the point in question on old curve 
    # qs -> s value for the point after that in question (old curve)
    #
    ### When pv == jv or jo == qv, exceptions happen! 

    if abs(jv) > abs(jo):
       sign0 = -1
    else:
       sign0 =  1
    kp = math.exp(sign0*(abs(ps)/dstep)/x) 
    kj = math.exp(sign0*(abs(js)/dstep)/x)
    g1 = (kj*jv - kp*pv) / (js - ps )
    g2 = (qv    - jo   ) / (qs - js )
    #print "derviation:",(kj*jv - jo)
    fx = jo     # -> current point - old 
    fp = qv     # -> past point - old
    fn = kp*pv  # -> before point - scaled 
    f2d = (fp - 2*fx + fn)/(dstep*dstep)
    #print "curvature",f2d
    if qv == jo:
       g2 = 1.0

    return (g2/g1,f2d)
    # Return also the curvature value as "f2d"

def detercutB(curvaturelist,thresh):
    # Determine when to cut via looking at the curvature value
    # The curvature characterize the change of gradients
    #
 
    gi  = -1
    log = []
    for i in range(len(curvaturelist)):
        log.append(abs(curvaturelist[i]))

    if min(log) > thresh:
       print ("Smoothening with regard to the 2nd derivative check failed...")
       print ("2nd derivative "+str(min(log))+" exceeding threshold of "+str(thresh)      ) 
       sys.exit()
    else:
       gi = log.index(min(log))
    return gi+1


def detercut(ratiolist,thresh):
    # determine where to cut for the scaling part 
    gi = -1
    log = []
    for i in range(len(ratiolist)):
        if ratiolist[i] > 0:
           if ratiolist[i] < 1:
              ratiolist[i] = 1.0 / ratiolist[i] 
    #print ratiolist  
    for i in range(len(ratiolist)):
        if ratiolist[i] >= 1.0:
           log.append( ratiolist[i] )
    #print log
    if  min(log) > thresh:
        print ("Smoothening with regard to gradient failed...")
        print ("Gradient ratio "+str(min(log))+" exceeding threshold of "+str(thresh)  )
        sys.exit()
    else:
        gi = ratiolist.index( min(log) )

    return gi+1
 


 
#############################################
# Following is an example which shows clearly 
# how to let Python to do optimization
#############################################
# from scipy import optimize
#
# def f(x,a,b):
#    return a*(x-1-a-b)**2
#
# a=10
# b=100
# res = optimize.minimize_scalar(f,args=(a,b))
# print res.x
############################################## 




def spline1dlist(lis1):
    x=[]
    y=[]
    for i in range(len(lis1)):
        x.append( lis1[i][0] )
        y.append( lis1[i][1] )
    f = interp1d(x,y,kind='cubic')

    return f 

def get1point(data,marks):
   pts=[]
   for i in range(len(marks)):
       if marks[i] == 1:
          pts.append(data[i])
          break
   return pts


def get2points(data,marks):
    # This function returns the two points to be predicted by spline fitting
    pts=[]
    for i in range(len(marks)):
        if marks[i] == 1:
           pts.append(data[i])
           break
    for i in range(len(marks)):
        j = len(marks) - i - 1
        if marks[j] == 1:
           pts.append(data[j])
           break

    return pts     
 
def checkerr(f,pts):
    # check the error between predicted value and train value
    x=[]
    y0=[]
    y1=[]
    err=[]
    for i in range(len(pts)):
        x.append(pts[i][0] )
        y0.append(pts[i][1])
    y1 = f(x)
    #print y0
    #print y1
    for i in range(len(y0)):
        err.append( abs(y0[i] - y1[i]) )
    print err
    return err

def adjustmark(marks,l,r):
    # this function could make 1 into 0
    if l == 1:
     for i in range(len(marks)):
        if marks[i] == 1:
            marks[i] = 0
            break
    if r == 1:
     for i in range(len(marks)):
        j = len(marks) - i - 1 
        if marks[j] == 1:
            marks[j] = 0
            break 
    return marks


def updateval( err,thresh,data,pts,f ):

    x=[]
    replaceflag = []
    for i in range(len(err)):
        if err[i] > thresh:
           replaceflag.append(1) 
        else:
           replaceflag.append(0) # modified!!!!

    if len(pts) == 1:
       for i in range(len(data)):
           if pts[0][0] == data[i][0]:
              if replaceflag[0] == 1:
                 #data[i][1] =
                 x.append( pts[0][0] )
                 m = f( x )
                 data[i][1] = m[0]
                 #print "w",f( pts[0][0] ) 
                 print "1:got it!"

    if len(pts) == 2:
       for i in range(len(pts)):
           for j in range(len(data)):
               if pts[i][0] == data[j][0]:
                  if replaceflag[i] == 1: 
                     #data[j][1] =
                     x.append( pts[i][0] ) 
                     m = f( x )
                     #print "m",type(m) ,m[0],type(m[0])
                     data[j][1] = m[0]    # f( pts[i][0] ) 
                     print "2:got it!",pts[i][0]
    return data              




 # SPLINE FITTING 
def itersplif(data,marks):   
    if len(data)!=len(marks):
       print ("Length not match...")
       sys.exit()
    
    if sum(marks) >= 2:
       run = 1

    while run == 1: 
     if sum(marks) >= 2:
       #print "not last" 
       train=[]
       for i in range(len(marks)):
           if marks[i] == 0:
              train.append( data[i] )

       f = spline1dlist(train) 
       pts = get2points(data,marks)
       err = checkerr(f,pts) # Set threshold before here.

       data = updateval( err,prederrthresh,data,pts,f )
       
       # If error larger than threshold, take the spline fitting value,
       # Otherwise, stay the same.
       # Remember to update the "data" when necessary

       #print marks
       marks = adjustmark(marks,1,1) 
       #print marks
       #run = 0 
       if sum(marks) == 0:
          run = 0

     if sum(marks) == 1:
       #print "last" 
       train=[]
       for i in range(len(marks)):
           if marks[i] == 0:
              train.append( data[i] )
       f = spline1dlist(train)
       pts = get1point(data,marks)
       err = checkerr(f,pts) 
       run = 0
       # This last point, we just take the spline fitting value!
       data = updateval( err,prederrthresh,data,pts,f )
       marks = adjustmark(marks,1,0)
        
    return data      



def semispline( forsf, RL ):
    # Semi-automatic spline fitting 

    for i in range(len(forsf)):
      if forsf[i][0] == 0.0:
         TSi = i
    right = RL[1]
    left =  RL[0]

    # Points that need to be predicted    
    excp=[]
    excp.append( TSi )
    for i in range(right):
        j = TSi + i + 1
        excp.append(j)
    for i in range(left):
        j = TSi - i - 1
        excp.append(j)

    #print excp
    train=[]
    preds =[]
    for k in range(len(forsf)):
        if k in excp:
           preds.append(forsf[k][0])
        else:
           train.append(forsf[k])
    #print train
    f = spline1dlist( train )
    
    for i in range(len(forsf)):
        if forsf[i][0] in preds:
           x =[]
           x.append(forsf[i][0])
           y = f(x)
           forsf[i][1] = y[0]  
        
    return forsf      
 

        
