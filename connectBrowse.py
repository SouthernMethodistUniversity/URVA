#!from optparse import OptionParser
#parser = OptionParser()


print "Warning: Be careful with the D and E format"
print   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
print   "@           Connecting Browsing File Tool (CBFT)             @"
print   "@                                                            @"
print   "@                        CATCO  Group                        @"
print   "@                         08/28/2015                         @"
print   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

RevFile="IRC.reverse"
ForFile="IRC.forward"
OutFile="IRC.browse"
RevHandle=open(RevFile,"r")
OutHandle=open(OutFile,"w")
RevLines=RevHandle.readlines()
RevBeginN=0
RevEndN=0
ForBeginN=0
ForEndN=0
RevBeginLN=[]
RevEndLN=[]
ForBeginLN=[]
ForEndLN=[]


for s in range(len(RevLines)):
     if RevLines[s]=="BEGIN\n":
         RevBeginN=RevBeginN+1
         RevBeginLN.append(s)
     if RevLines[s]=="END\n":
         RevEndN=RevEndN+1
         RevEndLN.append(s)

ForHandle=open(ForFile,"r")
ForLines=ForHandle.readlines()
for s in range(len(ForLines)):
     if ForLines[s]=="BEGIN\n":
         ForBeginN=ForBeginN+1
         ForBeginLN.append(s)
     if ForLines[s]=="END\n":
         ForEndN=ForEndN+1
         ForEndLN.append(s)

print "Forward Point # ",ForBeginN
print "Reverse Point # ",RevBeginN

if ForEndN!=ForBeginN:
         exit()
if RevEndN!=RevBeginN:
         exit()

#We need to modified the XXIRC value in reverse part


for m in range(RevEndN):
   sline= RevLines[RevEndLN[m]-1].split()
   xx= RevLines[RevEndLN[m]-1].split()[2]
   print type(xx)
   print xx
#   print float(xx)
#   kk=-1* float(xx)
   bb=float(xx)
   kk=-1*float(xx)
   ss=str(sline[0])+" "+str(sline[1])+" "+str(kk)+"\n"
   RevLines[RevEndLN[m]-1]=ss











for k in range(RevEndN):
    j=RevEndN - k - 1 
    for item in RevLines[RevBeginLN[j]:RevEndLN[j]+1]:
         OutHandle.write( "%s" % item ) 

for k in range(1,ForEndN):
   for item in  ForLines[ForBeginLN[k]:ForEndLN[k]+1]:
         OutHandle.write( "%s" % item )
 
print range(1,3)

for s in range(3):
   print ForBeginN


