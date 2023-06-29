print("Validata RMM matrix") 

# do you want to remove 9 Mjj masses?
isRemoveMjj=True
print("Remove Mjj=",isRemoveMjj)


import os,sys
if sys.version_info >= (3, 0):
    sys.stdout.write("Sorry, requires Python 2.x, not Python 3.x\n")
    sys.exit(1)


# Max events to look at
MaxEvents=100000

from ROOT import gROOT,gPad,gStyle,TCanvas,TSpline3,TFile,TLine,TLatex,TAxis,TLegend,TPostScript
from ROOT import TH2D,TF2,TArrow,TCut,TPad,TH1D,TF1,TObject,TPaveText,TGraph,TGraphErrors,TGraphAsymmErrors
from ROOT import TGraph2D,TTree,TMultiGraph,TBranch,gSystem,gDirectory
from ROOT import TPaveStats,TProfile2D 
import math
import ROOT
import sys
import numpy
import numpy as np
#from  global_module import *


#sys.path.append("nnet/")
#import bpnn

print ('Number of arguments:', len(sys.argv), 'arguments.') 
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)') 

channel=0
if (len(sys.argv) == 3):
   channel = sys.argv[1]
   myinput = sys.argv[2]

myinput="interactive"
if (len(sys.argv) ==2):
    channel = sys.argv[1]

print ("Mode=",myinput)
print("Validata RMM matrix for channel = ",channel)

out=[0]
gROOT.Reset()
figdir="figs/"
fname=os.path.basename(__file__)
epsfig=figdir+fname.replace(".py",".eps")
epsfig=figdir+fname.replace("data","data_"+str(channel))


name="projection"
#name="profile"
nameX=""
nameY=""
Ymin=0.0
Xmin=0
Xmax=6.0 
ZSmin=0.000001
ZSmax=0.5 

#from array import array
#colors = [0, 1, 2, 3, 4, 5, 6]; # #colors >= #levels - 1
#s = array('i', colors)
#gStyle.SetPalette(s, len(colors));
gStyle.SetNumberContours(500)
gStyle.SetPalette(1)

######################################################
gROOT.SetStyle("Plain");
gStyle.SetLabelSize(0.035,"xyz");
c1=TCanvas("c_massjj","BPRE",10,10,800,700);
c1.Divide(1,1,0.008,0.007);
ps1 = TPostScript( epsfig,113)
c1.SetGrid();

c1.cd(1);
gPad.SetLogy(0)
gPad.SetLogz(1)
gPad.SetTopMargin(0.05)
gPad.SetBottomMargin(0.1)
gPad.SetLeftMargin(0.1)
gPad.SetRightMargin(0.15)
gPad.SetGrid();
gPad.SetTicks(2,0)

# test data
# input data 
# proc=["analysis/out/sys0/mc_wzprime2021_anomaly1/801125.root","analysis/out/sys0/mc_wzprime2021/801125.root"]
# proc=["analysis/out/sys0/data/data_2015_2016.root"]

# proc=["analysis/out/sys0/mc_hplus2021/450004.root"]
# data
proc=["out/out_pythia8_sm14tev.root"]


rfile=[]
print ("Look at ", len(proc) ," processes using RMM", proc) 

for i in proc:
     rfile.append(ROOT.TFile.Open(i))
     print(i) 

dimensions=(rfile[0]).Get("dimensions");
h_cpucores=(rfile[0]).Get("cpucores");
print("Nr of used CPU=")
cpucores=int(h_cpucores.GetBinContent(2))
print("Nr of CPU used=",cpucores)

dimensions.Print("All")

maxNumber=int(dimensions.GetBinContent(2)/cpucores)
maxTypes=int(dimensions.GetBinContent(3)/cpucores)
mSize=int(dimensions.GetBinContent(4)/cpucores)
print ("maxNumber=",maxNumber," maxTypes=",maxTypes," mSize=",mSize)
mSize=maxTypes*maxNumber+1;
hhD = TProfile2D("profile", "profile", mSize, 0, mSize, mSize, 0, mSize, 0, 10);
names=["MET","j", "b", "#mu", "e", "#gamma"]


# dijet invariant mass
x=1+0*maxNumber+1  # X position  
y=1+0*maxNumber    # Y position 
mjj=(x,y) #  index of Mjj  matrix ellement 

# PT of first jet
x=1+0*maxNumber  # X position  
y=1+0*maxNumber  # Y position 
pT=(x,y) #  index of Mjj  matrix ellement 

#  bb mass 
x=1+1*maxNumber+1
y=1+1*maxNumber
mbb=(x,y)

#  bj mass 
x=1+1*maxNumber
y=1+0*maxNumber
mbj=(x,y)

# mu+mu 
x=1+2*maxNumber+1
y=1+2*maxNumber
mmumu=(x,y)

# e+e 
x=1+3*maxNumber+1
y=1+3*maxNumber
mee=(x,y)

# j+mu 
x=1+2*maxNumber
y=1+0*maxNumber
mjmu=(x,y)

# j+e 
x=1+3*maxNumber
y=1+0*maxNumber
mje=(x,y)

# j+gamma 
x=1+4*maxNumber
y=1+0*maxNumber
mjg=(x,y)

# b+mu 
x=1+2*maxNumber
y=1+1*maxNumber
mbmu=(x,y)

# b+e  
x=1+3*maxNumber
y=1+1*maxNumber
mbe=(x,y)

# b+gamma 
x=1+4*maxNumber
y=1+1*maxNumber
mbg=(x,y)


hmass=TH1D("Mass","mass",100,0,7000)

Names1=[]
Names1.append(names[0]);
for h in range(1,maxTypes+1,1):
       for i in range(1,maxNumber+1):
                 Names1.append(names[h]+str(i));
Names2=[]
for i in range(len(Names1)):
         Names2.append(Names1[i]);
Names1= Names1[::-1]

print(Names1) 

# first make empty matrix with the labels
for h in range(mSize):
      for w in range(mSize):
        i1=h;
        i2=w;
        hhD.Fill(Names2[i1],  Names1[i2], 0.0);
        #hhD.Fill(i1,  i2, 0.0);

kk=0
ntot=0
active=0
MaxTrim=5
zeroMatrix = numpy.zeros(shape=(mSize,mSize))
for h in range(mSize):
      for w in range(mSize):
        i1=h;
        i2=w;
        val=1.0 
        gamma=(maxTypes-1)*maxNumber
        electrons=(maxTypes-2)*maxNumber
        muons=(maxTypes-3)*maxNumber

        # vertical removals 
        if (h>gamma+MaxTrim): val=0 # photons
        if (h>electrons+MaxTrim and h<gamma+1): val=0 # muons
        if (h>muons+MaxTrim and h<electrons+1): val=0 # electrons 
 
        # horisontal removals 
        if (w<MaxTrim): val=0 # photons
        if (w>MaxTrim*2-1 and w<MaxTrim*3): val=0 # muons
        if (w>MaxTrim*4-1 and w<MaxTrim*5): val=0 # electrons 

        ntot=ntot+1
	if (val==1): 
               active=active+1

        #k1=w;
        #k2=mSize-h-1;

        #k1=mSize-w;
        #k2=mSize-h-1;

        k1=mSize-w-1;
        k2=h;

        # mjj
        if (isRemoveMjj==True):
          if (k1==mjj[0] and k2==mjj[1]): val=0
          if (k1==mbj[0] and k2==mbj[1]): val=0
          if (k1==mje[0] and k2==mje[1]): val=0
          if (k1==mjmu[0] and k2==mjmu[1]): val=0
          if (k1==mjg[0] and k2==mjg[1]): val=0
          if (k1==mbb[0] and k2==mbb[1]): val=0
          if (k1==mbe[0] and k2==mbe[1]): val=0
          if (k1==mbmu[0] and k2==mbmu[1]): val=0
          if (k1==mbg[0] and k2==mbg[1]): val=0
 
        # fill 0 matrix if there are values > 0 
        if (val>0): 
             zeroMatrix[k2][k1] = val
             hhD.Fill(Names1[w],  Names2[h], val);
             #hhD.Fill(i1,  i2, val);


print("number of active=",active)

# now create list with indexes for removal (0 trimming) 
idx2remove=[]
zeroRMM=(zeroMatrix.flatten()).tolist()
for k in range(len( zeroRMM )):
             if (zeroRMM[k]==0):
                   idx2remove.append(k+1)

print("Start validation:")
evt=0
for event in rfile[0].inputNN:
       NN=(event.proj).size()
       a=event.proj
       inx1=event.proj_index1
       inx2=event.proj_index2
       Trun = event.run
       Tevent=event.event
       #Tweight=event.weight # for MC with weigths
       weight=1.0; # Tweight;
       rmmMatrix = numpy.zeros(shape=(mSize,mSize))
       for i3 in range(NN):
              w=inx1[i3];
              h=inx2[i3];
              i1=w;
              i2=mSize-h-1;
              val=float(a[i3])
              rmmMatrix[w][h] = val
              #print Names2[h],  Names1[w],float(val)
              #hhD.Fill(Names1[i2],  Names2[i1], val); 
              #hhD.Fill(i1,  i2, val);
       evt=evt+1
       if (evt%10==0): print("Event=",evt)
       if (evt>MaxEvents): break

      
       # validate to make sure we do not remove actual data!!
       dataRMM=(rmmMatrix.flatten()).tolist()
       for k in range(len( dataRMM )):
                          x=k/mSize
                          y=k%mSize
                          i1=x;
                          i2=mSize-y-1;
                          val=float(dataRMM[k])
                          isremove=False
                          if (val>0):
                            if k+1 in idx2remove:
                                     hmass.Fill(val*13000);
                                     #print("Mass=",val*13000);
                                     #print("Error = data removed!",k+1, "value=",val, " i1=",i1," i2=",i2) 
                                     rmmMatrix[w][h] = val
                                     isremove=True
                          if (isremove==False):
                                hhD.Fill(Names1[i2],  Names2[i1], val);
                          #hhD.Fill(i1, i2, val);



"""
# validate
for k in range(len( dataRMM )):
                          x=k/mSize
                          y=k%mSize
                          i1=x;
                          i2=mSize-y-1;
                          val=float(dataRMM[k])

                          if k+1 in idx2remove:
                                 val=0
                          #hhD.Fill(Names2[i1],  Names1[i2],val);
                          #hhD.Fill(i1,i2, val);
"""

# write file
xfile=open("columns_with_0_10j10b5rest.txt","w")
for k in idx2remove:
       cell="V_"+str(k)
       xfile.write(cell+"\n")
xfile.close()

print("Total cells=",ntot," active=",active)
print("Made: columns_with_0_10j10b5rest.txt")
hhD.SetTitle("")
hhD.SetStats(0)
hhD.GetZaxis().SetRangeUser(ZSmin,ZSmax);
hhD.GetZaxis().SetLabelOffset(0.02)
hhD.GetYaxis().SetLabelOffset(0.02)
hhD.GetXaxis().SetLabelOffset(0.02)
gStyle.SetPaintTextFormat(".0e");

hhD.Draw("colz")
#hhD.Draw("text0 same")

# hmass.Draw();

x1=0
x2=mSize 
YRMIN=mSize
YRMAX=0

ar7=TLine(x1,YRMIN,x2,YRMAX);
ar7.SetLineWidth(2)
ar7.SetLineStyle(2)
ar7.SetLineColor(0)
ar7.Draw("same")


print (epsfig) 
gPad.RedrawAxis()
c1.Update()
ps1.Close()
if (myinput != "-b"):
              if (input("Press any key to exit") != "-9999"):
                         c1.Close(); sys.exit(1);




