#!/usr/bin/env python
'''
read, process absorbance data provided as csv file
20151203
'''
import math
import csv
import graphUtils,photons,ProcessMaterialProperty
from ROOT import TMultiGraph,TFile

fn = '/Users/djaffe/Documents/Neutrinos/LDRD2010/OneTonPrototypeIn2-224/Water/RO_Test_11-24.csv'
data = {}
wavelength = []
headers = []
ymin = 0.
with open(fn,'r') as csvfile:
    for l in csv.reader(csvfile):
        if len(headers)==0:
            headers = l
            for h in headers: data[h] = []
        else:
            for h,x in zip(headers,l):
                data[h].append(float(x))
                ymin = min(ymin,float(x))
                if h=='Wavelength':
                    wavelength.append(float(x))
                    

ymin = 1.001*ymin # no zeros
                    
gU = graphUtils.graphUtils()
ph = photons.photons()
pmp = ProcessMaterialProperty.ProcessMaterialProperty()

tmg = TMultiGraph()
tmg.SetName('tmg')
tmg.SetTitle('#splitline{Absorbance in 10cm cell}{Adjusted to be non-zero}')
amg = TMultiGraph()
amg.SetName('amg')
amg.SetTitle('Attenuation length (cm)')
bmg = TMultiGraph()
bmg.SetName('bmg')
bmg.SetTitle('#splitline{Effective effy for}{100cm path length}')
cmg = TMultiGraph()
cmg.SetName('cmg')
cmg.SetTitle('Attenuation in 100cm')
dmg = TMultiGraph()
dmg.SetName('dmg')
dmg.SetTitle('#splitline{Absorbance in 10cm cell.}{Adj to 0 at 600nm}') # SetTitle("#splitline{aaa}{bbb}")

MGlist = [tmg,amg,bmg,cmg,dmg]

Graphs = {}

X = wavelength
iend = wavelength.index(801.)  # only interested in 200-800nm
title = 'QE_R7723'
QE = []
for wl in X[:iend]:
    QE.append( ph.getQE(wl,material=title) )
    
name  = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X[:iend],QE,title,name)
gU.drawGraph(g)

pathLength = 100. # cm
pathLmm    = pathLength*10
pureWater = []
for wl in X[:iend]:
    pureWater.append(ph.getAtten(wl, pathLmm, medium='water'))
title = 'Pure water attenuation in'+str(int(pathLength))+'cm'
name  = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X[:iend],pureWater,title,name)
gU.drawGraph(g)
cmg.Add(g)

waterAbsLen = []
for wl in X[:iend]: # convert from mm to cm
    waterAbsLen.append( 0.1*pmp.getMaterialProperty(wl, material='water', property='ABSLENGTH') )
title = 'pure water abslength in cm'
name  = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X[:iend],waterAbsLen,title,name)
gU.drawGraph(g)
amg.Add(g)

E = [eff*math.exp(-pathLength/attl) for eff,attl in zip(QE,waterAbsLen[:iend])]
title = 'Eff Attl pure water'
name = title.replace(' ','_').replace('.','')
Graphs[name] = g = gU.makeTGraph(X[:iend],E[:iend],title,name)
#gU.color(g,icol,icol,setMarkerColor=True)
bmg.Add(g)


waterAbsorbance = []
wAat600 = None
for al,wl in zip(waterAbsLen,X[:iend]):
    q =  -math.log10( math.exp(-10./al) ) 
    waterAbsorbance.append( q )
    if wl==600. : wAat600 = q
title = 'pure water absorbance'
name  = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X[:iend],waterAbsorbance,title,name)
tmg.Add(g)

print 'wAat600',wAat600

Y = [a-wAat600 for a in waterAbsorbance]
title = 'pure water absorbance adjusted to zero at 600 nm'
name = title.replace(' ','_')
Graphs[name] = g = gU.makeTGraph(X[:iend],Y,title,name)
dmg.Add(g)

cellLength = 10. # cm
ln10 = math.log(10.)
print 'ymin',ymin
icol = 0
for h in headers:
    if h!='Wavelength':
        icol += 1
        Y = [a-ymin for a in data[h]]
        print h,'min(Y)',min(Y)
        title = h
        name = title.replace(' ','_').replace('.','')
        Graphs[name] = g = gU.makeTGraph(X[:iend],Y[:iend],title,name)
        gU.color(g,icol,icol,setMarkerColor=True)
        tmg.Add(g)

        title = 'raw ' + h
        name = title.replace(' ','_').replace('.','')
        Graphs[name] = g = gU.makeTGraph(X[:iend],data[h][:iend],title,name)
        gU.color(g,icol,icol,setMarkerColor=True,setMarkerType=False)
        dmg.Add(g)

        A = [cellLength/a/ln10 for a in Y]
        title = 'Attenuation length ' + h
        name = title.replace(' ','_').replace('.','')
        Graphs[name] = g = gU.makeTGraph(X[:iend],A[:iend],title,name)
        gU.color(g,icol,icol,setMarkerColor=True)
        amg.Add(g)

        C = [math.exp(-pathLength/attl) for attl in A[:iend]]
        title = 'Attenuation in '+str(int(pathLength))+'cm ' + h
        name = title.replace(' ','_').replace('.','')
        Graphs[name] = g = gU.makeTGraph(X[:iend],C[:iend],title,name)
        gU.color(g,icol,icol,setMarkerColor=True)
        cmg.Add(g)


        E = [eff*math.exp(-pathLength/attl) for eff,attl in zip(QE,A[:iend])]
        title = 'Eff Attl ' + h
        name = title.replace(' ','_').replace('.','')
        Graphs[name] = g = gU.makeTGraph(X[:iend],E[:iend],title,name)
        gU.color(g,icol,icol,setMarkerColor=True)
        bmg.Add(g)


        

for mg in MGlist:
    gU.drawMultiGraph(mg,abscissaIsTime=False)
    gU.drawMultiGraph(mg,SetLogy=True,abscissaIsTime=False)

outfilename = 'water.root'
ofn = TFile(outfilename,"RECREATE")
for g in Graphs:
    ofn.WriteTObject(Graphs[g])
for mg in MGlist:
    ofn.WriteTObject(mg)
ofn.Close()
print 'wrote',len(Graphs),'graphs to',outfilename
    
