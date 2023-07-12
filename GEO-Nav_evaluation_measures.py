import numpy as np
import math
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

model = "5EK0"

number="1"
# set number="1" or ="2" depending on the output of HOLE to be analyzed

threshold=1
sC=200
sH=200


print("\n--- Model: {} ---".format(model))

def FromCloudToOFF(vertices, offname):

    nv = vertices.shape[0]

    f = open('example/output/{}.off'.format(offname), "w")
    f.write("OFF\n\n")
    f.write("{} 0 0 \n".format(nv))

    for i in range(nv):
        f.write("{} {} {} \n".format(
            vertices[i][0], vertices[i][1], vertices[i][2]))

    f.close()

    return 0


def createColoredOFFfunction(offname, vertices, Fv):

    nv=len(Fv)
    min_val = min(Fv)

    max_val=0

    for d in Fv:
         if ( d>max_val and d != np.inf):
            max_val=d

    check=0

    for d in Fv:
        if d == np.inf:
            d = max_val*2
            check=1

    if check == 1:
        max_val = max_val

    # use the coolwarm colormap that is built-in, and goes from blue to red
    cmap = mpl.cm.coolwarm
    normalized = mpl.colors.SymLogNorm(linthresh=1, vmin=min_val, vmax=max_val)

    # convert your distances to color coordinates
    Cv = cmap(normalized(Fv))

    f = open('example/output/{}_color.off'.format(offname), "w")
    f.write("COFF\n\n")
    f.write("{} 0 0 \n".format(nv))

    for i in range(nv):
        f.write("{} {} {} {} {} {} {}\n".format(vertices[i][0], vertices[i][1], vertices[i][2], Cv[i,0],Cv[i,1],Cv[i,2],Cv[i,3]))

    f.close()

    return 0


def dist(u, v):

    return math.sqrt((u[0]-v[0])**2+(u[1]-v[1])**2+(u[2]-v[2])**2)

def closest_dist(u, U):
    Cdist=np.infty
    Idist=np.infty

    for i in range(np.shape(U)[0]):
        if dist(U[i], u) < Cdist:
            Cdist=dist(U[i], u)
            Idist=i

    return Cdist, Idist


def Length(data):

    l=0

    for i in range(np.shape(data)[0]-1):
        l=l+dist(data[i,:],data[i+1,:])

    return l

def distance(u, v, w):

    l = w[0]-v[0]
    m = w[1]-v[1]
    n = w[2]-v[2]
    t = (l*(u[0]-v[0])+m*(u[1]
         - v[1])+n*(u[2]-v[2]))/(l**2+m**2+n**2)
    xH = v[0]+t*l
    yH = v[1]+t*m
    zH = v[2]+t*n

    return dist(u, [xH,yH,zH])

def Turt(data):

    t = 0

    for i in range(1,np.shape(data)[0]-1):
        t = t+distance(data[0,:], data[i,:], data[np.shape(data)[0]-1,:])
        i = i+1

    t = t/np.shape(data)[0]

    return t

def Volume(data):

    vol=0

    hMax=0

    for i in range(data.shape[0]-1):

        h=dist(data[i], data[i+1])

        if h > hMax:
            hMax = h

        cone = (data[i][3]**2 + data[i][3]*data[i+1][3] + data[i+1][3]**2) * h

        vol = vol + cone

    vol = 1/3 * math.pi * vol

    return vol

### LOAD CHANALYZER ###

infileC = "example/input/centerline_Chanalyzer.csv".format(model)
# infileC = "results/dataset/Chanalyzer/{}/centerline_Chanalyzer.csv".format(model)

chan = np.loadtxt(infileC)

nC = np.shape(chan)

print("Number of Vertices Chanalyzer: {}".format(nC[0]))

### LOAD HOLE ###

hole=[]
hole=set()

infileH = open("example/input/hole_out_{}.sph.pdb".format(number), 'r')
# infileH = open("results/dataset/hole_test/{}/hole_out_{}.sph.pdb".format(model,number), 'r')

Lines = infileH.readlines()

for line in Lines:

    if line.split()[0] != "LAST-REC-END":

        if float(line.split()[-2]) == float(line.split()[-1]):

            hole.add((line.split()[-5],line.split()[-4],line.split()[-3],line.split()[-2]))

hole=list(hole)

infileH.close()

hole = np.asarray(hole, dtype=float)

nH = np.shape(hole)

print("Number of Vertices Hole {}: {}".format(number, nH[0]))

distMat=np.zeros((nH[0],nH[0]))

GH=nx.Graph()
for i in range(nH[0]):
    GH.add_node(i)
# for i in range(8):
    for j in range(i):
        GH.add_edge(i, j, weight=dist(hole[i],hole[j]) )
        GH.add_edge(j, i, weight=dist(hole[i],hole[j]) )
        distMat[i,j]=dist(hole[i],hole[j])

MaxEnt=np.unravel_index(np.argmax(distMat, axis=None), distMat.shape)

GH.add_edge(MaxEnt[0], MaxEnt[1], weight=0)

TSPHprimo=nx.approximation.traveling_salesman_problem(GH)

TSPH=[]

for i in range(len(TSPHprimo)):
    if TSPHprimo[i] not in TSPH:
        TSPH.append(TSPHprimo[i])

listH=[]

ind=TSPH.index(MaxEnt[0])

if TSPH.index(MaxEnt[0])>TSPH.index(MaxEnt[1]):

    for i in range(len(TSPH)):

        j=(ind+i)%(len(TSPH))

        listH.append([hole[TSPH[j],0],hole[TSPH[j],1],hole[TSPH[j],2],hole[TSPH[j],3]])

else:

    for i in range(len(TSPH)):

        j=(ind-i)%(len(TSPH))

        listH.append([hole[TSPH[j],0],hole[TSPH[j],1],hole[TSPH[j],2],hole[TSPH[j],3]])

hole = np.asarray(listH, dtype=float)

### OFF FILES ###

createColoredOFFfunction(model+"_chanalyzer", chan, chan[:,3])

createColoredOFFfunction(model+"_hole_"+number, hole, hole[:,3])

### OFF FILES (SPHERES) ###

chan_spheres = []
hole_spheres = []


for j in range(nC[0]):

    r=chan[j,3]

    for i in range(sC):
        A = random.uniform(0, 2*math.pi)
        z = random.uniform(-r,+r)

        v=[ math.sqrt(r**2-z**2)*math.cos(A)+chan[j,0], math.sqrt(r**2-z**2)*math.sin(A)+chan[j,1], z+chan[j,2], r ]

        toBeAddedC=True

        for k in range(nC[0]):
            if k!=j and dist(chan[k],v) <= r:
                toBeAddedC=False
                break
        if toBeAddedC==True:
            chan_spheres.append( v )

chan_spheres = np.asarray(chan_spheres, dtype=float)

for j in range(nH[0]):

    r=hole[j,3]

    for i in range(sH):
        A = random.uniform(0, 2*math.pi)
        z = random.uniform(-r,+r)

        v=[ math.sqrt(r**2-z**2)*math.cos(A)+hole[j,0], math.sqrt(r**2-z**2)*math.sin(A)+hole[j,1], z+hole[j,2], r ]

        toBeAdded=True

        for k in range(nH[0]):
            if k!=j and dist(hole[k],v) <= r:
                toBeAdded=False
                break
        if toBeAdded==True:
            hole_spheres.append( v )

hole_spheres = np.asarray(hole_spheres, dtype=float)

createColoredOFFfunction(model+"_spheres_chanalyzer", chan_spheres, chan_spheres[:,3])

createColoredOFFfunction(model+"_spheres_hole"+number, hole_spheres, hole_spheres[:,3])

### ANALYSIS ###

CXmap = [0]
CYmap = [chan[0,3]]

for i in range(1,nC[0]):

    CXmap.append(CXmap[i-1]+dist(chan[i-1,:],chan[i,:]))
    CYmap.append(chan[i,3])


hFirst=np.infty
disthFirst=np.infty

for i in range(nH[0]):

    if dist(hole[i,:], chan[0,:]) < disthFirst:
        hFirst=i
        disthFirst=dist(hole[i], chan[0])

if dist(hole[(hFirst+5)%(nH[0]-1),:], chan[5,:]) <= dist(hole[(hFirst-5)%(nH[0]-1),:], chan[5,:]):
    dir=1
else:
    dir=-1

#######

HXmap = []
HYmap = []
HXind = []

distPlus=0
distMinus=0

for i in range(hFirst):
    distPlus=distPlus+dist(hole[i,:],hole[i+1,:])

for i in range(hFirst, nH[0]-1):
    distMinus=distMinus+dist(hole[i,:],hole[i+1,:])

if dir>0:
    HXmap.append(-distPlus)
    HYmap.append(hole[0,3])
    HXind.append(0)
    for i in range(1,nH[0]):
        HXmap.append(HXmap[i-1]+dist(hole[i,:],hole[i-1,:]))
        HYmap.append(hole[i,3])
        HXind.append(i)

if dir<0:
    HXmap.append(-distMinus)
    HYmap.append(hole[nH[0]-1,3])
    k=0
    for i in range(nH[0]-2,-1,-1):
        HXmap.append(HXmap[k]+dist(hole[i,:],hole[i+1,:]))
        HYmap.append(hole[i,3])
        HXind.append(i)
        k=k+1

DHXmap=[]
DHYmap=[]

NumAlignC=0
NumAlignH=0
XAlignC=[]
XAlignH=[]

dRho=0

for i in range(nC[0]):
    DHXmap.append(CXmap[i])
    DHYmap.append(closest_dist(chan[i], hole)[0])
    if closest_dist(chan[i], hole)[0] < threshold:
        NumAlignC=NumAlignC+1
        XAlignC.append(CXmap[i])
        dRho=dRho+(chan[i,3]-hole[closest_dist(chan[i], hole)[1],3])**2

i=0
for j in HXind:
    DHXmap.append(HXmap[i])
    DHYmap.append(closest_dist(hole[j], chan)[0])
    if closest_dist(hole[j], chan)[0] < threshold:
        NumAlignH=NumAlignH+1
        XAlignH.append(HXmap[i])
    i=i+1

OneMin=max(min(XAlignC),min(XAlignH))
OneMax=min(max(XAlignC),max(XAlignH))


plt.scatter(CXmap, CYmap, s=3, c='tab:blue', label='Chanalyzer')

plt.scatter(HXmap, HYmap, s=3, c='tab:orange', label='Hole {}'.format(number))

if OneMin < np.infty and OneMax > -np.infty:
    plt.axvline(x=OneMin, color='y', linestyle=':')
    plt.axvline(x=OneMax, color='y', linestyle=':')

plt.title("Radius Function Model: {} Chanalyzer vs Hole {}".format(model, number))
plt.legend()
plt.savefig("example/output/Radius_Model{}_Hole{}.png".format(model, number))


print("\n- Length -")

print("Length Chanalyzer: {}".format(round(Length(chan),2)))
print("Length Hole {}: {}".format(number, round(Length(hole),2)))


print("\n- Straightness -")

print("Straightness Chanalyzer: {}".format(round(1/Turt(chan),2)))
print("Straightness Hole {}: {}".format(number, round(1/Turt(hole),2)))


print("\n- Alignement -")

print("Match Chanalyzer - Hole: {} %".format(round(100*NumAlignC/nC[0],2)))
print("Match Hole - Chanalyzer: {} %".format(round(100*NumAlignH/nH[0],2)))
print("Distance of Radius Functions in Aligned Portion: {}".format(round(dRho,2)))

print("\n")

plt.show()
