#Used to interface between MATLAB and the Python script that determines image charges
#Performs appropriate type conversions and formatting

import dielectricInterfacesGeneralized as dig
import sys

tempInt = 1
tempDouble = 1.0
if(type(epsilon) == type(tempInt) or type(epsilon) == type(tempDouble)):
    epsilon = [epsilon]
else:
    epsilon = list(epsilon)
if(type(boundaries) == type(tempInt) or type(boundaries) == type(tempDouble)):
    boundaries = [boundaries]
else:
    boundaries = list(boundaries)
if(type(regionIndices) == type(tempInt) or type(regionIndices) == type(tempDouble)):
    regionIndices = [int(regionIndices)]
else:
    regionIndices = list(regionIndices)
    regionIndices = [int(i) for i in regionIndices]
if(type(qinit) == type(tempInt) or type(qinit) == type(tempDouble)):
    qinit = [qinit]
else:
    qinit = list(qinit)
if(type(zinit) == type(tempInt) or type(zinit) == type(tempDouble)):
    zinit = [zinit]
else:
    zinit = list(zinit)

iterations=int(iterations)

d = dig.buildDomain(epsilon,boundaries,qinit,zinit,iterations)
q = list()
x = list()
y = list()
z = list()
netCharge = dig.getNetCharge(d,regionIndices);

for i in range(len(d.regions)):
    qt = list()
    xt = list()
    yt = list()
    zt = list()
    for p in d.regions[i].pointCharges:

        #Leading ion is removed to eliminate self-interaction
        if(isLeading):
            atInitPosition = False;
            for j in range(len(zinit)):
                if(p.z == zinit[j]):
                    qt.append(p.q - qinit[j])
                    atInitPosition = True
            if(not atInitPosition):
                qt.append(p.q)
        else:
            qt.append(p.q)

        xt.append(p.x)
        yt.append(p.y)
        zt.append(p.z)
    q.append(qt)
    x.append(xt)
    y.append(yt)
    z.append(zt)
