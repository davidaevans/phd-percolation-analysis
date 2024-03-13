import numpy as np
import sys
import os
from math import sqrt

datafilepaths = []
data = [] 

#print("Datafiles:\n")
for arg in sys.argv:
    #skip name of program
    if arg == sys.argv[0]: continue
    #add paths to list to store them
    datafilepaths.append(arg)
    #open data and add contents to np array in data list
    data.append(np.loadtxt(arg, delimiter=" "))
    #print(arg)

rhos = []
for datum in data:
    for i in range(np.size(datum,0)):
        if datum[i,0] not in rhos:
            rhos.append(datum[i,0])
rhos.sort()
#print(rhos)

# array is arr[rho][x,y,xy,xandy,cs,cs_expc,coord][run_i]
# set elements to -1 initially so that 
organisedData = np.empty((len(rhos), 7, len(data)))

#loop over each data file
for j in range(len(data)):
    datum = data[j]
    #loop over each density
    for i in range(len(rhos)):
        #loop through rows in file
        for k in range(np.size(datum,0)):
            if datum[k][0] == rhos[i]:
                try:
                    organisedData[i][0][j] = datum[i][2]
                    organisedData[i][1][j] = datum[i][4]
                    organisedData[i][2][j] = datum[i][6]
                    organisedData[i][3][j] = datum[i][8]
                    organisedData[i][4][j] = datum[i][10]
                    organisedData[i][5][j] = datum[i][11]
                    organisedData[i][6][j] = datum[i][12]
                except IndexError:
                    continue

print("#density x_mean x_sd x_sem y_mean y_sd y_sem xy_mean xy_sd xy_sem xandy_mean xandy_sd xand_sem cs_mean cs_sd cs_sem cs_expc_mean cs_expc_sd cs_expc_sem coordnum_ave coordnum_sd coordnum_sem ")
#start averaging for each density.
for i in range(len(rhos)):
    x=[]
    y=[]
    xy =[]
    xandy = []
    cs=[]
    cs_expc = []
    coord = []
    xsum = ysum = xysum = xandysum = cssum = cs_expcsum = coordsum = 0
    xcount = ycount = xycount = xandycount = cscount = cs_expccount = coordcount = 0

    for j in range(len(data)):
        if organisedData[i][0][j] != -1: 
            x.append(organisedData[i][0][j])
        if organisedData[i][1][j] != -1: 
            y.append(organisedData[i][1][j])
        if organisedData[i][2][j] != -1: 
            xy.append(organisedData[i][2][j])
        if organisedData[i][3][j] != -1:
            xandy.append(organisedData[i][3][j])
        if organisedData[i][4][j] != -1: 
            cs.append(organisedData[i][4][j])
        if organisedData[i][5][j] != -1:
            cs_expc.append(organisedData[i][5][j])
        if organisedData[i][5][j] != -1:
            coord.append(organisedData[i][6][j])

    xave = sum(x)/len(x)
    yave = sum(y)/len(y)
    xyave = sum(xy)/len(xy)
    xandyave = sum(xandy)/len(xandy)
    csave = sum(cs)/len(cs)
    cs_expcave = sum(cs_expc)/len(cs_expc)
    coordave = sum(coord)/len(coord)

    xsd = [(xi-xave)**2 for xi in x]
    xstddev = sqrt(sum(xsd)/(len(x)-1))
    xseom = xstddev/sqrt(len(x))

    ysd = [(yi-yave)**2 for yi in y]
    ystddev = sqrt(sum(ysd)/(len(y)-1))
    yseom = ystddev/sqrt(len(y))

    xysd = [(xyi-xyave)**2 for xyi in xy]
    xystddev = sqrt(sum(xysd)/(len(xy)-1))
    xyseom = xystddev/sqrt(len(xy))

    xandysd = [(xandyi-xandyave)**2 for xandyi in xandy]
    xandystddev = sqrt(sum(xandy)/(len(xandy)-1))
    xandyseom = xandystddev/sqrt(len(xandy))

    cssd = [(csi-csave)**2 for csi in cs]
    csstddev = sqrt(sum(cssd)/(len(cs)-1))
    csseom = csstddev/sqrt(len(cs))

    cs_expcsd = [(cs_expci-cs_expcave)**2 for cs_expci in cs_expc]
    cs_expcstddev = sqrt(sum(cs_expc)/(len(cs_expc)-1))
    cs_expcseom = cs_expcstddev/sqrt(len(cs_expc))

    coordsd = [(coordi-coordave)**2 for coordi in coord]
    coordstddev = sqrt(sum(coord)/(len(coord)-1))
    coordseom = coordstddev/sqrt(len(coord))

    print("{:.3f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.3f} {:.3f} {:.3f}".format(rhos[i], xave, xstddev, xseom, yave, ystddev, yseom, xyave, xystddev, xyseom, xandyave, xandystddev, xandyseom, csave, csstddev, csseom, cs_expcave, cs_expcstddev, cs_expcseom, coordave, coordstddev, coordseom))
