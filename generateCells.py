def generateCells(NX,NY,noisx,noisy,direc=1,rounding=1): 
    import numpy as np
    from geomProperties import PolyArea
    from scipy.spatial import Delaunay
    from scipy.spatial import Voronoi

    print('Generating cells')

    if direc == 1:
        nx = NX
        ny = NY
        randx = noisx
        randy = noisy
        Lx = nx*2/np.sqrt(3)
        Ly = ny
    else:
        nx = NY
        ny = NX
        randx = noisy
        randy = noisx
        Lx = nx
        Ly = ny*2/np.sqrt(3)
        
    x = np.linspace(1,nx,nx)
    y = np.linspace(1,ny,ny)
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
    xv = np.array(xv)
    yv = np.array(yv)
    
    # scale the length in the direction of shift
    xv = xv*2/np.sqrt(3)

    # shifting of the alternating rows of the cells 
    xv[:,::2] += 1/np.sqrt(3)
    
    xv = xv.flatten()
    yv = yv.flatten()
    xv = xv - np.mean(xv)
    yv = yv - np.mean(yv)
    # adding random noise to the cell centers
    noisx = randx*np.random.random_sample(len(xv)) - randx/2
    noisy = randy*np.random.random_sample(len(yv)) - randy/2
    xv = xv + noisx
    yv = yv + noisy

    if direc != 1:
        tmp = xv
        xv = yv
        yv = tmp

    if rounding ==1:
        print( 'Rounding the corners')
        # round the corner of tissue
        todrop = np.where((xv/(NX/2.1))**4 + (yv/(NY/2.1))**4 > 1)
        xv = np.delete(xv,todrop)
        yv = np.delete(yv,todrop)
    cellcount = len(xv)

    (ringx,ringy) = getOuterRing(xv,yv)

    xv = np.concatenate((xv , ringx))
    yv = np.concatenate((yv , ringy))
    
    print( 'Tesellating now')
    # voronoi tesellation
    cellcenters = zip(xv,yv)
    vor = Voronoi(cellcenters)

    v = vor.vertices
    
    # centers of the cells
    xc = xv[:cellcount]
    yc = yv[:cellcount]
    
    print( 'Generating cell edges')
    cells = []
    todrop  = []
    for idx, point in enumerate(vor.points[:cellcount]):
        vn = vor.regions[vor.point_region[idx]]
        x = v[vn,0]
        y = v[vn,1]
        ar = PolyArea(x,y)
        if ar < 0: # order is clock-wise in this case
            vn.reverse()

        cells.append(vn)
        if -1 in vn:
            todrop.append(idx)

    # deleting the cells with infinity as a vertex
    for index in sorted(todrop, reverse=True):
        cells.pop(index)
    
    print( 'Calculating cell neighbors')
    nbs = [[] for i in range(len(cells))]
    cellnum1 = 0
    for cell1 in cells:
        for i in range(len(cell1)):
            if i < len(cell1)-1:
                edge1 = [cell1[i+1],cell1[i]]
            else:
                edge1 = [cell1[0],cell1[i]]
            cellnum2 = 0
            nbedge = -1
            for cell2 in cells:
                for j in range(len(cell2)):
                    if j < len(cell2)-1:
                        edge2 = [cell2[j],cell2[j+1]]
                    else:
                        edge2 = [cell2[j],cell2[0]]
                    if edge1 == edge2:
                        nbedge = cellnum2
                cellnum2 += 1
            nbs[cellnum1].append(nbedge)
        cellnum1 += 1

    return v,cells,nbs

def generateCellsPeriodic(NX,NY,noisx,noisy,direc=1,scalingx=1,scalingy=1): 
    import numpy as np
    from geomProperties import PolyArea
    from scipy.spatial import Delaunay
    from scipy.spatial import Voronoi

    print('Generating cells')
    nx = NX; ny = NY; randx = noisx; randy = noisy;
#    Lx = scalingx*nx*2/np.sqrt(3); Ly = scalingy*ny;
    Lx = nx*2/np.sqrt(3); Ly = ny;
        
#    x = np.linspace(1,scalingx*nx,nx)
#    y = np.linspace(1,scalingy*ny,ny)
    x = np.linspace(1,nx,nx)
    y = np.linspace(1,ny,ny)
    xv, yv = np.meshgrid(x, y, sparse=False, indexing='ij')
    xv = np.array(xv)
    yv = np.array(yv)
    
    # scale the length in the direction of shift
    xv = xv*2/np.sqrt(3)

    # shifting of the alternating rows of the cells 
    xv[:,::2] += 1/np.sqrt(3)
    
    xv = xv.flatten()
    yv = yv.flatten()
    xv = xv - np.mean(xv)
    yv = yv - np.mean(yv)
    # adding random noise to the cell centers
    noisx = randx*np.random.random_sample(len(xv)) - randx/2
    noisy = randy*np.random.random_sample(len(yv)) - randy/2
    xv = xv + noisx
    yv = yv + noisy

    # Appending the periodic 
    xv = np.concatenate((xv,xv+Lx)) # EAST
    yv = np.concatenate((yv,yv   )) # EAST
    xv = np.concatenate((xv,xv-Lx)) # WEST
    yv = np.concatenate((yv,yv   )) # WEST
    xv = np.concatenate((xv,xv   )) # NORTH
    yv = np.concatenate((yv,yv+Ly)) # NORTH
    xv = np.concatenate((xv,xv   )) # SOUTH
    yv = np.concatenate((yv,yv-Ly)) # SOUTH
    xv = np.concatenate((xv,xv+Lx)) # NORTH-EAST
    yv = np.concatenate((yv,yv+Ly)) # NORTH-EAST
    xv = np.concatenate((xv,xv+Lx)) # NORTH-WEST
    yv = np.concatenate((yv,yv-Ly)) # NORTH-WEST
    xv = np.concatenate((xv,xv+Lx)) # SOUTH-EAST
    yv = np.concatenate((yv,yv-Ly)) # SOUTH-EAST
    xv = np.concatenate((xv,xv-Lx)) # SOUTH-WEST
    yv = np.concatenate((yv,yv-Ly)) # SOUTH-WEST

    cellcount = nx*ny

    print('Tesellating now')
    # voronoi tesellation
    cellcenters = zip(xv,yv)
    cellcenters = np.vstack((xv.T,yv.T)).T
    vor = Voronoi(cellcenters)

    v = vor.vertices

    v[:,0] *= scalingx; 
    v[:,1] *= scalingy
    Lx *= scalingx
    Ly *= scalingy

    # centers of the cells
    xc = xv[:cellcount]
    yc = yv[:cellcount]
    
    print('Generating cell edges')
    cells = []
    todrop  = []
    #for idx, point in enumerate(vor.points[:cellcount]):
    for idx, point in enumerate(vor.points):
        vn = vor.regions[vor.point_region[idx]]
        x = v[vn,0]
        y = v[vn,1]
        ar = PolyArea(x,y)
        if ar < 0: # order is clock-wise in this case
            vn.reverse()

        cells.append(vn)
        if -1 in vn:
            todrop.append(idx)

    # deleting the cells with infinity as a vertex
    #for index in sorted(todrop, reverse=True):
    #    cells.pop(index)
    
    print('Calculating cell neighbors')
    nbs = [[] for i in range(len(cells[:cellcount]))]
    cellnum1 = 0
    for cell1 in cells[:cellcount]:
        for i in range(len(cell1)):
            if i < len(cell1)-1:
                edge1 = [cell1[i+1],cell1[i]]
            else:
                edge1 = [cell1[0],cell1[i]]
            cellnum2 = 0
            nbedge = -1
            for cell2 in cells:
                for j in range(len(cell2)):
                    if j < len(cell2)-1:
                        edge2 = [cell2[j],cell2[j+1]]
                    else:
                        edge2 = [cell2[j],cell2[0]]
                    if edge1 == edge2:
                        nbedge = cellnum2
                        break
                if nbedge == -1:
                    cellnum2 += 1
                else:
                    break
            cellno = nbedge % (cellcount)
            nbs[cellnum1].append(cellno)
        cellnum1 += 1
    cells = cells[:cellcount]

#    # identifying the twin vertices
#    brdrpoints = np.array([])
#    for cl in range(len(cells)):
#        nb = nbs[cl]
#        vr = np.array(cells[cl])
#        brdr = np.where(np.array(nb)==-1)[0]
#        vbrdr = np.concatenate((vr[brdr], np.roll(vr,-1)[brdr]))
#        vbrdr = np.unique(vbrdr)
#        
#        brdrpoints = np.concatenate((brdrpoints,vbrdr))
#    brdrpoints = np.unique(brdrpoints)
#    brdrpoints = brdrpoints.astype(int)
#    twins_ew = []
#    twins_ns = []
#    for vbrdr in brdrpoints:
#        d = np.sqrt((v[brdrpoints,0]-Lx-v[vbrdr,0])**2 + (v[brdrpoints,1]-v[vbrdr,1])**2)
#        loc = np.where(d<1e-12)[0]
#        if loc.size == 1:
#            if v[int(vbrdr),0] < v[brdrpoints[loc[0]],0]:
#                twins_ew.append([int(vbrdr),brdrpoints[loc[0]]])
#            else:
#                twins_ew.append([brdrpoints[loc[0]],int(vbrdr)])
#        d = np.sqrt((v[brdrpoints,0]+Lx-v[vbrdr,0])**2 + (v[brdrpoints,1]-v[vbrdr,1])**2)
#        loc = np.where(d<1e-12)[0]
#        if loc.size == 1:
#            if v[int(vbrdr),0] < v[brdrpoints[loc[0]],0]:
#                twins_ew.append([int(vbrdr),brdrpoints[loc[0]]])
#            else:
#                twins_ew.append([brdrpoints[loc[0]],int(vbrdr)])
#
#        d = np.sqrt((v[brdrpoints,0]-v[vbrdr,0])**2 + (v[brdrpoints,1]-Ly-v[vbrdr,1])**2)
#        loc = np.where(d<1e-12)[0]
#        if loc.size  == 1:
#            if v[int(vbrdr),1] < v[brdrpoints[loc[0]],1]:
#                twins_ns.append([int(vbrdr),brdrpoints[loc[0]]])
#            else:
#                twins_ns.append([brdrpoints[loc[0]],int(vbrdr)])
#        d = np.sqrt((v[brdrpoints,0]-v[vbrdr,0])**2 + (v[brdrpoints,1]+Ly-v[vbrdr,1])**2)
#        loc = np.where(d<1e-12)[0]
#        if loc.size  == 1:
#            if v[int(vbrdr),1] < v[brdrpoints[loc[0]],1]:
#                twins_ns.append([int(vbrdr),brdrpoints[loc[0]]])
#            else:
#                twins_ns.append([brdrpoints[loc[0]],int(vbrdr)])
#    
#    print 'brdrpoints, twins_ew, twins_ns, ',len(brdrpoints), len(twins_ew),len(twins_ns)
                

    if direc != 1:
        vx = 1.0*v[:,0]
        vy = 1.0*v[:,1]
        v[:,0] = 1.0*vy
        v[:,1] = 1.0*vx
        temp0 = 1.0*Lx
        Lx = 1.0*Ly
        Ly = 1.0*temp0

    return v,cells,nbs,Lx,Ly


    
def getOuterRing(xv,yv):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.spatial import ConvexHull

    cellcenters = zip(xv,yv)
    hull = ConvexHull(cellcenters)
    hull = hull.vertices # counterclock-wise hull
    
    x = xv[hull]
    y = yv[hull]
    x = np.append(x,x[0])
    y = np.append(y,y[0])
    
    (hullx,hully) = refineCurve(x,y)

    ringx = []
    ringy = []
    for i in range(len(hullx)-1):
        v1x = hullx[i]
        v1y = hully[i]
        v2x = hullx[i+1]
        v2y = hully[i+1]
        midx = 0.5*(v1x+v2x)
        midy = 0.5*(v1y+v2y)
        rx = v2x - v1x
        ry = v2y - v1y
        r = np.sqrt(rx*rx+ry*ry)
        rx = rx/r
        ry = ry/r
        vxnew = midx + ry
        vynew = midy - rx
        
        ringx.append(vxnew)
        ringy.append(vynew)
    
    # close the ring
    ringx.append(ringx[0])
    ringy.append(ringy[0])

    (ringx,ringy) = refineCurve(ringx,ringy)
    return (ringx,ringy)

def refineCurve(x,y):
    import numpy as np
    hullx = []
    hully = []
    for i in range(len(x)-1):
        rx = x[i+1]-x[i]
        ry = y[i+1]-y[i]
        d = np.sqrt(rx*rx+ry*ry)
        d = np.ceil(d)
        t = np.linspace(0,1,d+1)
        xx = (1-t)*x[i] + t*x[i+1]
        yy = (1-t)*y[i] + t*y[i+1]
        
        hullx += xx[:-1].tolist()
        hully += yy[:-1].tolist()  
    # close the hull  
    hullx += [xx[-1]]
    hully += [yy[-1]]
    return (hullx,hully)

if __name__ == "__main__":
#    import times
    import numpy as np
    import matplotlib.pyplot as plt
    
    from geomProperties import PolyArea
#    t = time.time() # codetesting

    #v,cells,nbs = generateCells(4,4,0.6,0.6,direc=1)
#    v,cells,nbs = generateCells(20,20,0.8,0.8,direc=1,rounding=1)
    v,cells,nbs,Lx,Ly = generateCellsPeriodic(4,10,0.3,0.3,direc=1)
#    print Lx,Ly
#    print( time.time() - t)
    ar = []
    for i in range(len(cells[:100])):
        xv = v[cells[i],0]
        xv = np.append(xv,xv[0]) # to close the polygon
        yv = v[cells[i],1]
        yv = np.append(yv,yv[0]) # to close the polygon
        plt.plot(xv,yv,'k'); plt.plot(xv+Lx,yv); plt.plot(xv,yv+Ly); plt.plot(xv-Lx,yv); plt.plot(xv,yv-Ly)
        plt.plot(xv-Lx,yv-Ly,'k'); plt.plot(xv-Lx,yv+Ly); plt.plot(xv+Lx,yv-Ly); plt.plot(xv+Lx,yv+Ly);
        ar.append(PolyArea(xv[:-1],yv[:-1]))
    
    cn = 9
    xv = v[cells[cn],0]
    yv = v[cells[cn],1]
    plt.plot(np.mean(xv),np.mean(yv),'go')
    nb = nbs[cn]
#    print nb
    for n in nb:
        xv = v[cells[n],0]
        yv = v[cells[n],1]
        plt.plot(np.mean(xv),np.mean(yv),'ro')
    plt.axis('equal')
#    ax.set_xlim([0,11])
    plt.show()

    print( np.std(ar))
