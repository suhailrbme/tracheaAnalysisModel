class singleCell():
    def __init__(self,xcell=None,ycell=None,Lx=None,dh=0.1,nl=5, x3d=None,y3d=None,z3d=None,conn=None,cyl=False,R=None):
        import numpy as np
        def positive(i,N):
            return i*(i>=0) + (i+N)*(i<0)
        
        if x3d is None: # generate cell from scratch with given cross-section contour
            # caps- apical and basal
            from mesh2dPoly import create_mesh
            points, cells = create_mesh(np.array(list(zip(xcell,ycell))),edgelength=0.5, max_area=0.01)

            if cyl:
                # transforming from plane to cylinder
#                R = 2*Lx/(2*np.pi)
                theta = xcell/R
                xcell = R*np.cos(theta)
                zcell = R*np.sin(theta)
                # ycell remains unchanged
                # normals to the cylinder surface
                norx = xcell/R 
                norz = zcell/R
                nory = 0*ycell
            else: # not cylinder
                norx = 0*xcell
                nory = 0*ycell
                norz = 1.0+0*ycell
                zcell = 0*xcell

            Nm = len(xcell)
            self.x3d = []
            self.y3d = []
            self.z3d = []
            self.conn = [] # connectivity matrix
            for li in range(nl):
                self.x3d += (xcell + li*dh*norx).tolist()
                self.y3d += ycell.tolist()
                self.z3d += (zcell + li*dh*norz).tolist()
                
                if li > 0: # second layer onwards
                    for nmi in range(Nm):
                        if nmi == 0:
                            self.conn.append([li*Nm, li*Nm+Nm-1,(li-1)*Nm+Nm-1])
                            self.conn.append([li*Nm, (li-1)*Nm+Nm-1, (li-1)*Nm])
                        else:
                            self.conn.append([positive(li*Nm+nmi,Nm), positive(li*Nm+nmi-1,Nm), positive((li-1)*Nm+nmi-1,Nm)])
                            self.conn.append([positive(li*Nm+nmi,Nm), positive((li-1)*Nm+nmi-1,Nm), positive((li-1)*Nm+nmi,Nm)])

            # apical cap
            topx = points[len(xcell):,0]
            topy = points[len(xcell):,1]
            if cyl:
                theta = topx/R
                topx = (R+(nl-1)*dh)*np.cos(theta)
                topz = (R+(nl-1)*dh)*np.sin(theta)
            else:
                topz = (nl-1)*dh + 0*topx

            cells_a = cells + len(self.x3d) - len(xcell)
            self.x3d += topx.tolist()
            self.y3d += topy.tolist()
            self.z3d += topz.tolist()
            self.conn += cells_a.tolist()
            
            # basal cap
            cells_b = 1*cells
            cells_b[cells_b>=len(xcell)] += len(self.x3d) - len(xcell)
            if cyl:
                topx = R*np.cos(theta)
                topz = R*np.sin(theta)            
            else:
                topz = 0*topz
            self.x3d += topx.tolist()
            self.y3d += topy.tolist()
            self.z3d += topz.tolist()
            cells = 1*cells_b
            cells_b[:,0] = cells[:,1] # these two lines are to switch the orientation of the triangles
            cells_b[:,1] = cells[:,0]
            self.conn += cells_b.tolist()

            self.x3d = np.array(self.x3d)
            self.y3d = np.array(self.y3d)
            self.z3d = np.array(self.z3d)

        if xcell is None: # construct cell after reading data from vtk file
            self.x3d = x3d
            self.y3d = y3d
            self.z3d = z3d
            self.conn = conn

        self.R = R
    
    # geometric properties calculation
    # volume, curvature tensor, mean curvature, gaussian curvature
    def surfaceNormals_TotalArea_GaussianMeanCurvatures_Volume(self):
        import numpy as np
        Arcell = 0.0
        Vcell = 0.0
        ar = 0*self.x3d
        nx = 0*self.x3d
        ny = 0*self.y3d
        nz = 0*self.z3d
        K = 2*np.pi + 0*self.x3d
        Hx = 0*self.x3d
        Hy = 0*self.y3d
        Hz = 0*self.z3d
        for tri in self.conn:
            xs = self.x3d[tri]
            ys = self.y3d[tri]
            zs = self.z3d[tri]
            
#            xm = np.mean(xs)
#            ym = np.mean(ys)
#            zm = np.mean(zs)
#            rm2 = xm*xm + zm*z
#            if abs(rm2-self.R*self.R) > 1e-2:
#                continue

            dx10 = xs[1]-xs[0]; dx20 = xs[2]-xs[0]
            dy10 = ys[1]-ys[0]; dy20 = ys[2]-ys[0]
            dz10 = zs[1]-zs[0]; dz20 = zs[2]-zs[0]

            dx01 = -dx10; dx21 = xs[2]-xs[1]
            dy01 = -dy10; dy21 = ys[2]-ys[1]
            dz01 = -dz10; dz21 = zs[2]-zs[1]

            dx02 = -dx20; dx12 = -dx21
            dy02 = -dy20; dy12 = -dy21
            dz02 = -dz20; dz12 = -dz21

            dax = dy10*dz20 - dy20*dz10
            day = dz10*dx20 - dz20*dx10
            daz = dx10*dy20 - dx20*dy10
            da = np.sqrt(dax*dax + day*day + daz*daz)

            nx[tri] += dax/da
            ny[tri] += day/da
            nz[tri] += daz/da

            da = da/2.0
            Arcell += da

            Vcell += (1.0/6.0)*(xs[0]*dax + ys[0]*day + zs[0]*daz)

            # gaussian curvature
            d10 = np.sqrt(dx10*dx10 + dy10*dy10 + dz10*dz10)
            d20 = np.sqrt(dx20*dx20 + dy20*dy20 + dz20*dz20)
            ang0 = np.arccos((dx10*dx20+dy10*dy20+dz10*dz20)/(d10*d20))
            K[tri[0]] -= ang0

            d01 = np.sqrt(dx01*dx01 + dy01*dy01 + dz01*dz01)
            d21 = np.sqrt(dx21*dx21 + dy21*dy21 + dz21*dz21)
            ang1 = np.arccos((dx01*dx21+dy01*dy21+dz01*dz21)/(d01*d21))
            K[tri[1]] -= ang1

            dx1 = xs[0]-xs[2]; dx2 = xs[1]-xs[2]
            dy1 = ys[0]-ys[2]; dy2 = ys[1]-ys[2]
            dz1 = zs[0]-zs[2]; dz2 = zs[1]-zs[2]

            d02 = np.sqrt(dx02*dx02 + dy02*dy02 + dz02*dz02)
            d12 = np.sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12)
            ang2 = np.arccos((dx12*dx02+dy12*dy02+dz12*dz02)/(d12*d02))
            K[tri[2]] -= ang2

            def cot(ang):
                return 1.0/np.tan(ang)
                
            # mean curvature
            alpi = np.arccos((dx01*dx21 + dy01*dy21 + dz01*dz21)/(d01*d21))
            beti = np.arccos((dx02*dx12 + dy02*dy12 + dz02*dz12)/(d02*d12))
            cotalpi0 = cot(alpi)
            cotbeti0 = cot(beti)
            Hx[tri[0]] += dx20*cotalpi0 + dx10*cotbeti0
            Hy[tri[0]] += dy20*cotalpi0 + dy10*cotbeti0
            Hz[tri[0]] += dz20*cotalpi0 + dz10*cotbeti0

            alpi = np.arccos((dx10*dx20 + dy10*dy20 + dz10*dz20)/(d10*d20))
            beti = np.arccos((dx02*dx12 + dy02*dy12 + dz02*dz12)/(d02*d12))
            cotalpi1 = cot(alpi)
            cotbeti1 = cot(beti)
            Hx[tri[1]] += dx21*cotalpi1 + dx01*cotbeti1
            Hy[tri[1]] += dy21*cotalpi1 + dy01*cotbeti1
            Hz[tri[1]] += dz21*cotalpi1 + dz01*cotbeti1

            alpi = np.arccos((dx10*dx20 + dy10*dy20 + dz10*dz20)/(d10*d20))
            beti = np.arccos((dx01*dx21 + dy01*dy21 + dz01*dz21)/(d01*d21))
            cotalpi2 = cot(alpi)
            cotbeti2 = cot(beti)
            Hx[tri[2]] += dx12*cotalpi2 + dx02*cotbeti2
            Hy[tri[2]] += dy12*cotalpi2 + dy02*cotbeti2
            Hz[tri[2]] += dz12*cotalpi2 + dz02*cotbeti2

            if (dx10*dx20+dy10*dy20+dz10*dz20)/(d10*d20) < 0: # angle at 0 is obtuse
                ar[tri[0]] += da/2.0
                ar[tri[1:2]] += da/4.0
            elif (dx01*dx21+dy01*dy21+dz01*dz21)/(d01*d21) < 0: # angle at 0 is obtuse
                ar[tri[1]] += da/2.0
                ar[tri[0]] += da/4.0
                ar[tri[2]] += da/4.0
            elif (dx02*dx12+dy02*dy12+dz02*dz12)/(d02*d12) < 0: # angle at 0 is obtuse
                ar[tri[0:1]] += da/4.0
                ar[tri[2]] += da/2.0
            else:
                ar[tri[0]] += (1.0/8.0)*(d02*d02*cotalpi0 + d01*d01*cotbeti0)
                ar[tri[1]] += (1.0/8.0)*(d21*d21*cotalpi1 + d01*d01*cotbeti1)
                ar[tri[2]] += (1.0/8.0)*(d12*d12*cotalpi2 + d02*d02*cotbeti2)

        H = np.sqrt(Hx*Hx + Hy*Hy + Hz*Hz)/2.0
        n = np.sqrt(nx*nx + ny*ny + nz*nz)
        
        return (nx/n,ny/n,nz/n,ar,K/ar,H/ar,Arcell,Vcell)

    def LaplaceBeltrami(self,Q):
    # calculates the LB of Q on the mesh.
    # See section Mesh laplacian at https://en.wikipedia.org/wiki/Discrete_Laplace_operator
        import numpy as np
        lbQ = 0*self.x3d
        ar = 0*self.x3d
        for tri in self.conn:
            xs = self.x3d[tri]
            ys = self.y3d[tri]
            zs = self.z3d[tri]
            
            dx10 = xs[1]-xs[0]; dx20 = xs[2]-xs[0]
            dy10 = ys[1]-ys[0]; dy20 = ys[2]-ys[0]
            dz10 = zs[1]-zs[0]; dz20 = zs[2]-zs[0]

            dx01 = -dx10; dx21 = xs[2]-xs[1]
            dy01 = -dy10; dy21 = ys[2]-ys[1]
            dz01 = -dz10; dz21 = zs[2]-zs[1]

            dx02 = -dx20; dx12 = -dx21
            dy02 = -dy20; dy12 = -dy21
            dz02 = -dz20; dz12 = -dz21

            dax = dy10*dz20 - dy20*dz10
            day = dz10*dx20 - dz20*dx10
            daz = dx10*dy20 - dx20*dy10
            da = np.sqrt(dax*dax + day*day + daz*daz)/2.0

            d10 = np.sqrt(dx10*dx10 + dy10*dy10 + dz10*dz10)
            d20 = np.sqrt(dx20*dx20 + dy20*dy20 + dz20*dz20)
            ang0 = np.arccos((dx10*dx20+dy10*dy20+dz10*dz20)/(d10*d20))

            d01 = np.sqrt(dx01*dx01 + dy01*dy01 + dz01*dz01)
            d21 = np.sqrt(dx21*dx21 + dy21*dy21 + dz21*dz21)
            ang1 = np.arccos((dx01*dx21+dy01*dy21+dz01*dz21)/(d01*d21))

            dx1 = xs[0]-xs[2]; dx2 = xs[1]-xs[2]
            dy1 = ys[0]-ys[2]; dy2 = ys[1]-ys[2]
            dz1 = zs[0]-zs[2]; dz2 = zs[1]-zs[2]

            d02 = np.sqrt(dx02*dx02 + dy02*dy02 + dz02*dz02)
            d12 = np.sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12)
            ang2 = np.arccos((dx12*dx02+dy12*dy02+dz12*dz02)/(d12*d02))

            def cot(ang):
                return 1.0/np.tan(ang)
                
            alpi = np.arccos((dx01*dx21 + dy01*dy21 + dz01*dz21)/(d01*d21))
            beti = np.arccos((dx02*dx12 + dy02*dy12 + dz02*dz12)/(d02*d12))
            cotalpi0 = cot(alpi)
            cotbeti0 = cot(beti)
            lbQ[tri[0]] += cotalpi0*(Q[tri[2]]-Q[tri[0]]) + cotbeti0*(Q[tri[1]]-Q[tri[0]])

            alpi = np.arccos((dx10*dx20 + dy10*dy20 + dz10*dz20)/(d10*d20))
            beti = np.arccos((dx02*dx12 + dy02*dy12 + dz02*dz12)/(d02*d12))
            cotalpi1 = cot(alpi)
            cotbeti1 = cot(beti)
            lbQ[tri[1]] += cotalpi1*(Q[tri[2]]-Q[tri[1]]) + cotbeti1*(Q[tri[0]]-Q[tri[1]])

            alpi = np.arccos((dx10*dx20 + dy10*dy20 + dz10*dz20)/(d10*d20))
            beti = np.arccos((dx01*dx21 + dy01*dy21 + dz01*dz21)/(d01*d21))
            cotalpi2 = cot(alpi)
            cotbeti2 = cot(beti)
            lbQ[tri[2]] += cotalpi2*(Q[tri[1]]-Q[tri[2]]) + cotbeti2*(Q[tri[0]]-Q[tri[2]])

            if (dx10*dx20+dy10*dy20+dz10*dz20)/(d10*d20) < 0: # angle at 0 is obtuse
                ar[tri[0]] += da/2.0
                ar[tri[1:2]] += da/4.0
            elif (dx01*dx21+dy01*dy21+dz01*dz21)/(d01*d21) < 0: # angle at 0 is obtuse
                ar[tri[1]] += da/2.0
                ar[tri[0]] += da/4.0
                ar[tri[2]] += da/4.0
            elif (dx02*dx12+dy02*dy12+dz02*dz12)/(d02*d12) < 0: # angle at 0 is obtuse
                ar[tri[0:1]] += da/4.0
                ar[tri[2]] += da/2.0
            else:
                ar[tri[0]] += (1.0/8.0)*(d02*d02*cotalpi0 + d01*d01*cotbeti0)
                ar[tri[1]] += (1.0/8.0)*(d21*d21*cotalpi1 + d01*d01*cotbeti1)
                ar[tri[2]] += (1.0/8.0)*(d12*d12*cotalpi2 + d02*d02*cotbeti2)

        return lbQ/(2.0*ar)

    def bendingForce(self,H,K,nx,ny,nz,kb):
        import numpy as np
        # kb is the bending modulus of the surface

        # assuming that the intrinsic curvature is zero
        # See equation (5) in Farutin etal., 3D numerical simulations of vesicle and inextensible capsule dynamics, 2014, J. Comp. Phys.

        # We are considering the cell to be a compressible capsule. So we can ingore the membrane tension term from the bending force expression. 
        lbH = self.LaplaceBeltrami(Q)
        fmag = -kb*(2*H*(2*H*H-2*K) + 2*lbH)

        return (fmag*nx,fmag*ny,fmag*nz)

    def discreteMetricTensor(self):
        import numpy as np
        g11 = np.zeros((len(self.conn),3))
        g22 = np.zeros((len(self.conn),3))
        g12 = np.zeros((len(self.conn),3))

        xt = self.x3d
        yt = self.y3d
        zt = self.z3d
        tri0 = self.conn[:,0]
        tri1 = self.conn[:,1]
        tri2 = self.conn[:,2]

        g11[:,0] = (xt[tri1] - xt[tri0])*(xt[tri1] - xt[tri0]) + (yt[tri1] - yt[tri0])*(yt[tri1] - yt[tri0]) + (zt[tri1] - zt[tri0])*(zt[tri1] - zt[tri0])
        g22[:,0] = (xt[tri2] - xt[tri0])*(xt[tri2] - xt[tri0]) + (yt[tri2] - yt[tri0])*(yt[tri2] - yt[tri0]) + (zt[tri2] - zt[tri0])*(zt[tri2] - zt[tri0])
        g12[:,0] = (xt[tri1] - xt[tri0])*(xt[tri2] - xt[tri0]) + (yt[tri1] - yt[tri0])*(yt[tri2] - yt[tri0]) + (zt[tri1] - zt[tri0])*(zt[tri2] - zt[tri0])

        g11[:,1] = (xt[tri2] - xt[tri1])*(xt[tri2] - xt[tri1]) + (yt[tri2] - yt[tri1])*(yt[tri2] - yt[tri1]) + (zt[tri2] - zt[tri1])*(zt[tri2] - zt[tri1])
        g22[:,1] = (xt[tri0] - xt[tri1])*(xt[tri0] - xt[tri1]) + (yt[tri0] - yt[tri1])*(yt[tri0] - yt[tri1]) + (zt[tri0] - zt[tri1])*(zt[tri0] - zt[tri1])
        g12[:,1] = (xt[tri2] - xt[tri1])*(xt[tri0] - xt[tri1]) + (yt[tri2] - yt[tri1])*(yt[tri0] - yt[tri1]) + (zt[tri2] - zt[tri1])*(zt[tri0] - zt[tri1])

        g11[:,2] = (xt[tri0] - xt[tri2])*(xt[tri0] - xt[tri2]) + (yt[tri0] - yt[tri2])*(yt[tri0] - yt[tri2]) + (zt[tri0] - zt[tri2])*(zt[tri0] - zt[tri2])
        g22[:,2] = (xt[tri1] - xt[tri2])*(xt[tri1] - xt[tri2]) + (yt[tri1] - yt[tri2])*(yt[tri1] - yt[tri2]) + (zt[tri1] - zt[tri2])*(zt[tri1] - zt[tri2])
        g12[:,2] = (xt[tri0] - xt[tri2])*(xt[tri1] - xt[tri2]) + (yt[tri0] - yt[tri2])*(yt[tri1] - yt[tri2]) + (zt[tri0] - zt[tri2])*(zt[tri1] - zt[tri2])

        return (g11,g22,g12)

    def elasticForce(self,g110,g220,g120,detg0,nx,ny,nz,ar,mus,mue=0,capsule='NeoHookean'):
        import numpy as np
        # mus is the shear modulus of the capsule material
        # mue is the bulk modulus (used in Skalak model)
        # g0 is the metric tensor for the reference configuration

        # See section-7.3 in Farutin etal., 3D numerical simulations of vesicle and inextensible capsule dynamics, 2014, J. Comp. Phys.

        fxe = 0.0*self.x3d
        fye = 0.0*self.y3d
        fze = 0.0*self.z3d

        (g11c,g22c,g12c) = self.discreteMetricTensor()
        detgc = g11c*g22c - g12c*g12c# determinant of metric tensor
        I1 = (g11c*g220 + g22c*g110 - 2.0*g12c*g120)/detg0 - 2.0 # first invariant of deformation
        I2 = detgc/detg0 - 1.0                                  # second invariant of deformation

        if capsule is 'NeoHookean': # Neo-Hookean
            fxe = 1
#        if capsule is 'Skalak': # Skalak

        return 1
            
    

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''
    import numpy as np

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.3*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def generate3DTissue(Nx,Ny,noisx,noisy,direc,dl,dh,nl,cyl,cylCoverage,makePlotHere,datadir):
    import numpy as np
    from generateCells import generateCellsPeriodic
    from plottingUtils import plotCellsTriangular
    from generateCellMembrane import generateCellMembrane
    from vtkReaderWriter import writeVTK_cell, writeVTK_tissue


    # nuclei
    from genSphere import genNucleus
    nrad = 0.1
    (xn,yn,zn,conn_n) = genNucleus(radius=nrad,xc=0.0,yc=0.0,zc=0.0)
    print( 'Area and volume of a nucleus: ', 4.0*np.pi*nrad*nrad, 4.0*np.pi*nrad*nrad*nrad/3.0)
    ###########

    if makePlotHere:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        fig = plt.figure(); ax = fig.gca(projection='3d')

    v,cells,nbs,Lx,Ly = generateCellsPeriodic(Nx,Ny,noisx,noisy,direc=direc,scalingx=1.5, scalingy=1.2)
    compfac = 1.0; v[:,0] /= compfac; Lx /= compfac
    (xcells,ycells) = generateCellMembrane(v,cells,dl=dl,makeplot=0)

    R = Lx/(2*np.pi*cylCoverage)
    
    scList = [] # cell membrae
    ncList = [] # nuclei
    for ci in range(len(xcells)):
        sc = singleCell(xcell=xcells[ci],ycell=ycells[ci],Lx=Lx,dh=dh,nl=nl,cyl=cyl,R=R)
        nc = singleCell(x3d=xn+np.mean(sc.x3d),y3d=yn+np.mean(sc.y3d),z3d=zn+np.mean(sc.z3d),conn=conn_n)
        scList += [sc];
        ncList += [nc];
        writeVTK_cell(sc,datadir+'cellsWT/cell'+str(ci)+'.vtk');
        writeVTK_cell(nc,datadir+'cellsWT/nucleus'+str(ci)+'.vtk');
        if makePlotHere:
            for tri in sc.conn:
                x = [sc.x3d[tri[0]],sc.x3d[tri[1]],sc.x3d[tri[2]]]
                y = [sc.y3d[tri[0]],sc.y3d[tri[1]],sc.y3d[tri[2]]]
                z = [sc.z3d[tri[0]],sc.z3d[tri[1]],sc.z3d[tri[2]]]

    #            ax.plot(x,y,z,'k')
                verts = [np.array(list(zip(x, y,z)))]
                ax.add_collection3d(Poly3DCollection(verts,facecolor='r',alpha=0.7))
#    set_axes_equal(ax)
    print( 'Total number of nodes: ',(len(sc.x3d)+len(xn))*len(cells)); #exit()
    if makePlotHere:
        ax.set_xlim([-5,5])
        ax.set_ylim([-5,5])
        ax.set_zlim([-5,5])
        plt.show()

    writeVTK_tissue(scList,datadir+'tissue_memWT.vtk')
    writeVTK_tissue(ncList,datadir+'tissue_nucWT.vtk')
    return scList

if __name__ == "__main__":
    import sys
    sys.path.append('../')

    makePlotHere = 0
    Nx = 2 # width and
    Ny = 12 # length of the tissue in the epithelial plane
    noisx = 0.2
    noisy = 0.2
    direc = 1
    dl = 0.1 # length discretization in plane
    dh = 0.05 # height discretzation 
    nl = 2 # number of height elements
    cylindrical=True
    cylCoverage = 1.0
    datadir = 'cellsGeom/'
    scList = generate3DTissue(Nx,Ny,noisx,noisy,direc,dl,dh,nl,cylindrical,cylCoverage,makePlotHere,datadir)
    exit()
    
    
