def writeVTK_cell(sc,filename,geomProps=False):
    import numpy as np
    
    x = sc.x3d
    y = sc.y3d
    z = sc.z3d
    conn = sc.conn
    conn2 = np.zeros((len(conn),4))
    conn2[:,0] = 3
    conn2[:,1:] = conn
    
    Np = len(x)
    fl = open(filename,'w')
    header = '# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS '+str(Np)+' float\n'
    fl. write(header)
    fl.close()
    
    fl=open(filename,'a')
    np.savetxt(fl,np.array(list(zip(x,y,z))),fmt='%10.7f')

    size = 4*len(conn)
    fl.write('POLYGONS '+str(len(conn))+' '+str(size)+' \n')
    np.savetxt(fl,conn2,fmt='%i')

    if geomProps:
        from timeit import default_timer as timer
        startt = timer()
        (nx,ny,nz,ar,K,H,A,V) = sc.surfaceNormals_TotalArea_GaussianMeanCurvatures_Volume()
        print('Time taken in curvature etc calculation: ', timer()-startt)
        header = 'POINT_DATA '+str(Np)+'\nVECTORS normals float\n'
        fl.write(header)        
        np.savetxt(fl,np.array(list(zip(nx,ny,nz))),fmt='%10.7f')
        header = 'SCALARS gaussian_curvature float\nLOOKUP_TABLE default\n'
        fl.write(header)        
        np.savetxt(fl,K,fmt='%10.7f')
        header = 'SCALARS mean_curvature float\nLOOKUP_TABLE default\n'
        fl.write(header)        
        np.savetxt(fl,H,fmt='%10.7f')

    fl.close()

def writeVTK_tissue(scList,filename):
    import numpy as np
    
    x = np.empty(0)
    y = np.empty(0)
    z = np.empty(0)
    nxs = np.empty(0); nys = np.empty(0); nzs = np.empty(0)
    Ks = np.empty(0); Hs = np.empty(0); lbHs = np.empty(0)
    fxs = np.empty(0); fys = np.empty(0); fzs = np.empty(0)
    for sci,sc in enumerate(scList):
        (nx,ny,nz,ar,K,H,A,V) = sc.surfaceNormals_TotalArea_GaussianMeanCurvatures_Volume()
        lbH = sc.LaplaceBeltrami(H)
#        (fx,fy,fz) = sc.intracellularForce(H,K,lbH,nx,ny,nz)
        print( 'Cell area and volume: ',A,V)
        x = np.concatenate((x,sc.x3d),axis=0)
        y = np.concatenate((y,sc.y3d),axis=0)
        z = np.concatenate((z,sc.z3d),axis=0)
        nxs = np.concatenate((nxs,nx),axis=0)
        nys = np.concatenate((nys,ny),axis=0)
        nzs = np.concatenate((nzs,nz),axis=0)
        Ks = np.concatenate((Ks,K),axis=0)
        Hs = np.concatenate((Hs,H),axis=0)
        lbHs = np.concatenate((lbHs,lbH),axis=0)

#        fxs = np.concatenate((fxs,fx),axis=0)
#        fys = np.concatenate((fys,fy),axis=0)
#        fzs = np.concatenate((fzs,fz),axis=0)

        conn = sc.conn
        conn2 = np.zeros((len(conn),4))
        conn2[:,0] = 3
        conn2[:,1:] = conn
        if sci==0: # first cell
            connt = 1*conn2
        else:         
            conn2[:,1:] += len(x) - len(sc.x3d) 
            connt = np.concatenate((connt,conn2),axis=0)

    print( 'Total number of faces: ', len(connt))
    Np = len(x)
    fl = open(filename,'w')
    header = '# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS '+str(Np)+' float\n'
    fl. write(header)
    fl.close()
    
    fl=open(filename,'a')
    np.savetxt(fl,np.array(list(zip(x,y,z))),fmt='%10.7f')

    size = 4*len(connt)
    fl.write('POLYGONS '+str(len(connt))+' '+str(size)+' \n')
    np.savetxt(fl,connt,fmt='%i')

    header = 'POINT_DATA '+str(Np)+'\nVECTORS normals float\n'
    fl.write(header)        
    np.savetxt(fl,np.array(list(zip(nxs,nys,nzs))),fmt='%10.7f')

    header = 'SCALARS gaussian_curvature float\nLOOKUP_TABLE default\n'
    fl.write(header)        
    np.savetxt(fl,Ks,fmt='%10.7f')

    header = 'SCALARS mean_curvature float\nLOOKUP_TABLE default\n'
    fl.write(header)        
    np.savetxt(fl,Hs,fmt='%10.7f')

    header = 'SCALARS laplaceBeltrami_mean_curvature float\nLOOKUP_TABLE default\n'
    fl.write(header)        
    np.savetxt(fl,lbHs,fmt='%10.7f')
    
    fl.close()

def readVTK_cell(filename):
    import numpy as np

    fl = open(filename)
    lines = fl.readlines()
    for lind,line in enumerate(lines):
        if 'POINTS' in line:
            Np = int(line.split(' ')[1])
            stp = lind+1

        if 'POLYGONS' in line:
            Nc = int(line.split(' ')[1])
            stc = lind+1
            
    points = np.loadtxt(filename,skiprows=stp,max_rows=Np)
    conn = np.loadtxt(filename,skiprows=stc,max_rows=Nc).astype(int)
    conn = conn[:,1:]

    x3d = points[:,0]
    y3d = points[:,1]
    z3d = points[:,2]

    return (x3d,y3d,z3d,conn)

if __name__ == "__main__":
    import numpy as np
    from generate3DCylindrical import *

#    (x3d,y3d,z3d,conn) = readVTK_cell('shape_10000.vtk')
#    sc = singleCell(x3d=x3d,y3d=y3d,z3d=z3d,conn=conn)
#    writeVTK_cell(sc,'cellsGeom/cells/testcell.vtk',geomProps=1)
#    exit/()

    
    (x3d,y3d,z3d,conn) = readVTK_cell('cellsGeom/tissue_mem.vtk')
    sc = singleCell(x3d=x3d,y3d=y3d,z3d=z3d,conn=conn)
    writeVTK_cell(sc,'cellsGeom/tissue_mem2.vtk',geomProps=1)

