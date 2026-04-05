def genNucleus(radius=1.0,xc=0.0,yc=0.0,zc=0.0):
    import trimesh.creation as trim

    mesh = trim.icosphere(1) # 1 is subdivision
    vers = mesh.vertices
    cells = mesh.faces

    x = radius*vers[:,0] + xc
    y = radius*vers[:,1] + yc
    z = radius*vers[:,2] + zc
    
    return (x,y,z,cells)
    
def genRandomPointsSphere(N):
    import numpy as np

    d = 2.0/np.sqrt(N) # approximate radius of a patch on sphere
    ths = np.zeros(N) # goes from 0 to 2 pi
    phs = np.zeros(N) # goes from 0 to pi
    for i in range(N):
        while(True):
            th = 2*np.pi*np.random.rand(1)
            ph = np.arccos(2*np.random.rand(1)-1)
            if i == 0:
                break
            else:
                dijs = np.sqrt(2.0-2.0*( np.cos(ph)*np.cos(phs[:i]) + np.sin(ph)*np.sin(phs[:i])*np.cos(th-ths[:i]) ))
                if sum((dijs>1.3*d).astype(int)) == i:
                    break
        ths[i] = th
        phs[i] = ph

    x = np.sin(phs)*np.cos(ths)
    y = np.sin(phs)*np.sin(ths)
    z = np.cos(phs)
    return (x,y,z)

def sphericalBasalLayer(radius=3.0,Ncell=4):
    import numpy as np
    import trimesh.creation as trim
    import matplotlib.pyplot as plt 
    from scipy.spatial import SphericalVoronoi

    (x,y,z) = genRandomPointsSphere(Ncell)
    points = np.zeros((Ncell,3))
    points[:,0] = x;
    points[:,1] = y;
    points[:,2] = z;
    sv = SphericalVoronoi(points, radius=1, center=np.array([0.0,0.0,0.0]))
    sv.sort_vertices_of_regions()

#    from mpl_toolkits.mplot3d import proj3d
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    # plot the unit sphere for reference (optional)
#    u = np.linspace(0, 2 * np.pi, 100)
#    v = np.linspace(0, np.pi, 100)
#    x = np.outer(np.cos(u), np.sin(v))
#    y = np.outer(np.sin(u), np.sin(v))
#    z = np.outer(np.ones(np.size(u)), np.cos(v))
#    ax.plot_surface(x, y, z, color='y', alpha=1.0)
#    # plot generator points
#    ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')
#    
#    for region in sv.regions:
#        n = len(region)
#        x = sv.vertices[region,0]; x = x.tolist() + [x[0]]
#        y = sv.vertices[region,1]; y = y.tolist() + [y[0]]
#        z = sv.vertices[region,2]; z = z.tolist() + [z[0]]
#        ax.plot(x,y,z)
#    fig.set_size_inches(10, 10)
#    plt.show()

    return (radius*sv.vertices,sv.regions)

if __name__ == "__main__":
    from generate3DCylindrical import *
    from vtkReaderWriter import writeVTK_cell
    (x,y,z,conn) = genNucleus(radius=1.0,xc=0.0,yc=0.0,zc=0.0)
    

    sphericalBasalLayer(radius=3.0)
    
