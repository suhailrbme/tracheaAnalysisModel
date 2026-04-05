def smooth(xi,yi,iters=3,frac=0.5,zi=None):
    import numpy as np
    
    d = 0.005
    xc = np.mean(xi)
    yc = np.mean(yi)

    xi = (1.0-d)*xi + d*xc
    yi = (1.0-d)*yi + d*yc
        
    if zi is not None:    
        zc = np.mean(zi)
        zi = (1.0-d)*zi + d*zc

    for i in range(iters):
        xi = frac*xi + (0.5-frac/2)*np.roll(xi,-1) + (0.5-frac/2)*np.roll(xi,+1)
        yi = frac*yi + (0.5-frac/2)*np.roll(yi,-1) + (0.5-frac/2)*np.roll(yi,+1)
        if zi is not None:            
            zi = frac*zi + (0.5-frac/2)*np.roll(zi,-1) + (0.5-frac/2)*np.roll(zi,+1)
    if zi is None:
        return (xi,yi)
    else:
        return (xi,yi,zi)

def generateCellMembrane(v,cells,dl=0.1,makeplot=1):
    import numpy as np
    import matplotlib.pyplot as plt

    xcells = []
    ycells = []    
    for c in cells:
        x = v[c,0]; xs = x.tolist() + [x[0]]
        y = v[c,1]; ys = y.tolist() + [y[0]]
        
        xinow = []
        yinow = []
        for vi in range(len(xs)-1):
            l = np.sqrt((xs[vi]-xs[vi+1])**2 + (ys[vi]-ys[vi+1])**2)
            t = np.linspace(0,1,int(l/dl))[:-1]
            xi = (1.0-t)*xs[vi] + t*xs[vi+1];
            yi = (1.0-t)*ys[vi] + t*ys[vi+1]
    
            xinow += xi.tolist()
            yinow += yi.tolist()
        xinow = np.array(xinow)
        yinow = np.array(yinow)

#        (xinow,yinow) = smooth(xinow,yinow,iters=1)
        
        xcells.append(xinow)
        ycells.append(yinow)

        if makeplot:
            plt.fill(xinow,yinow,'-')
    if makeplot:
        plt.axis('equal')
        plt.show()
    return (xcells,ycells)

def generateCellMembrane3D(v,cells,dl=0.1,makeplot=1):
    import numpy as np
    import matplotlib.pyplot as plt

    if makeplot:
        from mpl_toolkits.mplot3d import proj3d
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    xcells = []
    ycells = []    
    zcells = []    
    for c in cells:
        x = v[c,0]; xs = x.tolist() + [x[0]]
        y = v[c,1]; ys = y.tolist() + [y[0]]
        z = v[c,2]; zs = z.tolist() + [z[0]]
        
        xinow = []
        yinow = []
        zinow = []
        for vi in range(len(xs)-1):
            l = np.sqrt((xs[vi]-xs[vi+1])**2 + (ys[vi]-ys[vi+1])**2 +  + (zs[vi]-zs[vi+1])**2)
            t = np.linspace(0,1,int(l/dl))[:-1]
            xi = (1.0-t)*xs[vi] + t*xs[vi+1];
            yi = (1.0-t)*ys[vi] + t*ys[vi+1]
            zi = (1.0-t)*zs[vi] + t*zs[vi+1]
        
            xinow += xi.tolist()
            yinow += yi.tolist()
            zinow += zi.tolist()

        xinow = np.array(xinow)
        yinow = np.array(yinow)
        zinow = np.array(zinow)

        (xinow,yinow,zinow) = smooth(xinow,yinow,iters=1,frac=0.5,zi=zinow)
        if makeplot:
            ax.plot(xinow,yinow,zinow,'-')

        xcells.append(xinow)
        ycells.append(yinow)
        zcells.append(zinow)

    if makeplot:
        fig.set_size_inches(4, 4)
        plt.show()
    return (xcells,ycells,zcells)

