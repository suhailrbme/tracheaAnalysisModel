def plotCellsTriangular(v,cells,ts):
    import sys
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    sys.path.append('../../../')
    from geomProperties import PolyArea

    minarea = 0.25
    maxarea = 0.6

    patches = []
    
    jt = plt.get_cmap('jet')# colormap
    colors = []
    for c in cells:
        c.append(c[0])
        x = v[c,0]
        y = v[c,1]

        ar = PolyArea(x[:-1],y[:-1])
        cellcolor = (ar-minarea)/(maxarea-minarea)
        if cellcolor < 0.0:
            cellcolor = 0.0
        if cellcolor > 1.0:
            cellcolor = 1.0
        cellcolor = jt(cellcolor)
        xy = []
        for nn in range(len(x)):
            xy.append([y[nn],x[nn]])
        polygon = Polygon(xy,closed=True)
        patches.append(polygon)
        #plt.gca().add_artist(polygon)
        colors.append(cellcolor)

    p = PatchCollection(patches,cmap=matplotlib.cm.jet,alpha=1)
        #plt.plot(x,y,'k',linewidth=0.4)
    plt.gca().add_collection(p)
    p.set_color(colors)
    p.set_edgecolor('black')
    plt.gca().autoscale_view()

    #plt.axis('off')
    plt.axis('equal')
    #plt.ylim([-45,45])
    plt.savefig('withFj_cellsize.eps',format='eps', bbox_inches='tight')
    plt.show()

def plotCellsTriangularEdges(v,cells,x0,y0,p,th,q,ph,protein,losscolor='y',savefile=''):
    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    plt.figure(figsize=(5,5))
    for i, c in enumerate(cells):
        xv = v[c,0]
        yv = v[c,1]
        plt.plot(xv.tolist()+[xv[0]],yv.tolist()+[yv[0]],'k',linewidth=0.4)
        if protein[i] == 0:
            plt.fill(xv,yv,losscolor,alpha=0.2)

    wd = 0.002; scl = 20
    #q2 = plt.quiver(x0,y0,p*np.cos(th),p*np.sin(th),angles='xy',width=wd,scale=scl,color='r')
    q2 = plt.quiver(x0,y0,q*np.cos(ph),q*np.sin(ph),angles='xy',width=wd,scale=scl,color='g')

    plt.axis('off')
    plt.axis('equal')
    plt.margins(0.1,0.1)
    #plt.ylim([-25,25])
    if savefile != '':
        plt.savefig(savefile,format='png', bbox_inches='tight', pad_inches = 0)
    plt.show()

def plotCellsHeterodimer(v,cells,FtDs,DsFt):
    import sys
    import numpy as np
    
    for i, c in enumerate(cells):
        c.append(c[0])
        xv = v[c,1]
        yv = v[c,0]
        plt.plot(xv,yv,'k',linewidth=0.4)
    
    for i, c in enumerate(cells):
        c.append(c[0])
        xv = v[c,1]
        yv = v[c,0]
        
        fd = np.array(FtDs[i])
        df = np.array(DsFt[i])
        
        tot = fd+df
        mn = np.mean(tot)
        if min(tot)>0.2:
            hig = np.where(tot>0)
        else:
            hig = np.where(tot>mn+0.02)
        
        xv.tolist().append(xv[0])
        yv.tolist().append(yv[0])
        for h in hig:
            plt.plot([xv[h],xv[h+1]],[yv[h],yv[h+1]],'r',linewidth=0.4)
        

   
    #plt.gca().invert_yaxis()
    #plt.gca().set_aspect('equal',adjustable='box')
    plt.axis('equal')
    #plt.axis('off')
    plt.savefig('withFj.eps',format='eps', bbox_inches='tight')
    plt.show()


def plotCellsFill(v,cells,prot):
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    from matplotlib.collections import LineCollection

    ## plotting Ft, Ds, Fj
    fig,ax = plt.subplots()
    patches = []
    jt = plt.get_cmap('Greens')
    colors = []
    for idx,c in enumerate(cells):
        clr = prot[idx]
        cellcolor = jt(clr)
        xy = []
        for nn in range(len(c)):
            xy.append([v[c[nn],1],v[c[nn],0]])
        polygon = Polygon(xy,closed=True)
        patches.append(polygon)
        #ax.add_artist(polygon)
        colors.append(cellcolor)
        
    p = PatchCollection(patches,cmap=matplotlib.cm.jet,alpha=1)

    ax.add_collection(p)
    p.set_color(colors)
    p.set_edgecolor('black')
    ax.autoscale_view()
    #plt.gca().invert_yaxis()
    #plt.gca().set_aspect('equal',adjustable='box')
    plt.axis('equal')
    #plt.axis('off')
    plt.show()

if __name__ == "__main__":
    import sys
    import pickle
    import numpy as np
    import matplotlib.pyplot as plt
    sys.path.append('../../../FtDsFj/')
    
    from getFtDsFj import getFtDsFj
    from equilibriumHeterodimer import equilibriumHeterodimer
    from nextHeterodimer import nextHeterodimer
    from heterodimer2Adhesion import heterodimer2Adhesion

    datadir = 'data/'
    ts = int(sys.argv[1])
    print( ts)
#    for ts in range(999,1000):
    v = pickle.load(open(datadir+'ve'+str(10*ts)+'.dat', 'rb'))
    cells = pickle.load(open(datadir+'cells'+str(10*ts)+'.dat', 'rb'))
    nbs = pickle.load(open(datadir+'nbs'+str(10*ts)+'.dat', 'rb'))
#    plotCellsTriangularEdges(v,cells,ts)
#    plotCellsTriangular(v,cells,ts)
#    exit()

    v0 = pickle.load(open(datadir+'ve_init.dat', 'rb'))
    cells0 = pickle.load(open(datadir+'cells_init.dat', 'rb'))
    nbs0 = pickle.load(open(datadir+'nbs_init.dat', 'rb'))
    (Ft,Ds,Fj) = getFtDsFj(v0,cells0,40,25)
    (FtDs,DsFt) = equilibriumHeterodimer(v,cells,nbs,Ft,Ds,Fj)
    plotCellsHeterodimer(v,cells,FtDs,DsFt)
