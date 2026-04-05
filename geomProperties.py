def PolyArea(x,y):    
    import numpy as np
    # x and y should NOT be closed

    area = 0.0
    for i in range(len(x)-1):
        area += x[i]*y[i+1] - x[i+1]*y[i]
    area += x[-1]*y[0] - x[0]*y[-1]
    area = area/2.0
    
    return area

def PolyPerimeter(x,y):    
    import numpy as np
    # x and y should NOT be closed

    perimeter = 0.0
    for i in range(len(x)-1):
        perimeter += np.sqrt((x[i]-x[i+1])*(x[i]-x[i+1]) + (y[i]-y[i+1])*(y[i]-y[i+1]))
    perimeter += np.sqrt((x[-1]-x[0])*(x[-1]-x[0]) + (y[-1]-y[0])*(y[-1]-y[0]))
    
    return perimeter

def cellCenters(v,cells):
    import numpy as np
    # obtain the initial position of cell centers
    x0 = []
    y0 = []
    for c in cells:
        x0.append(np.mean(v[c,0]))
        y0.append(np.mean(v[c,1]))
    x0 = np.array(x0)
    y0 = np.array(y0)

    return (x0,y0)
