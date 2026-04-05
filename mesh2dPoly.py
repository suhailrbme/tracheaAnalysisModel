def create_mesh(boundary_points,edgelength=1.0, max_area=0.05):
    import numpy as np
    import meshpy.triangle
    # dimensions of the rectangle
    lx = edgelength
    ly = edgelength

    info = meshpy.triangle.MeshInfo()
    info.set_points(boundary_points)

    def _round_trip_connect(start, end):
        result = []
        for i in range(start, end):
            result.append((i, i + 1))
        result.append((end, start))
        return result

    info.set_facets(_round_trip_connect(0, len(boundary_points) - 1))

    def _needs_refinement(vertices, area):
        return bool(area > max_area)

    meshpy_mesh = meshpy.triangle.build(info, refinement_func=_needs_refinement)

    # append column
    pts = np.array(meshpy_mesh.points)
    points = np.c_[pts[:, 0], pts[:, 1], np.zeros(len(pts))]

    return points, np.array(meshpy_mesh.elements)

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    
    boundary_points = np.array([[2*np.cos(th),np.sin(th)] for th in np.linspace(0,2*np.pi,20)[:-1]])
    points, cells = create_mesh(boundary_points)

    for cl in cells:
        x = points[cl,0].tolist() + [points[cl[0],0]]
        y = points[cl,1].tolist() + [points[cl[0],1]]
        plt.plot(x,y,'k-')
    plt.plot(points[:19,0],points[:19,1],'k-')
    plt.axis('equal')
    plt.show()
