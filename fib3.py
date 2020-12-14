from stl import mesh
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import math
phi = 1.61803398875
'''3D Fibonacci, yes pls!'''
def MakeVertices(n, T, r):
    # return: V
    V = []
    x = 0
    y = 0
    z = 0
    l = 1
    for t in range(T):
        dx = -1*(l*phi**(t/10))*math.cos(phi*t/10)
        dy = -1*(l*phi**(t/10))*math.sin(phi*t/10)
        dz = 0
        x = dx + x
        y = dy + y
        z = dz + z
        #r = r + .1

        #x_tt = -1*(l*phi**((t+1)/10))*math.cos(phi*(t+1)/10) + x
        #y_tt = -1*(l*phi**((t+1)/10))*math.sin(phi*(t+1)/10) + y
        #dxt = x_tt - x
        #dyt = y_tt - y
        theta = math.atan2(dy, dx)+90*math.pi/180
        theta2 = math.atan2(dz,dx)

        for i in range(n):
            ang = (360/n * i)*math.pi/180
            v = np.array([r*math.cos(ang), 0, r*math.sin(ang)])

            A = np.array([
                [math.cos(theta), -math.sin(theta), 0],
                [math.sin(theta), math.cos(theta), 0],
                [0, 0, 1]
            ])
            v =np.matmul(A,v)

            A = np.array([
                [math.cos(theta), -math.sin(theta), 0],
                [0, 0, 1],
                [math.sin(theta), math.cos(theta), 0]

            ])

            #v =np.matmul(A,v)
            v = v + np.array([x, y, z])
            V.append(v)
    return V

n = 6
T = 120
r = 100
# Define the vertices of the curve
V = MakeVertices(n, T, r)
vertices = np.array(V)
# Define the 12 triangles composing the cube
#faces = [[x for x in range(0, n)],[y for y in range(n,2*n)]]
def MakeFaces(n, T):
    faces = []
    print('faces ',faces)
    for t in range(T-1):
        for j in range(2*n):
            i = j%n

            print(i)
            if j > i:
                #2nd
                if i+1 == n:
                    f = [t*n+i, t*n, (t+1)*n]
                else:
                    f = [t*n+i, t*n+i+1, (t+1)*n+i+1]
            else:
                #first
                if i+n+1 == 2*n:
                    f = [t*n+i, (t+1)*n+i, (t+1)*n]
                else:
                    f = [t*n+i, (t+1)*n+i, (t+1)*n+i+1]
            #f = []
            faces.append(f)
    return faces
faces = MakeFaces(n, T)
faces = np.array(faces)

# Create the mesh
curve = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        curve.vectors[i][j] = vertices[f[j],:]

# Write the mesh to file "cube.stl"
curve.save('models/fib.stl')
# Create a new plot
figure = pyplot.figure()
axes = mplot3d.Axes3D(figure)

# Load the STL files and add the vectors to the plot
curve_mesh = mesh.Mesh.from_file('models/fib.stl')
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(curve_mesh.vectors))

# Auto scale to the mesh size
scale = curve_mesh.points.flatten()
axes.auto_scale_xyz(scale, scale, scale)

# Show the plot to the screen
pyplot.show()
