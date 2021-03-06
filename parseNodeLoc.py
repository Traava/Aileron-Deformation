import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

node = np.genfromtxt('B737.inp', delimiter = ',', skip_header = 9, skip_footer = 14594-6598)
elem = np.genfromtxt('B737.inp', delimiter = ',', skip_header = 6598, skip_footer = 14594-13233)

nodes = {}
elems = {}
elem_dict = {}
elem_lst = []

for i in range(len(node)):
    nodes[node[i,0]] = node[i,1:]
for i in range(len(elem)):
    elems[elem[i,0]] = elem[i,1:]

elemxz = []

for i in elems:
    
    coords = []
    nod = elems[i]
    for n in nod:
        coords.append(nodes[n])
    coords = np.array(coords)
    x = sum(coords[:,0])/4
    y = sum(coords[:,1])/4
    z = sum(coords[:,2])/4
    
    if i == 7:
        print(coords[:,1])
        print(y)
    
    elem_dict[i] = np.array([x,y,z])
    elem_lst.append(np.array([x,y,z]))
    
    elemxz.append(np.array([x,z]))

elemxz = np.array(elemxz)
elem_lst = np.array(elem_lst)

from scipy.spatial import distance

def closest_node(node, nodes):
    elem_num = distance.cdist([node], nodes).argmin() +1
    return elem_num




    
#x = np.array([[elem_coords[k][0]] for k in elem_coords])
#y = np.array([[elem_coords[k][1]] for k in elem_coords])
#z = np.array([[elem_coords[k][2]] for k in elem_coords])




cases = ['Bending', 'Jam-Bent','Jam-Straight']

f = open('B737.rpt')

data = [[],[],[],[],[],[] ]
pos=0
sizes = [5778, 856,5778, 856,5778, 856]
for line in f:
    if line.startswith('   Element Label'):
        next(f)
        next(f)
        while len(data[pos]) < sizes[pos]:
            data[pos].append(np.array(next(f).split()))
        pos +=1
        
#5577
data = np.array(data)

elem_vm_dict = {}

for i in data[0]:
    elem_vm_dict[float(i[0])] = sum([float(i[2]),float(i[3])])/2
    
    




grain = 100

x = np.linspace(10, 2650, grain)
z = np.linspace(-490, 101, grain)
y = np.array([elem_vm_dict.get(closest_node([i,j],elemxz), 0) for j in z for i in x])

print(len(y))

X, Z = np.meshgrid(x, z)
Y = y.reshape(grain, grain)

plt.pcolor(X, Z, Y)
plt.show()















#
#plt.contour([a[:,1],a[:,2]], a[:,3])
#plt.plot()






## converts quad elements into tri elements
#def quads_to_tris(quads):
#    tris = [[None for j in range(3)] for i in range(2*len(quads))]
#    for i in range(len(quads)):
#        j = 2*i
#        n0 = quads[i][0]
#        n1 = quads[i][1]
#        n2 = quads[i][2]
#        n3 = quads[i][3]
#        tris[j][0] = n0
#        tris[j][1] = n1
#        tris[j][2] = n2
#        tris[j + 1][0] = n2
#        tris[j + 1][1] = n3
#        tris[j + 1][2] = n0
#    return tris
#
## plots a finite element mesh
#def plot_fem_mesh(nodes_x, nodes_y, elements):
#    for element in elements:
#        x = [nodes_x[element[i]] for i in range(len(element))]
#        y = [nodes_y[element[i]] for i in range(len(element))]
#        plt.fill(x, y, edgecolor='black', fill=False)
#
### FEM data
#nodes_x = [0.0, 0.5,1.,0.0, 0.5,1.,0.0, 0.5,1.]
#nodes_y = [0.0, 0.0, 0.0 ,0.5,0.5,0.5,1.,1.,1.]
#nodal_values = [1.0,1.0,1.0,1.0,1.0,1.0,2.0,1.0,1.0]
#elements_tris = []
#elements_quads = [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]]
#elements = elements_tris + elements_quads
#        
##nodes_x = a[:,1]
##nodes_y = a[:,2]
#
#
## convert all elements into triangles
#elements_all_tris = elements_tris + quads_to_tris(elements_quads)
#
## create an unstructured triangular grid instance
#triangulation = tri.Triangulation(nodes_x, nodes_y, elements_all_tris)
#
## plot the finite element mesh
#plot_fem_mesh(nodes_x, nodes_y, elements)
#
## plot the contours
#plt.tricontourf(triangulation, nodal_values)
#
## show
#plt.colorbar()
#plt.axis('equal')
#plt.show()


