import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.spatial import distance

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

elemxz = [] #remove values where y<0

for i in elems:
    
    coords = []
    nod = elems[i]
    for n in nod:
        coords.append(nodes[n])
    coords = np.array(coords)
    x = sum(coords[:,0])/4
    y = sum(coords[:,1])/4
    z = sum(coords[:,2])/4
    

    
    elem_dict[i] = np.array([x,y,z])
    elem_lst.append(np.array([x,y,z]))
    if y>0:
        elemxz.append(np.array([x,z]))

elemxz = np.array(elemxz)
elem_lst = np.array(elem_lst)



def closest_node(node, nodes):
    elem_num = distance.cdist([node], nodes).argmin() +1
    return elem_num

####---------Parsing stresses---------#################

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
        
data = np.array(data)

elem_vm_dict = {}
elem_vm = []

for i in data[0]:
    elem_vm_dict[float(i[0])] = sum([float(i[2]),float(i[3])])/2
    elem_vm.append(sum([float(i[2]),float(i[3])])/2)
    


###------------Pygame----------------
    
def color(val):# zero to 0.362
    R = round(255 -val*255/0.362)
    B = round(val*255/0.362)
    return (R,0,B)


    
import pygame

pygame.init()
screen = pygame.display.set_mode((2661, 605))
done = False

#while not done:
for i in elems:
    verts = []
    for j in elems[i]:
        verts.append((nodes[j][0],nodes[j][2]+502.5))
    
    col = color(elem_vm_dict.get(int(i-1), 0))
    #TODO: elem_vm has less elements than elems. Probably need to use elem_vm_dict so we know what elem we are using.
    print(col)
    pygame.draw.polygon(screen,col, verts)


for event in pygame.event.get():
        if event.type == pygame.QUIT:
                done = True
                


pygame.display.flip()


#########---------PLOTTING------------------
#grain = 200
#
#x = np.linspace(10, 2650, grain)
#z = np.linspace(-490, 101, grain)
#y = np.array([elem_vm_dict.get(closest_node([i,j],elemxz), 0.5) for j in z for i in x])
#
#
#X, Z = np.meshgrid(x, z)
#Y = y.reshape(grain, grain)
#
#plt.pcolor(X, Z, Y)
##plt.imshow(Y,origin='lower',interpolation='bilinear')
#
#
#plt.show()





##----TRYING TO FIND ALL ELEMS ACCOCIATED WITH EACH NODE---
#nodes_elem = np.zeros((len(nodes),4))
#for n in range(len(nodes)):
#    els = []
#    
#    for e in elem:
#
#        if n+1 in e[1:]:
#            print(n+1, ' is in ',e[1:] )
#            els.append(e[0])
#    print(els)
#    nodes_elem[n] = np.array(els)
#    
#    


















#
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
#elements_quads = [[0, 1, 4,3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]]
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


