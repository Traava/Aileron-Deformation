import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.spatial import distance
#import MatrixSystem_validation as ms

###-Parsing node locations and elements
node = np.genfromtxt('B737.inp', delimiter = ',', skip_header = 9, skip_footer = 14594-6598)
elem = np.genfromtxt('B737.inp', delimiter = ',', skip_header = 6598, skip_footer = 14594-13233)

####---------Parsing stresses & deflections---------#################


f = open('B737.rpt')

rawStress = [[],[],[],[],[],[] ] 
rawDef = [[] for _ in range(3)]
pos=0
node_label_count = 0
sizes = [5778, 856,5778, 856,5778, 856]
sizes2 = 6588
for line in f:
    if line.startswith('   Element Label'):
        next(f)
        next(f)
        while len(rawStress[pos]) < sizes[pos]:
            rawStress[pos].append(np.array(next(f).split()))
        pos +=1
    if line.startswith('      Node Label'):
        node_label_count +=1
        if (node_label_count % 2 == 1) and (node_label_count // 2 <= 2):
            next(f)
            next(f)
            while len(rawDef[node_label_count//2]) < sizes2:
                rawDef[node_label_count//2].append(np.array(next(f).split()))
                

        
rawStress = np.array(rawStress)
rawDef = np.array(rawDef)


###################################################
###----organising data into dictionaries-----

stressdata = []

for j in range(6):
    elem_vm_dict = {}
    elem_s_dict = {}
    for i in rawStress[j]:
        elem_vm_dict[int(i[0])] = sum([float(i[2]),float(i[3])])/2
        elem_s_dict[int(i[0])] = sum([float(i[4]),float(i[5])])/2
    stressdata.append(elem_vm_dict)
    stressdata.append(elem_s_dict)
    
defdata = []

for j in range(3):
    node_def = {}
    for i in rawDef[j]:
        node_def[int(i[0])] = np.array([float(i[1]), float(i[2]), float(i[3]), float(i[4])])
    defdata.append(node_def)


#####################
    
nodes = {}
elems = {}

for i in range(len(node)):
    nodes[node[i,0]] = node[i,1:]
for i in range(len(elem)):
    elems[elem[i,0]] = elem[i,1:]
###-------------Splitting elements into sections, and finding the max stress at each section-------------

def deflections(def_data):
    hinge_nodes= []
    
    xpos = []
    ydef = []
    zdef = []
    
    for nod in nodes:
        if nodes[nod][1] == 0 and nodes[nod][2] == 0:
            hinge_nodes.append(nod)
    
    hinge_nodes = np.array(hinge_nodes)
    
    for nod in hinge_nodes:
        xpos.append(nodes[nod][0])
        ydef.append(def_data[nod][2])
        zdef.append(def_data[nod][3])
        
    xpos = np.array(xpos)/1000
    ydef = -np.array(ydef)/1000# converting to meters and swapping signs to fit with our coord system
    zdef = -np.array(zdef)/1000
    
        
    for nod in hinge_nodes:
        print(nodes[nod])
        print(def_data[nod][2:4])
        print()
        
    print('max ydef: ', max(ydef), ' max zdef: ', max(zdef))
    return xpos, ydef, zdef


    
def maxSectionStress(sections, stress_data ):
    #sections = 100
    sec_bounds = np.linspace(0, 2661, sections+1)
    elem_sections = [{} for _ in range(sections)]
    maxstr_sections = []
    
    elem_dict = {}

    
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
    
        for j in range(len(sec_bounds)-1):
            if elem_dict[i][0] > sec_bounds[j] and elem_dict[i][0] < sec_bounds[j+1]:
                elem_sections[j][i] = elem_dict[i]
        
    
    for sec in elem_sections:
        stress = []
        for el in sec:
            stress.append(stress_data.get(el, 0))
        maxstr = max(abs(np.array(stress)))
        maxstr_sections.append(maxstr)
        
    sec_mid = [sec_bounds[i]+sec_bounds[i+1]/2 for i in range(sections)]
    
    return sec_mid, maxstr_sections




cases = ['von Mises, Bending (reg. 1)','Shear, Bending (reg. 1)', 'von Mises, Bending (reg. 2)','Shear, Bending (reg. 2)', 'Von Mises, Jam-Bent (reg. 1)','Shear, Jam-Bent (reg. 1)','Von Mises, Jam-Bent (reg. 2)','Shear, Jam-Bent (reg. 2)', 'Von Mises, Jam-Straight (reg. 1)','Shear, Jam-Straight (reg. 1)', 'Von Mises, Jam-Straight (reg. 2)','Shear, Jam-Straight (reg. 2)']                                                 

#for i in range(len(stressdata)):
#    sec_mid, maxstr_sections = maxSectionStress(100, stressdata[i])
#    plt.plot(sec_mid, maxstr_sections, label = cases[i])    
#plt.legend()
#plt.show()

#----Parsing txt files of our deflection results -----
deflect = [[] for _ in range(6)]

xx = np.linspace(0,2.661, 100)


for i in range(3):
    f= open('Case_'+str(i+1)+'.txt', 'r')
    j =1
    for line in f:
        deflect[2*i].append(float(line.split()[1]))
        deflect[2*i+1].append(float(line.split()[2]))


        
deflection_labels  = []
plt.figure()
j = 0
for i in range(len(defdata)):
    xpos, ydef, zdef = deflections(defdata[i])
    plt.subplot(int('32'+str(i+j+1)))
    plt.scatter(xpos, ydef, label = 'Validation model', marker ='+', color='#ff7f0e')
    plt.plot(xx, deflect[j+i], label = 'Numerical model')
    plt.xlabel("Spanwise position [m]").set_size(14)
    plt.ylabel("Y-deflection [m]").set_size(14)
    plt.legend()
    plt.grid()


    j+=1##z:
    plt.subplot(int('32'+str(i+j+1)))
    plt.scatter(xpos, zdef, label = 'Validation model', marker ='+', color='#ff7f0e')
    plt.plot(xx, deflect[j+i], label = 'Numerical model')
    plt.xlabel("Spanwise position [m]").set_size(14)
    plt.ylabel("Z-deflection [m]").set_size(14)

    plt.legend()
    plt.grid()


    
plt.legend()
plt.show()
