import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.spatial import distance


###-Parsing node locations and elements
node = np.genfromtxt('B737.inp', delimiter = ',', skip_header = 9, skip_footer = 14594-6598)
elem = np.genfromtxt('B737.inp', delimiter = ',', skip_header = 6598, skip_footer = 14594-13233)

####---------Parsing stresses---------#################


f = open('B737.rpt')

rawdata = [[],[],[],[],[],[] ]
pos=0
sizes = [5778, 856,5778, 856,5778, 856]
for line in f:
    if line.startswith('   Element Label'):
        next(f)
        next(f)
        while len(rawdata[pos]) < sizes[pos]:
            rawdata[pos].append(np.array(next(f).split()))
        pos +=1
        
rawdata = np.array(rawdata)

###################################################
###----

data = []

for j in range(6):
    elem_vm_dict = {}
    elem_s_dict = {}
    for i in rawdata[j]:
        elem_vm_dict[float(i[0])] = sum([float(i[2]),float(i[3])])/2
        elem_s_dict[float(i[0])] = sum([float(i[4]),float(i[5])])/2
    data.append(elem_vm_dict)
    data.append(elem_s_dict)
    


#####################
###-------------Splitting elements into sections, and finding the max stress at each section-------------
    
    
def maxSectionStress(sections, stress_data ):
    #sections = 100
    sec_bounds = np.linspace(0, 2661, sections+1)
    elem_sections = [{} for _ in range(sections)]
    maxstr_sections = []
    
    nodes = {}
    elems = {}
    elem_dict = {}
    elem_lst = []
    
    for i in range(len(node)):
        nodes[node[i,0]] = node[i,1:]
    for i in range(len(elem)):
        elems[elem[i,0]] = elem[i,1:]
    
    
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
    
        for j in range(len(sec_bounds)-1):
            if elem_dict[i][0] > sec_bounds[j] and elem_dict[i][0] < sec_bounds[j+1]:
                elem_sections[j][i] = elem_dict[i]
        
    elem_lst = np.array(elem_lst)
    
    for sec in elem_sections:
        stress = []
        for el in sec:
            stress.append(stress_data.get(el, 0))
        maxstr = max(abs(np.array(stress)))
        maxstr_sections.append(maxstr)
        
    sec_mid = [sec_bounds[i]+sec_bounds[i+1]/2 for i in range(sections)]
    
    return sec_mid, maxstr_sections


cases = ['von Mises, Bending (reg. 1)','Shear, Bending (reg. 1)', 'von Mises, Bending (reg. 2)','Shear, Bending (reg. 2)', 'Von Mises, Jam-Bent (reg. 1)','Shear, Jam-Bent (reg. 1)','Von Mises, Jam-Bent (reg. 2)','Shear, Jam-Bent (reg. 2)', 'Von Mises, Jam-Straight (reg. 1)','Shear, Jam-Straight (reg. 1)', 'Von Mises, Jam-Straight (reg. 2)','Shear, Jam-Straight (reg. 2)']                                                 

for i in range(len(data)):
    sec_mid, maxstr_sections = maxSectionStress(100, data[i])
    plt.plot(sec_mid, maxstr_sections, label = cases[i])    
plt.legend()
plt.show()
    
