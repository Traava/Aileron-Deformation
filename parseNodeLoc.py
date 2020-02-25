import numpy as np

a = np.genfromtxt('B737.inp', delimiter = ',', skip_header = 9, skip_footer = 14594-6598)






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

        
    
