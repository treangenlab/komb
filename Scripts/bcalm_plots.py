from collections import defaultdict
import time
import matplotlib.pyplot as plt
import numpy as np
import json
with open('/Users/advaitbalaji/Downloads/result_shells.json') as handle:
    temp = json.load(handle)
num_nodes = []
for i in xrange(1,len(temp)+1):
    num_nodes.append(len(temp[str(i)]))
print(num_nodes)
x = [i+1 for i in range(len(num_nodes))]
plt.bar(x,num_nodes,align='center') # A bar chart
#plt.xticks(np.arange(min(x), max(x)+1, 3))
plt.xticks(np.arange(min(x), max(x)+1, 15))
plt.yticks(np.arange(min(num_nodes), max(num_nodes), 50))
plt.title('Nodes per shell: Haloferax volcanii DS2')
plt.xlabel('Shell')
plt.ylabel('Number of nodes')
plt.savefig("/Users/advaitbalaji/Desktop/Yersinia3", dpi = 200, format='png')
plt.show()
