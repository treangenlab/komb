import glob
from collections import defaultdict
import json
import numpy as np
import matplotlib.pyplot as plt

def makeTotalNode():
	res  = glob.glob('SORT_*')
	total_nodes = defaultdict(list)

	for res_dir in res:
		with open(res_dir+'/result_shells.json') as handle:
			temp_dict = json.load(handle)
			for k in temp_dict:
				total_nodes[int(k)] = total_nodes[int(k)] + temp_dict[k]
		print('Finished '+ res_dir)

	with open('total_nodes.json', 'w') as fp:
		json.dump(total_nodes, fp)


with open('result_shells.json') as handle:
	shells_sg = json.load(handle)

with open('total_nodes.json') as handle:
	shells_mg = json.load(handle)

node_count_sg = []
node_count_mg = []

for i in xrange(1,301):
	node_count_sg.append(len(shells_sg[str(i)]))
	node_count_mg.append(len(shells_mg[str(i)]))

print node_count_sg
print node_count_mg

x = [i+1 for i in xrange(300)]

plt.bar(x,node_count_sg,width = 0.8, color = 'b', alpha = 0.7, align= 'center', label = 'Single genome mode')
plt.bar(x,node_count_mg, width =0.6, color = 'r', alpha = 0.9, align = 'center', label = 'Metagenome mode')
plt.axvline(x = 1.33, color='black', lw=1, linestyle='--', label = 'Average genomes per genus')
lgd = plt.legend(bbox_to_anchor=(1.01, 1), loc='best')
plt.yticks(np.arange(0, 50001, 5000))
plt.xticks(np.arange(min(x), max(x)+1, 20))
plt.title('Mode comparison')
plt.xlabel('Shell')
plt.ylabel('Number of nodes')
plt.savefig('Shakya.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 200, format='png')
plt.show()

#Metagenome genus level and family level
#Just keep the yersinia plots
#Try do a stackled bar plot
#Kraken krona plot 
#Real metagenome dataset
#Show whether kraken recovers everything - Table 
#


