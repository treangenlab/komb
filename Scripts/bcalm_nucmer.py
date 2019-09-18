import matplotlib.pyplot as plt
from sklearn.metrics import auc as ska
import numpy as np
import math
import time
from collections import defaultdict
import json

true_unitigs = defaultdict(int)
true_neg = []
shell_dict = defaultdict(list)

with open('result_shells.json') as handle:
    dict_shells = json.load(handle)
for k in dict_shells:
    shell_dict[int(k)] = dict_shells[k]

print ('Loaded shells')
num_shells = len(shell_dict)

with open('/Users/advaitbalaji/Downloads/repeats.txt','r') as repfile:
    for line in repfile:
        unitig = int(line.split('\t')[-1].split('\n')[0])
        identity = float(line.split('\t')[-3])
        if identity > 99:
            true_unitigs[unitig]+=1
            
for key in list(true_unitigs):
    if true_unitigs[key] < 2:
        true_neg.append(key)
        del true_unitigs[key]
        
all_unitigs = []
with open('/Users/advaitbalaji/Downloads/final.unitigs.fa','r') as allfile:
    i = 0
    for line in allfile:
        if i%2 == 0:
            all_unitigs.append(int(line.split(' ')[0].split('>')[1]))
        i+=1

for i in all_unitigs:
    if not i in true_unitigs:
        true_neg.append(i)

print("Number of shells: "+str(len(shell_dict)))

top = [i for i in range(1,len(shell_dict)-13)]

sens = []
prec =[]
for val in top:
    shell = len(shell_dict) - val + 1
    #pred_pos = dict_ft[shell]
    pred_pos = []
    for i in xrange(shell,len(shell_dict)+1):
        pred_pos = pred_pos + shell_dict[i]
    pred_neg = [unitig for unitig in all_unitigs if not unitig in pred_pos]

    tp = float(len([unitig for unitig in pred_pos if unitig in true_unitigs]))
    fp = float(len([unitig for unitig in pred_pos if not unitig in true_unitigs]))
    fn = float(len([unitig for unitig in true_unitigs if not unitig in pred_pos]))
    true_n = float(len([unitig for unitig in pred_neg if unitig in true_neg]))
    
    #print(fp)
    
    sensitivity = float(tp/(tp+fn))
    precision = float(tp/(tp+fp))
    sens.append(sensitivity)
    prec.append(precision)

print(sens)
print(prec)

top_mod = [val-0.3 for val in top]
plt.bar(top_mod,sens,width = 0.3, color = 'b', align='center', label = 'Sensitivity')
plt.bar(top,prec,width = 0.3, color = 'r', align='center', label = 'Precision')
lgd = plt.legend(bbox_to_anchor=(1.01, 1), loc='best')
plt.xticks(np.arange(min(top), max(top)+1, 5))
plt.title('Validation: Francisella')
plt.xlabel('Top n shells considered')
plt.ylabel('Value')
plt.savefig("/Users/advaitbalaji/Desktop/Validation_francisella.png", bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 200, format='png')
plt.show()

    