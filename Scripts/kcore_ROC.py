import matplotlib.pyplot as plt
from sklearn.metrics import auc as ska
import numpy as np
import math
import time
from collections import defaultdict
import json
import sys

tpr = []
fpr = []

shell_dict = defaultdict(list)

with open('result_shells.json') as handle:
    dict_shells = json.load(handle)
for k in dict_shells:
    shell_dict[int(k)] = dict_shells[k]

print ('Loaded shells')
num_shells = len(shell_dict)

all_unitigs = defaultdict(int)
true_unitigs = defaultdict(int)
true_neg = []
all_unitigs = []

with open('/Users/advaitbalaji/Downloads/repeats.txt','r') as repfile:
    for line in repfile:
        unitig = int(line.split('\t')[-1].split('\n')[0])
        identity = float(line.split('\t')[-3])
        if identity == 100:
            true_unitigs[unitig]+=1
            
for key in list(true_unitigs):
    if true_unitigs[key] < 2:
        true_neg.append(key)
        del true_unitigs[key]

print("Read Truth Set")

with open('/Users/advaitbalaji/Downloads/final.unitigs.fa','r') as allfile:
    i = 0
    for line in allfile:
        if i%2 == 0:
            all_unitigs.append(int(line.split(' ')[0].split('>')[1]))
        i+=1

for i in all_unitigs:
    if not i in true_unitigs:
        true_neg.append(i)
        
print("Read untig file")
        
top = [i for i in range(1,len(shell_dict)-13)]

sens = []
prec =[]
for val in top:
    shell = len(shell_dict) - val + 1
    #pred_pos = dict_ft[shell]
    pred_pos = []
    for i in xrange(shell,len(shell_dict)+1):
        pred_pos = pred_pos+shell_dict[i]
    
    pred_neg = [unitig for unitig in all_unitigs if not unitig in pred_pos]
    
    tp = float(len([unitig for unitig in pred_pos if unitig in true_unitigs.keys()]))
    fp = float(len([unitig for unitig in pred_pos if not unitig in true_unitigs.keys()]))
    fn = float(len([unitig for unitig in true_unitigs.keys() if not unitig in pred_pos]))
    true_n = float(len([unitig for unitig in pred_neg if unitig in true_neg]))
    

    sensitivity = float(tp/(tp+fn))
    specificity = float(true_n/(true_n+fp))
    false_pr = 1 - specificity
    tpr.append(sensitivity)
    fpr.append(specificity)
    #print("Sensitivity or TPR: "+str(sensitivity))
    #print("Specificity: " + str(specificity))
    #print("FPR: "+str(false_pr))
    print("Cutoff: "+str(shell))
    

print(tpr)
print(fpr)
auc = ska(fpr,tpr)
fpr = [0]+fpr+[1]
tpr = [0]+tpr+[1]
print(auc)

data_n = "RB_400_50x_100pe"
plt.figure()
plt.plot(fpr,tpr,color='darkorange',lw=2, linestyle='-', label='ROC curve (area = %0.2f)' % auc)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic for '+data_n)
plt.legend(loc="upper left")
plt.savefig("/Users/advaitbalaji/Desktop/ROC_"+data_n+".png", dpi = 200, format='png')
plt.show()