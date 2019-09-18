
import json

with open('result_shells.json') as handle:
	shells = json.load(handle)

top_shell = set(shells[str(len(shells))])
flag = False

with open('final.unitigs.fa','r') as f, open('top_unitigs.fa','w+') as outf:
	i = 0
	for line in f:
		if i%2 == 0:
			token  = line.split(' ')[0].split('>')[1]
			if token in top_shell:
				print('YES')
				flag = True
				outf.write(line)
				outf.write('\n')
			i+=1
		else:
			if flag:
				outf.write(line)
				outf.write('\n')
				flag = False
			i+=1

