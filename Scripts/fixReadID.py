import sys

'''
Usage: python fixReadID.py [READ FILE TO FIX]

'''


i = 0;
with open(sys.argv[1],"r") as f, open(sys.argv[1].split(".fq")[0]+"_cor"+".fq","w+") as outf:
	for line in f:
		if(i%4 == 0):
			modline = '.'.join(line.split(' ')[0].split('.')[:-1]) + ' '+' '.join(line.split(' ')[1:])
			outf.write(modline)
			i+=1
		else:
			outf.write(line)
			i+=1