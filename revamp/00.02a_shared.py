import re
fin = open("shared-bed-files.txt",'r')
fout = open("shared_processed.txt",'w')
count = 0

for entry in fin:
	f = entry.strip()
	print f
	fi = open(f,'r')
	for line in fi:
		l = line.strip().split()
		if len(l) ==4:
			fout.write("\t".join(l)+"\n")
		else:
			count += 1
	fi.close()
print count
fin.close()
fout.close()
