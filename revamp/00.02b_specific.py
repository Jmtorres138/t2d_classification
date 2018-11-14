import re
# Need to run this first to create file
# ls -d -1 $PWD/*.* | grep specific > specific-bed-files.txt
fin = open("/well/mccarthy/users/jason/projects/t2d_classification/revamp/analysis_files/specific-bed-files.txt",'r')
fout = open("/well/mccarthy/users/jason/projects/t2d_classification/revamp/analysis_files/specific_processed.txt",'w')
count = 0

for entry in fin:
	f = entry.strip()
	print f
	fi = open(f,'r')
	for line in fi:
		l = line.strip().split()
		if len(l) == 4:
			fout.write("\t".join(l)+"\n")
		else:
			count += 1
	fi.close()
print count
fin.close()
fout.close()
