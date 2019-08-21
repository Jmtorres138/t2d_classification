#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for running atac-seq peak annotation enrichments
Usage:
module add bedtools
python script fname out_file
'''

import os,sys,gzip
import subprocess as sp
import numpy

work_dir = "/well/mccarthy/users/jason/projects/t2d_classification/revamp/"
genome_file = work_dir + "hg19.chrom.sizes"
annot_dir = "/well/mccarthy/users/jason/datasets/Roadmap/H3K27ac/gappedPeak/"


def get_intersect(snp_bed,annot_bed,temp_name,proc_dir):
	command = ["/apps/well/bedtools/2.24.0/bedtools", "intersect", "-wa","-a",snp_bed,"-b",annot_bed,"|","uniq",">",proc_dir+"inter.temp."+temp_name+".bed"]
	sp.check_call(" ".join(command),shell=True)
	command = ["cat", proc_dir+"inter.temp."+temp_name+".bed","|","wc -l"]
	p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
	output, err = p.communicate(b"input data that is passed to subprocess' stdin")
	count = output.strip()
	os.remove(proc_dir+"inter.temp."+temp_name+".bed")
	return count

def enrich(fname,snp_bed,proc_dir):
	temp_name = "enrichment-tests"
	fin = gzip.open(fname,'rb')
	#fin.readline() # use if there is a header
	fout = open(proc_dir+"temp1."+temp_name+".bed",'w')
	for line in fin:
		l = line.strip().split()
		chrom,start,end,name =l[0],l[1],l[2],l[3]
		fout.write("\t".join([chrom,start,end,name])+"\n")
	fin.close()
	fout.close()
	command1 = ["/apps/well/bedtools/2.24.0/bedtools", "sort -i", proc_dir+"temp1."+temp_name+".bed", ">", proc_dir+"temp2."+temp_name+".bed"]
	command2 = ["/apps/well/bedtools/2.24.0/bedtools", "merge -i", proc_dir+"temp2."+temp_name+".bed",">", proc_dir+"temp."+temp_name+".bed"]
	sp.check_call(" ".join(command1),shell=True)
	sp.check_call(" ".join(command2),shell=True)
	os.remove(proc_dir+"temp1."+temp_name+".bed")
	os.remove(proc_dir+"temp2."+temp_name+".bed")
	observed = int(get_intersect(snp_bed,proc_dir+"temp."+temp_name+".bed",temp_name,proc_dir))
	iter_list = []
	for i in range(1,1001):
		sys.stdout.write("\r%d"%i)
		sys.stdout.flush()
		command = ["/apps/well/bedtools/2.24.0/bedtools", "shuffle -i", proc_dir+"temp."+temp_name+".bed", "-g",genome_file, "-chrom", ">", proc_dir+"shuffle.temp."+temp_name+".bed"]
		sp.check_call(" ".join(command),shell=True)
		rand_count = get_intersect(snp_bed,proc_dir+"shuffle.temp."+temp_name+".bed",temp_name,proc_dir)
		iter_list.append(int(rand_count))
		os.remove(proc_dir+"shuffle.temp."+temp_name+".bed")
	print ("\nObserved: %d" % observed)

	pval =  (sum([x>=observed for x in iter_list])+1) / float(1000+1)
	print ("Pvalue: %f" % pval)
	return [observed,numpy.mean(iter_list), observed/numpy.mean(iter_list), pval]

def run_enrichments(file_list,out_file,proc_dir,shared=False):
	afile_list = sorted([f for f in os.listdir(annot_dir) if ".gappedPeak.gz" in f])
	aname_list = [s.split("-")[0] for s in afile_list]
	print afile_list

	fout = open(out_file,'w')
	fout.write("\t".join(["tissue_toa","tissue_H3K27ac","observed.count","null.mean","enrich.factor","p.val"])+"\n")
	for f in file_list:
		tname = f.split("_")[1].split(".")[0]
		print tname
		for i in range(0,len(aname_list)):
			aname = aname_list[i]
			a = annot_dir + afile_list[i]
			print(aname)
			out_list = enrich(a,proc_dir+f,proc_dir)
			out_list = [str(x) for x in out_list]
			write_list = [tname,aname] + out_list
			fout.write("\t".join(write_list)+"\n")
	fout.close()

def main():
	#query_dir = work_dir + "enrichment_files/enrich_H3K27ac/thresh00/"
	#file_list = sorted(os.listdir(query_dir))
	#l00 = [x for x in file_list if "thresh00" in x]
	#run_enrichments(l00,query_dir+"enrich_H3K27ac_t00.txt",query_dir)

	#query_dir = work_dir + "enrichment_files/enrich_H3K27ac/thresh20/"
	#file_list = sorted(os.listdir(query_dir))
	#l20 = [x for x in file_list if "thresh20" in x]
	#run_enrichments(l20,query_dir+"enrich_H3K27ac_t20.txt",query_dir)

	#query_dir = work_dir + "enrichment_files/enrich_H3K27ac/thresh50/"
	#file_list = sorted(os.listdir(query_dir))
	#l50 = [x for x in file_list if "thresh50" in x]
	#run_enrichments(l50,query_dir+"enrich_H3K27ac_t50.txt",query_dir)

	query_dir = work_dir + "enrichment_files/enrich_H3K27ac/thresh80/"
	file_list = sorted(os.listdir(query_dir))
	l80 = [x for x in file_list if "thresh80" in x]
	run_enrichments(l80,query_dir+"enrich_H3K27ac_t80.txt",query_dir)




if (__name__=="__main__"): main()
