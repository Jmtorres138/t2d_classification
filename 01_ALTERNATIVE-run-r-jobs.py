#!/usr/bin/python -O
# Jason Matthew Torres
'''

'''
# libraries
import sys,os,gzip
import subprocess as sp


# globals

cur_dir = "/well/mccarthy/users/jason/projects/t2d_classification/"
rscript = "/apps/well/R/3.3.1/bin/Rscript"
log_dir = cur_dir + "logs/"
job_dir = cur_dir + "jobs/"

num_genes = 56202

def get_sets(num_genes):
    set_list = []
    count=1
    while count<num_genes:
        count +=100
        if count > num_genes:
            count=num_genes
        my_end = count
        my_start = count - 100
        set_list.append([my_start,my_end])
    return set_list




def run_job(my_start,my_end):

    command = rscript + " --vanilla " + cur_dir + "01_expression-specificity-scores.R " + str(my_start) + " " + str(my_end)

    job_file = job_dir+"job_indices_"+str(my_start)+"-"+str(my_end)+".sh"
    fout=open(job_file,'w')
    script='''
#$ -N indices_%s_%s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %sindices_%s_%s.error
#$ -o %sindices_%s_%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % (str(my_start),str(my_end),
               log_dir, str(my_start),str(my_end),
               log_dir, str(my_start),str(my_end),
               command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)



def main():
    set_list = get_sets(num_genes)
    for l in set_list:
        print l
        run_job(l[0],l[1])

if (__name__=="__main__"): main()
