#!/usr/bin/env python

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=ps6-p2_gen      ### Job Name
#SBATCH --output=ps6.out         ### File in which to store job output
#SBATCH --error=ps6.err          ### File in which to store job error messages
#SBATCH --time=0-01:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp      ### Account used for job submission

DIR = "/projects/bgmp/shared/Bi621/"
IN_FILE_PE1 = "/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_1"
IN_FILE_PE2 = "/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq_2"
IN_FILE_SR = "/projects/bgmp/shared/Bi621/800_3_PE5_interleaved.fq.unmatched"

# given fosmid library specifications
fosmid_len = 40000
num_fosmids = 50
tot_nt = fosmid_len * num_fosmids
#print(tot_nt)

# calculating number nt sequenced in all three files
seq_len_PE1 = []
with open(IN_FILE_PE1, "r") as PE1file:
    LN1 = 0
    for line in PE1file:
        LN1 += 1
        line = line.strip("\n")
        if LN1 % 4 == 2:
            #print(line)
            seq_len_PE1.append(len(line))

seq_len_PE2 = []
with open(IN_FILE_PE2, "r") as PE2file:
    LN2 = 0
    for line in PE2file:
        LN2 += 1
        line = line.strip("\n")
        if LN2 % 4 == 2:
            #print(line)
            seq_len_PE2.append(len(line))

seq_len_SR = []
with open(IN_FILE_SR, "r") as SRfile:
    LN3 = 0
    for line in SRfile:
        LN3 += 1
        line = line.strip("\n")
        if LN3 % 4 == 2:
            #print(line)
            seq_len_SR.append(len(line))

#summing total number nt seq for each file
num_nt_sq_PE1 = sum(seq_len_PE1)
num_nt_sq_PE2 = sum(seq_len_PE2)
num_nt_sq_SR = sum(seq_len_SR)

print(num_nt_sq_PE1)
print(num_nt_sq_PE2)
print(num_nt_sq_SR)

# summing number of reads
cnt = (LN1 + LN2 + LN3) / 4
print("cnt:", cnt)
# summing total num of nt seq
tot_nt_seq = num_nt_sq_PE1 + num_nt_sq_PE2 + num_nt_sq_SR
# calc depth of coverage
tot_depth_cov = float(tot_nt_seq) / float(tot_nt)
# calc Lmean
avg_seq_len = ((float(tot_nt_seq)) / cnt)


print("avg_seq_len:", avg_seq_len)

print("total cov depth:", tot_depth_cov)
print("total nt seq:", tot_nt_seq)
print("avg_seq_len:", avg_seq_len)

# calc kmer covg for velvetg input
kmer_cov_31 = (tot_depth_cov * (avg_seq_len - 31 + 1)) / avg_seq_len
kmer_cov_41 = (tot_depth_cov * (avg_seq_len - 41 + 1)) / avg_seq_len
kmer_cov_49 = (tot_depth_cov * (avg_seq_len - 49 + 1)) / avg_seq_len

print("kmer31:", kmer_cov_31)
print("kmer41:", kmer_cov_41)
print("kmer49:", kmer_cov_49)
