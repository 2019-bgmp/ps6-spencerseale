#!/usr/bin/env python3

import re

FILE = "contigs.fa"

# creating an empty list to store fasta header id
fa_id = []
with open(FILE, "r") as opfile:
    for line in opfile:
        line = line.strip("\n")
        if ">" in line:
            fa_id.append(line)
#print(fa_id)

# creating empty lists to store kmer length and coverage
kmer_length = []
kmer_coverage = []
for i in fa_id:
    kmer_len = re.findall("h_[0-9]+", i)
    kmer_length.append(kmer_len[0].strip("h_"))
    kmer_cov = re.findall("v_[0-9]+.[0-9]+", i)
    kmer_coverage.append(kmer_cov[0].strip("v_"))
#print(kmer_coverage)
#print(kmer_length)

# ps6 states to set k to 49
k = 49
# calculating physical length and changing kmer_length array to physical length as ps stated
for idx, kcnt in enumerate(kmer_length):
    kmer_length[idx] = int(kcnt) + k - 1
# print(kmer_length) # used to determine if the list was correctly updated

# calculating number of contigs using header list fa_id
num_contig = len(kmer_length)
print("# of contigs:", num_contig)

# calculating max contig length
max_contig_len = 0
for c in kmer_length:
    if c > max_contig_len:
        max_contig_len = c
print("max contig length:",max_contig_len)

# calculating mean contig length
mean_contig_len = (sum(kmer_length)) / num_contig
print("mean contig length:",mean_contig_len)

# calculating genome length
gen_len = 0
for g in kmer_length:
    gen_len += g
print("genome length:",gen_len)

# calculating mean kmer depth of coverage
kmer_depth = 0
for i in kmer_coverage:
    kmer_depth += float(i)
mean_kmer_depth = kmer_depth / num_contig
print("This is the mean kmer depth of coverage:", mean_kmer_depth)

# calculating mean depth of coverage
depth_of_coverage = []
count = 0
for depth_cov in kmer_coverage:
    ck = float(depth_cov)
    l = kmer_length[count]
    c = (ck * l) / (l - k + 1)
    depth_of_coverage.append(c)
    count += 1

# calculating mean depth of coverage
mean_depth_coverage = (sum(depth_of_coverage)/num_contig)
print("This is the mean depth of coverage:", mean_depth_coverage)

# determining N50
kmer_length.sort()
#print(kmer_length) # used to check if the kmer length was actually sorted

tot_sum = 0
for z in range(len(kmer_length)):
    tot_sum += kmer_length[z]
check = tot_sum / 2
#print(check)
#print(tot_sum)

# determining which contig contains the N50
idx_counter = 0
list_sum_counter = 0
list_sum_counter += (kmer_length[idx_counter])
while check >= list_sum_counter:
    idx_counter += 1
    list_sum_counter += (kmer_length[idx_counter])
#print(idx_counter)
N50 = kmer_length[idx_counter]
print("N50:", N50)

# calculating distribution of contig lengths and adding to buckets
# taking each value in list, dividing by 100 and truncating at decimal, then mult by 100
# then adding that as key and incrementing value
contig_dist = {}
for value in kmer_length:
    calc = (value // 100) * 100
    if calc in contig_dist:
        contig_dist[calc] += 1
    else:
        contig_dist[calc] = 1
# print(contig_dist)

# turning into list and printing out
print("# Contig length\tNumber of contigs in this category")
for i, v in sorted(contig_dist.items()):
    print(i,"\t",v)
