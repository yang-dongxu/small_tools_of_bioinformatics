#! /usr/bin/env python3

import os
import sys
import csv

USAGE = '{} <proj_file.tsv> <output_file.md5>'.format(os.path.basename(sys.argv[0]))

proj_file = sys.argv[1]
output_file = sys.argv[2]

def md5sum_c_files(rows):
    for row in rows:
        for md5, fastq_file in zip(row.get('fastq_md5', None).split(';'), row.get('fastq_ftp', None).split(';')) :
            fastq_file = fastq_file.rsplit('/',1)[1]
            yield f'{md5}  {fastq_file}'


with open(proj_file) as input_fhd, \
     open(output_file, 'w') as output_fhd:
    intput_csv = csv.DictReader(input_fhd, delimiter='\t')
    output_fhd.write('\n'.join(md5sum_c_files(intput_csv))+'\n')
