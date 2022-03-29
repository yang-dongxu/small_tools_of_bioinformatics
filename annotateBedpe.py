#!/usr/env/bin python
import os
import sys
import io
import subprocess
import click
import logging
from tempfile import NamedTemporaryFile
import pandas as pd

def split_bedpe(bedpe):
    bed1 = NamedTemporaryFile("w",delete=False)
    bed2 = NamedTemporaryFile("w",delete=False)

    bed1_name = bed1.name
    bed2_name = bed2.name

    bed1.close()
    bed2.close()

    cmd = '''cat %s | cut -f 1-3 | awk 'BEGIN{OFS="\t"}{$(NF+1) =NR ; print $0}' > %s   '''
    cmd1 = cmd%(bedpe,bed1_name)
    os.system(cmd1)

    cmd = '''cat %s | cut -f 4-6 | awk 'BEGIN{OFS="\t"}{$(NF+1) =NR ; print $0}' > %s   '''
    cmd2 = cmd%(bedpe,bed2_name)
    os.system(cmd2)
    return bed1_name, bed2_name

def combine_bed_to_bedpe(bed1,bed2,feature_num,oname):
    cmd = '''join -j 4 %s %s | \
        awk -v f=%d 'BEGIN{OFS="\t"}
        {printf("%%s\t%%d\t%%d\t%%s\t%%d\t%%d",$2,$3,$4,$(5+f),$(6+f),$(7+f)); 
        for(i=1;i<=f;i++){
            r1=$(4+i);
            r2=$(7+f+i);
            if (r1 >0 && r2>0){o=3;}
            else if (r1 >0 && r2 == 0){o=1;}
            else if (r1 ==0 && r2>0){o=2;}
            else if (r1==0 && r2==0){o=0;}
            printf("\t%%d",o);
        };
        printf("\\n");
        }
        '  '''%(bed1,bed2,feature_num)
    subprocess.run(cmd,shell=True,stdout=oname)
    os.remove(bed1)
    os.remove(bed2)
    return oname

def stat_annotates(ofile,ostat,titles):
    if ofile.name == "<stdout>":
        return None
    cmd = f'''cat {ofile.name} | cut -f 7- | sort | uniq -c | tr -s ' ' |sed -e 's/^ *//'  -e 's/ /\t/g'  '''
    cmd_out = subprocess.run(cmd,shell=True,capture_output=True,text=True)
    df = pd.read_csv(io.StringIO(cmd_out.stdout),sep="\t",names=["count"]+list(titles))
    ostat.write(df.to_csv(sep="\t",index=False))

def annotate_bed_by_feature(bed,*features):
    bed_new = NamedTemporaryFile("w",delete=False)
    bed_new_name = bed_new.name
    bed_new.close()
    
    assert len(features) > 0
    cmd = f'cat {bed} | sort -k1,1 -k2,2n | bedtools intersect -c -a - -b {features[0]}'

    if len(features) > 1:
        for b in features[1:]:
            cmd = f'''{cmd} | bedtools intersect -c  -a - -b {b} '''
    
    cmd = f'''{cmd} | sort -k4n > {bed_new_name}'''
    os.system(cmd)
    os.remove(bed)

    return bed_new_name

@click.command()
@click.option("-b","--bedpe", type=click.Path(exists=True), help="input bedpe")
@click.option("-r","--region", type=click.Path(exists=True),multiple=True, required=True,help="regions to annotate")
@click.option("-t","--title",type=str, multiple=True, required=False,help="title for your regions")
@click.option("-o","--ofile",type=click.File("w"),default=sys.stdin, help="where to output, default is stdout.\n info: in feature col, 0:no anchor; 1:R1 anchor; 2:R2 anchor; 3:R1+R2 anchor ")
@click.option("-s","--ostat",type=click.File("w"),default = sys.stderr, help="where to output stats info.default is stderr")
def process(bedpe,region,title,ofile,ostat):

    flag = True
    if len(title) == 0:
        flag = False
    elif len(title) != len(region):
        flag = False
    if not flag:
        title = [f"attr_{i}" for i in range(len(region))]
    
    bed1, bed2 = split_bedpe(bedpe)
    bed1 = annotate_bed_by_feature(bed1,*region)
    bed2 = annotate_bed_by_feature(bed2,*region)
    combine_bed_to_bedpe(bed1,bed2,len(region),ofile)
    stat_annotates(ofile,ostat,title)

if __name__ == '__main__':
    process()