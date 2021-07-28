#!/usr/env/bin python
import os
import sys
import re

import numpy as np
import pandas as pd
import argparse


class Motif:
    def __init__(self,meme:str):
        self.text=meme.strip()
        self.header="\n".join(self.text.split("\n")[0:3])
        self.id=self.header.split()[1]
        self.project=self.id.split("-")[-1]
        self.length=int(re.findall("w= (\d+)",self.header)[0])
        self.nsite=int(re.findall("nsites= (\d+)",self.header)[0])
        self.E=float(re.findall("[PE]= (.*)\s",self.header)[0])
        self.pwm_str="\n".join(self.text.split("\n")[3:]).strip()
        self.pwm=self.parser_pwm()
        
        self.get_score()
    
    def parser_pwm(self):
        pwm=[]
        for line in self.pwm_str.split("\n"):
            if len(line)==0:
                continue
            a=[float(i) for i in line.split()]
            pwm.append(a)
        pwm=np.array(pwm)
        return pwm
    
    def get_score(self):
        p=0
        for i in range(self.pwm.shape[0]):
            line=self.pwm[i,:]
            line=line+0.001
            
            line=line/np.sum(line)
            p+=np.log2(np.max(line)/0.25)
            self.pwm[i,:]=line
        self.score=p
        return p
    
    def out(self,alias=""):
        if len(alias) == 0:
            alias=self.id
        header=f">{self.project}\t{alias}\t{self.score:f}\n"
        body=pd.DataFrame(self.pwm).to_csv(sep="\t",index=False,header=False)+"\n"
        #body=self.pwm_str.replace(" ","\t")+"\n"
        return header+body

class MotifFile:
    def __init__(self, iname, oname,alias):
        self.iname=iname
        self.oname=oname
        self.motifs=[]

        if len(alias.strip())==0:
            self.alias=os.path.split(self.iname)[-1].split(".")[0]
        else:
            self.alias=alias

        self.parse()

    def parse(self):
        with open(self.iname) as f:
            text=f.read()
        motifs=re.findall("(MOTIF[\s\S]+?)(?=MOTIF|URL|\*+)",text)
        for i in motifs:
            try:
                self.motifs.append(Motif(i))
            except Exception as e:
                print(f"##error at \n{i}".replace("\n","\n#"))
                raise(e)
        #self.motifs=[Motif(i) for i in motifs]

    def out(self, pvalue= 1e-5):
        if self.oname == "std":
            f=sys.stdout
        else:
            f=open(self.oname,'w')
        
        c=0
        for i in self.motifs:
            if i.E < pvalue:
                c+=1
                f.write(i.out(f"{c}-{self.alias}").strip())
                f.write("\n") ## no blank lines are allowd!
        f.close()

if __name__ == "__main__":
    opt=argparse.ArgumentParser("trans-format of meme-format motif to homer-format")
    opt.add_argument("-i", "--ifile", dest="ifile", type=str, action="store", required=True, help="the input meme-format file\n")
    opt.add_argument("-o", "--ofile", dest="ofile", type=str, action="store", required=True, help= "output file name")

    opt.add_argument("-p","--pvalue", dest="pvalue", type=float, action="store", required=False, default=1e-5, help="q/p vaule you want to filter motfis")
    opt.add_argument("-a", "--alias", dest="alias", type=str, action="store", required=False, default="", help="the alias of each motif. the 2nd part of header. ")

    arg=opt.parse_args()

    m=MotifFile(arg.ifile, arg.ofile,arg.alias)
    m.out(arg.pvalue)