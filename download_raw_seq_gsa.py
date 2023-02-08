import os 
import sys 
import click 
import pathlib

import pandas as pd
import numpy as np

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())
logger.handlers[0].setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s : %(message)s"))

readme = """this script aims to download fastq files from GSA database(https://ngdc.cncb.ac.cn/) """
GSA_BASE_URL = "download.cncb.ac.cn:gsa"

def load_gsa_meta(metafile):
    sample_table = pd.read_excel(metafile, sheet_name="Sample", header=0, index_col=0)
    experiment_table = pd.read_excel(metafile, sheet_name="Experiment", header=0, index_col=0)
    run_table = pd.read_excel(metafile, sheet_name="Run", header=0, index_col=0)

    return sample_table, experiment_table, run_table

def process_metatable_pe(sample_table, experiment_table, run_table, gsa_accession, accessions = None):
    if accessions == None:
        pass
    else:
        if isinstance(accessions, list):
            experiment_table = experiment_table.query("Accession in @accessions")
        else:
            raise ValueError("accessions should be a list of accession, but got {}".format(type(accessions)))
    experiment_table = experiment_table.query("Layout == 'PAIRED'")

    table_tmp = experiment_table[["Accession", "Experiment title", "BioSample accession"]]\
        .rename(
            columns = {
                "Accession": "Experiment accession",
            }
        )
    table_tmp = table_tmp.merge(run_table).rename(
        columns = {
            "Accession": "Run accession",
        }
    ).drop(
        columns = [
            "Run title",
            "BioProject accession",
        ]
    )

    table_tmp["file1"] = table_tmp["File name 1"].apply(lambda x: x.split("(")[0].strip() + ".gz")
    table_tmp["file2"] = table_tmp["File name 2"].apply(lambda x: x.split("(")[0].strip() + ".gz")

    # prepare md5 info
    files = []
    md5s = []
    for i, row in table_tmp.iterrows():
        files.append(row["file1"])
        files.append(row["file2"])
        md5s.append(row["MD5 checksum 1"])
        md5s.append(row["MD5 checksum 2"])
    table_md5 = pd.DataFrame({"md5": md5s, "file": files})

    # prepare seqpair info
    infos = []
    for g, df in table_tmp.groupby("Experiment accession"):
        info = dict()
        info["Experiment accession"] = g
        info["Experiment title"] = ",".join(list(set(df["Experiment title"])))
        info["file1"] = ",".join(list(set(df["file1"])))
        info["file2"] = ",".join(list(set(df["file2"])))
        infos.append(info)
    table_seqpair = pd.DataFrame(infos)

    # prepare download cmd 
    files = []
    urls = []
    for i, row in table_tmp.iterrows():
        url = f"{GSA_BASE_URL}/{gsa_accession}/{row['BioSample accession']}/{row['Experiment accession']}/{row['Run accession']}/{row['file1']}"
        files.append(row["file1"])
        urls.append(url)
        url = f"{GSA_BASE_URL}/{gsa_accession}/{row['BioSample accession']}/{row['Experiment accession']}/{row['Run accession']}/{row['file2']}"
        files.append(row["file2"])
        urls.append(url)
    table_download = pd.DataFrame({"file": files, "url": urls})


    return table_md5, table_seqpair, table_download


def wrap_download(url, odir, ascp, idfile, user):
    # ascp -P 33001 -i ./source/aspera01.openssh -QT -l100m -k1 -d aspera01@
    if not odir.endswith("/"):
        odir = odir.strip() + "/"
    return f"{ascp} -P 33001 -i {idfile} -QT -l 300m -k1 -d {user}@{url} {odir}"
    
def generate_cmd(table_download, odir, ascp, idfile, user, md5):
    cmds = []
    for i, row in table_download.iterrows():
        cmd = wrap_download(row["url"], odir, ascp, idfile, user)
        cmds.append(cmd)
    cmds.append("echo 'download finished'")
    cmds.append("echo 'start to check md5'")
    md5log = md5 + ".log"
    md5log = os.path.abspath(md5log)
    cmd = f"cd {odir} && md5sum -c {md5} > {md5log} 2>&1  && cd - "
    cmds.append(cmd)
    return cmds

@click.command()
@click.option("-m", "--metafile", type=click.Path(exists=True), required=True, help="the metafile you want to download, used to down metafile \n the file name has to be gsa accession + '.xlsx' \n")
@click.option("-a", "--accession", type=click.File(mode="r"), required=False, default = None, help="accessions you want to download, used to download metafile \nThe parameter should be a filename, with one accession per line. \n")
@click.option("-i", "--idfile", type=click.Path(exists=True), required=False, default="~/aspera.openssh", help="the aspera key file you want to use, used to download metafile \n")
@click.option("-o", "--outdir", type=click.Path(resolve_path=True), required=False, default="rawdata", help="where to save downloaded files\n")
@click.option("-c", "--cmd", type=click.Path(), required=False, default="download.sh", help="where you want to output download cmd scripts\n")
@click.option("-d", "--md5", type=click.Path(), required=False, default="md5s.txt", help="where to save your md5 info")
@click.option("-s","--seqpair", type = click.Path(dir_okay=False), required = False, default = None, help = "the file to save seqpair info. Default is stem of metafile + '.seqpair.txt'")
@click.option("--ascp", type=click.Path(exists=False), required=False, default="ascp", help="the ascp path you want to use, used to download metafile \n")
@click.option("--user", type=str, required=False, default="aspera01", help="the ascp user you want to use, used to download metafile \n")
def main(metafile, accession, idfile, outdir, cmd, md5, seqpair, ascp, user):
    """this script aims to download fastq files from GSA database(https://ngdc.cncb.ac.cn/) """
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    logger.info("loading metafile...\n")

    if accession != None:
        accessions = [i.strip() for i in accession.readlines()]
    else:
        accessions = None
    logger.info("selected accessions loaded: %s" % accessions)
    
    gsa_accession = pathlib.Path(metafile).stem
    sample_table, experiment_table, run_table = load_gsa_meta(metafile)

    table_pe_md5, table_pe_seqpair, table_pe_download = process_metatable_pe(sample_table, experiment_table, run_table, gsa_accession, accessions = accessions)
    logger.info("pe metafile processed")

    table_pe_seqpair["idir"] = outdir
    if seqpair == None:
        seqpair = pathlib.Path(metafile).stem + ".seqpair.txt"
        logger.info("seqpair file not specified, use default: %s" % seqpair)

    table_pe_seqpair[["Experiment accession", "Experiment title", "idir", "file1", "file2"]].to_csv(seqpair, sep=";", index=False, header=False)
    logger.info("pe seqpair info saved to %s" % seqpair)

    
    table_pe_md5.to_csv(md5, sep="\t", index=False, header=False)
    logger.info("md5 info saved to %s" % md5)

    logger.info(f"use user: {user} to download files...")
    cmds = generate_cmd(table_pe_download, outdir, ascp, idfile, user, md5)
    with open(cmd, "w") as f:
        f.write("\n\n".join(cmds))
    
    logger.info("download cmd saved to %s" % cmd)


if __name__ == "__main__":
    main()
