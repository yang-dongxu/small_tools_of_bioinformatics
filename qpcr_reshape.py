#python


import numpy as np
import pandas as pd 
from matplotlib import pyplot as plt
import seaborn as sns
import logging
from scipy import stats as st
import argparse

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO,format='HUHU~: %(asctime)s - %(name)s - %(levelname)s - %(message)s')


def opt():
    args = argparse.ArgumentParser()
    args.add_argument("-i", "--ifile",default = "sample.xls", required=True, help="input file, excel format!")
    args.add_argument("-o", "--ofile",default=None, help="output file, excel format!")
    args.add_argument("-t","--threshold",default=0.3,type = float, help="threshold!")
    args.add_argument("-r","--reference",default="GAPDH", help="inner reference target name")
    args.add_argument("-d","--delete_samples",dest="ds",default= [], action="append", help="delete samples. If multi samples need to be deleted, use -d multiple times, comma separated like: -d s1,s2")
    args.add_argument("-s","--suffix-length", default=1, action="store", help="length of suffix to distinguish replicates. For example: sA_1, sA_2 are replicats, so -s is 1 or 2.")

    return args.parse_args()


# @click.command()
# @click.option("-i", "--ifile", required=True,type=click.Path(exists=True), help="input file, excel format!")
# @click.option("-o", "--ofile",default=None, help="output file, excel format!")
# @click.option("-t","--threshold",default=0.3, help="threshold!")
# def opt(ifile,ofile,threshold):
#     args = namedtuple("args",["ifile","ofile","threshold"])
#     args.ifile = ifile
#     args.ofile = ofile
#     args.threshold = threshold
#     return args


def filter(values, positions, title, thre=0.5):
    info = dict(zip(positions, values))
    assert len(values) == len(positions)
    if len(positions) <= 1:
        logger.error(f"Not enough data for {title}!")
        return np.nan, ["Not enough data"]
    for p in positions:
        p_remain = list(set(positions) - set([p]))
        v_remain = [info[p] for p in p_remain]
        v_remain_mean = np.mean(v_remain)
        v = info[p]
        if max(v_remain) - min(v_remain) < thre:
            if abs(v_remain_mean - v) < thre:
                continue
            else:
                logger.warning(f"{title} {p} were filtered out!")
                return filter(values=v_remain, positions=p_remain, title=title, thre=thre)
        else:
            continue
    return np.mean(values), positions


def load(ifile):
    df = pd.read_excel(ifile, sheet_name="Results", skiprows=42)[
        ["Well", "Well Position", "Sample Name", "Target Name", "CT"]
        ].dropna(subset=["Target Name","CT"], how="all")
    if len(df) < 1:
        raise Exception("No data found in Results sheet of {}".format(ifile))
    df["Well_row"] = pd.Categorical(df["Well Position"].apply(lambda x:x[0]))
    df["CT"] = pd.to_numeric(df["CT"],errors="coerce")
    df["Well_col"] = df["Well Position"].apply(lambda x:int(x[1:]))
    df["Sample_id"] =df["Well_row"].astype(str) + "_" + df["Sample Name"] + "_" + df["Target Name"]
    
    return df


def samples_delete(df, ds):
    ds = [d.split(",") for d in ds]
    ds = [j.strip() for d in ds for j in d]
    ds = [d for d in ds if len(d) > 0]
    for d in ds:
        logger.warning(f"Ignore sample {d}!")
    df = df[~df["Sample Name"].isin(ds)].copy()
    return df


def filter_df(df, thre):
    list_df_sub = []
    for g, df_sub in df.groupby("Sample_id"):
        df_s = df_sub.sort_values("Well_col")
        mean, positions = filter(
            values=df_s["CT"].values,
            positions=df_s["Well Position"].values,
            title=g,
            thre=thre
            )
        row = dict(
            Sample_id = g,
            Sample_name = df_s["Sample Name"].values[0],
            Target_name = df_s["Target Name"].values[0],
            Mean = mean,
            Positions = "; ".join(positions)
        )
        list_df_sub.append(row)
    df_mean = pd.DataFrame(list_df_sub)
    return df_mean

def normalize(df_mean_wide,ref = "GAPDH"):
    assert ref in df_mean_wide.columns
    df = df_mean_wide.copy()
    for col in df.columns:
        df[col] = (df_mean_wide[ref]-df_mean_wide[col] ).apply(lambda x: 2**x)
    return df


def visualize(df_mean_wide, oname, ref = "GAPDH", figsize=(1,1),vmin = -3, vmax = 3):
    fig_width = df_mean_wide.shape[1] * figsize[0]
    fig_height = df_mean_wide.shape[0] * figsize[1]
    figsize = (fig_width, fig_height)
    df = df_mean_wide.fillna(0).drop(columns = [ref])
    g = sns.clustermap(df, cmap = "coolwarm", figsize = figsize, vmin = vmin, vmax = vmax,z_score = 1, yticklabels = True)

    ax = g.ax_heatmap
    fig = ax.figure
    ax.set_xlabel("Targets")
    ax.set_ylabel("Samples")
    fig.suptitle("Normalized Expression\n(z-score by col, na filled with 0)")
    g.savefig(oname)

    return None

def format_pvalue(p):
    if p< 0.001:
        f = "***"
    elif p< 0.01:
        f = "**"
    elif p< 0.05:
        f = "*"
    else:
        f = "n.s."
    return f

def clean_reps(name, sl=1):
    n = name[:-sl]
    if len(n) ==0:
        n = "ctrl"
    return n

def stats(df_mean_wide, ref = "GAPDH", sl=1):
    df_a = df_mean_wide.copy()
    list_stat = []
    for col in df_a.columns:
        if col == ref:
            continue
        df = df_a[[col]].reset_index().dropna()
        df["sample"] = df["Sample_name"].apply(lambda x:clean_reps(x,sl=sl))

        samples = df["sample"].unique().tolist()
        for s1 in samples:
            for s2 in samples:
                if s1 == s2:
                    continue
                v1 = df.query("sample == @s1")[col]
                v2 = df.query("sample == @s2")[col]
                fold_change = np.nanmean(v1) / np.nanmean(v2)
                fold_median = np.nanmedian(v1) / np.nanmedian(v2)

                p = st.ttest_ind(v1, v2)[1]
                mark = format_pvalue(p)
                sign = f"{fold_change:.2f}:{fold_median:.2f} ({mark})"
                if fold_change < 1:
                    mark = f"- {mark}"

                list_stat.append(dict(
                    target = col,
                    sample1 = s1,
                    sample2 = s2,
                    p = p,
                    exp1 = np.nanmean(v1),
                    exp2 = np.nanmean(v2),
                    fold_mean = fold_change,
                    fold_median = fold_median,
                    mark = mark,
                    sign = sign
                ))
    df_stat = pd.DataFrame(list_stat)
    # df_stat["mark"] = df_stat["p"].apply(lambda x: format_pvalue(x))

    df_stat_sign_wide = df_stat.set_index(["target","sample1","sample2"])["sign"].unstack().reset_index()
    # df_stat_p_wide = df_stat_p_wide.fillna(1)

    df_stat_mark_wide = df_stat.set_index(["target","sample1","sample2"])["mark"].unstack().reset_index()
    # df_stat_mark_wide = df_stat_mark_wide.fillna("n.s.")

    return df_stat, df_stat_sign_wide, df_stat_mark_wide


def run(args):
    if args.ofile is not None:
        ofile = args.ofile
    else:
        ofile = args.ifile.replace(".xls","_preprocessed.xls").replace(".xls",".xlsx").replace(" ","_")

    df = load(args.ifile)
    df = samples_delete(df,args.ds)

    df_design_target = df.set_index(["Well_row","Well_col"])["Target Name"].unstack().fillna("missing")
    df_design_sample = df.set_index(["Well_row","Well_col"])["Sample Name"].unstack().fillna("missing")

    df_mean = filter_df(df,args.threshold)

    df_value_wide = df_mean.set_index(["Sample_name","Target_name"])["Mean"].unstack()

    df_exp = normalize(df_value_wide, args.reference)

    df_stat, df_stat_sign_wide, df_stat_mark_wide = stats(df_exp, args.reference, sl=args.suffix_length)


    visualize(df_exp, oname=ofile.replace(".xlsx", "_heatmap.pdf"), ref=args.reference)

    with pd.ExcelWriter(ofile,mode = "w") as f:
        df_design_target.to_excel(f, sheet_name="Design_target")
        df_design_sample.to_excel(f, sheet_name="Design_sample")
        df_mean.to_excel(f, sheet_name="CT_detail")
        df_value_wide.to_excel(f, sheet_name="CT_wide")
        df_exp.to_excel(f, sheet_name="Exp")
        df_stat.to_excel(f, sheet_name="Stat",index=False)
        df_stat_sign_wide.to_excel(f, sheet_name="Stat_sign",index=False)
        df_stat_mark_wide.to_excel(f, sheet_name="Stat_mark",index=False)
    logger.info("Done! Results saved in {}".format(ofile))
    logger.info("See you ~")


if __name__ == "__main__":
    args = opt()
    run(args)
