from NanoPreP.preptools.Optimizer import Optimizer
from NanoPreP.seqtools.FastqIO import FastqIndexIO
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import json
import os

def my_jointplot(data, x, y, hue, palette, title):
    # set figure size/dpi
    plt.figure(
        figsize=(10, 10),
        dpi=300
    )
    
    # create jointgrid
    g = sns.JointGrid(
        data=data,
        x=x,
        y=y,
        hue=hue,
        space=.1,
        ratio=2,
        palette=palette,
        marginal_ticks=True,
    )

    # plot stripplot to joint margin
    sns.stripplot(
        data=data,
        x=x,
        y=y,
        hue=hue,
        alpha=1,
        jitter=.1,
        s=1,
        ax=g.ax_joint,
        palette=palette,
        legend=False,
        dodge=False,
    )

    # plot histogram to x margin
    sns.histplot(
        data=data,
        x=data[x] + 1,
        hue=hue,
        ax=g.ax_marg_x,
        legend=False,
        log_scale=True,
        palette=palette,
    )

    # plot count plot to y margin
    sns.countplot(
        data=data,
        y=data[y],
        hue=hue,
        ax=g.ax_marg_y,
        dodge=False,
        alpha=.6,
        palette=palette,
        hue_order=[k for k, v in palette.items()],
        width=.5,
        edgecolor="black",
    )
    g.ax_marg_y.get_legend().remove()

    # set title
    g.ax_marg_x.set(title=title)
    
    return g

# get command line arguments
FILE_IN=sys.argv[1]
DIR_OUT=sys.argv[2]
P5_SENSE=sys.argv[3]
P3_SENSE=sys.argv[4]

# set global variables
PLENS = [20, 40, 60]
N_IQR = [1, 2, 3, 4, 5]
N_PROCESS = 100
BETAS = [.1, .5, 1, 2]
BASENAME = FILE_IN.split("/")[-1].split(".")[0]

# set sample color
palette = {
    "positive": "#FFAF33",
    "negative": "#2190FF",
    .1: "darkred",
    .5: "red",
    1: "black",
    2: "lawngreen"
}

# initial optimizer
optimizer = Optimizer(
    p5_sense=P5_SENSE,
    p3_sense=P3_SENSE
)

# initial fastq index (interable)
fq_iter = FastqIndexIO(FILE_IN)

# get optimzied parameters under different `beta`
optm_params = {}
for beta in BETAS:
    optm_params[beta] = optimizer.optimze(
        fq_iter = fq_iter,
        plens = PLENS,
        n_iqr=N_IQR,
        processes=N_PROCESS,
        target="fscore",
        beta=beta
    )
    
# make output directory
os.makedirs(DIR_OUT, exist_ok=True)    

# write optimized parameters to json file
with open(f"{DIR_OUT}/{BASENAME}.json", "w") as f:
    f.write(json.dumps(optm_params, indent=4))
    
# plot loc/pid distribution under each primer length 
# with addition of cutoff lines optimized under different `beta`
for plen in PLENS:
    # plot left
    data = optimizer.getPLC(fq_iter, plen, "left", N_PROCESS)
    data = data.sort_values(["pid"], ascending=False)
    data["pid_str"] = [str(i) for i in data["pid"]]
    g = my_jointplot(
        data,
        x="loc",
        y="pid_str",
        hue="cls",
        palette=palette,
        title=f"{BASENAME}_left_plen{plen}"
    )
    
    # check if there is optimized cutoff
    for beta, optm_param in optm_params.items():
        if optm_param["left"]["plen"] == plen:
            g.refline(
                x=optm_param["left"]["loc"],
                y=str(optm_param["left"]["pid"]),
                color=palette[beta],
                linewidth=1,
                linestyle="--"
            )
    g.savefig(f"{DIR_OUT}/{BASENAME}_left_plen{plen}.png")
    plt.close()

    # plot right
    data = optimizer.getPLC(fq_iter, plen, "right", N_PROCESS)
    data = data.sort_values(["pid"], ascending=False)
    data["pid_str"] = [str(i) for i in data["pid"]]
    g = my_jointplot(
        data,
        x="loc",
        y="pid_str",
        hue="cls",
        palette=palette,
        title=f"{BASENAME}_right_plen{plen}"
    )
    for beta, optm_param in optm_params.items():
        if optm_param["right"]["plen"] == plen:
            g.refline(
                x=optm_param["right"]["loc"],
                y=str(optm_param["right"]["pid"]),
                color=palette[beta],
                linewidth=1,
                linestyle="--"
            )
    g.savefig(f"{DIR_OUT}/{BASENAME}_right_plen{plen}.png")
    plt.close()