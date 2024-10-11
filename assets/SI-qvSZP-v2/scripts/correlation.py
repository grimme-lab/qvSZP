"""
Script to calculate the correlation between atomic charges calculated with different methods
and plot the correlation between the reference and the calculated charges.
"""

import argparse
from itertools import combinations

import numpy as np
import seaborn as sns  # type: ignore
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

# Statistic parameters
OUTLIER_THRESHOLD = 1.0
SIMILARITY_THRESHOLD = 0.1

# Constants
CM_TO_INCH = 1 / 2.54
FFORMAT = "svg"

# Checking if the input file was supplied
parser = argparse.ArgumentParser(
    description="Plotting the correlation between Reference and PTB atomic charges."
)
parser.add_argument(
    "-i",
    "--input_file",
    required=True,
    type=str,
    help="The input file containing the DFT and PTB atomic charges.",
)
parser.add_argument(
    "--methods",
    "-m",
    nargs="+",
    required=False,
    type=str,
    help="The methods to be compared.",
)
parser.add_argument(
    "-acs",
    "--actinoids",
    default=False,
    action="store_true",
    help="Plot actinoids charges.",
)
parser.add_argument(
    "--clean",
    "-c",
    default=False,
    action="store_true",
    help="Remove all rows where all charge values have a deviation less than 0.1.",
)
parser.add_argument(
    "--plotall",
    default=False,
    action="store_true",
    help="Plot all data points (e.g., for 'mindless' comparison).",
)
parser.add_argument(
    "--debug",
    default=False,
    action="store_true",
    help="Print debug information.",
)
args = parser.parse_args()


# The input charge file should have the following format:
"""
CID,Method,Charge,Atom number,Atom type
6695,EEQ,-0.32361661784948015,1,8
6695,EEQ,-0.48401285989710696,2,8
6695,EEQ,-0.40118257329635776,3,8
...
6695,GFN1-xTB,-0.06406547372206295,23,6
6695,GFN1-xTB,0.01746306275701065,24,6
6695,GFN1-xTB,-0.11005238125478034,25,6
6695,GFN1-xTB,-0.23135872399329793,26,6
...
9993533,wB97M-V,0.042646,36,1
9993533,wB97M-V,0.042513,37,1
"""

REF_COL = "wB97M-V"
if args.actinoids:
    REF_COL = "HirshfeldAPC"

df = pd.read_csv(args.input_file, comment="#")
# get all unique Method values
methods = df["Method"].unique()
df = df.pivot(
    index=["CID", "Atom number", "Atom type"], columns="Method", values="Charge"
).reset_index()
# in the actinoid mode, keep only lines with actinoids (Z > 88)
if args.actinoids:
    df = df[df["Atom type"] > 88]
df.to_csv("charges_pivoted.csv", index=False)

# create a deep copy of the dataframe to keep the original data
df_noclean = df.copy(deep=True)

print(f"Number of raw data points: {len(df_noclean)}")

# if method is EspalomaCharge, delete all rows with an EspalomaCharge of -10.0 (dummy value)
if "EspalomaCharge" in args.methods:
    df = df[df["EspalomaCharge"] > -9.0]

# create combinations of args.methods and REF_COL
methods_to_plot: list[str] = [
    method for method in methods if method in args.methods
] + [REF_COL]
plot_combs = list(combinations(methods_to_plot, 2))

outliers: dict = {}  # store the outliers in a dictionary with the method as key and the outliers (CID) as value
for method in args.methods:
    outliers[method] = []
if args.clean:
    print("\nCleaning the data...")
    # remove all rows where all values of the remaining method columns
    # are between 0.0 and 0.5 AND the respecive deviation is less than 0.1
    # include a tqdm progress bar
    with tqdm(total=len(df)) as pbar:
        for index, row in df.iterrows():
            pbar.update(1)
            # check for outliers (large discrepancies between the reference and the calculated charges)
            for method in args.methods:
                if abs(row[method] - row[REF_COL]) > OUTLIER_THRESHOLD:
                    print(
                        f"Outlier found: CID {row['CID']}, atom {int(row['Atom number'])} of type {int(row['Atom type'])} {method}"
                    )
                    print(
                        f"Reference: {row[REF_COL]:.4f}; Calculated: {row[method]:.4f}"
                    )
                    outliers[method].append(row["CID"])
            # if the row is an outlier for any of the methods, remove it
            if any(
                abs(row[method] - row[REF_COL]) > OUTLIER_THRESHOLD
                for method in args.methods
            ):
                df.drop(index, inplace=True)
                continue
            # remove rows with atoms of type "1" (hydrogen) as they are not interesting
            if row["Atom type"] == 1:
                df.drop(index, inplace=True)
                continue
            # if the deviation between all combined methods is less than threshold, remove the row
            if all(
                abs(row[comb[0]] - row[comb[1]]) < SIMILARITY_THRESHOLD
                for comb in plot_combs
            ):
                df.drop(index, inplace=True)
                continue
    df.reset_index(drop=True, inplace=True)
    print(f"\nNumber of data points after filtering: {len(df)}")

# for each method, add a new column to the data frame containing the deviation from the reference column
for method in args.methods:
    df[method + "_dev"] = df[method] - df[REF_COL]
    # print the 10 largest deviations
    print(f"\n10 largest deviations for {method}:")
    print(df.nlargest(10, method + "_dev"))
    # print information about the outliers
    if outliers[method]:
        print(f"\nOutliers for {method}:")
        for cid in set(outliers[method]):
            print(f"  CID {cid}")
        print(f"Number of outliers: {len(outliers[method])}")

# do some statistics with the chosen methods
pearsons = {}
active_methods = [method for method in methods if method in args.methods]
for method in active_methods:
    print(f"\nStatistics for {method}:")
    # if method is EspalomaCharge, delete all rows with an EspalomaCharge of -10.0 (dummy value)
    tmp_df = df_noclean.copy(deep=True)
    if method == "EspalomaCharge":
        tmp_df = tmp_df[tmp_df["EspalomaCharge"] > -9.0]
    # remove outliers of the specific method (get outliers from the dictionary outliers)
    tmp_df = tmp_df[~tmp_df["CID"].isin({item for item in outliers[method]})]

    ## > Statistics
    # Pearson correlation
    pearson = tmp_df[[REF_COL, method]].corr(method="pearson")
    print(
        3 * " "
        + f"Pearson correlation coefficient (PCC) for {REF_COL} vs {method}: {pearson.iloc[0, 1]:.4f}"
    )
    pearsons[method] = pearson.iloc[0, 1]
    # Spearman correlation
    spearman = tmp_df[[REF_COL, method]].corr(method="spearman")
    print(
        3 * " "
        + f"Spearman correlation coefficient (SCC) for {REF_COL} vs {method}: {spearman.iloc[0, 1]:.4f}"
    )
    # MAE
    mae = (tmp_df[REF_COL] - tmp_df[method]).abs().mean()
    print(3 * " " + f"MAE for {REF_COL} vs {method}: {mae:.4f}")
    # RMSE
    rmse = np.sqrt(((tmp_df[REF_COL] - tmp_df[method]) ** 2).mean())
    print(3 * " " + f"RMSE for {REF_COL} vs {method}: {rmse:.4f}")

    print(
        "LaTeX table row for the header: 'Method & PCC & SCC & MAE & RMSE & \\# outliers \\\\'"
    )
    print(
        f"{method} & {pearson.iloc[0, 1]:.3f} & {spearman.iloc[0, 1]:.3f} & {mae:.3f} & {rmse:.3f} & {len(outliers[method])} \\\\"
    )
    if args.debug:
        # save the tmp_df of the method for debugging the statistics
        tmp_df.to_csv(f"debug_stats_{method}.csv", index=False)

plt.figure(figsize=(8 * CM_TO_INCH, 6 * CM_TO_INCH))
sns.set(style="darkgrid")

# define Roboto Condensed as the default font
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Roboto Condensed"
plt.rcParams["font.weight"] = "regular"

ALPHA_DOT = 0.3
POINTSIZE = 2
if args.plotall:
    ALPHA_DOT = 0.75
    POINTSIZE = 4

list_method_specific_settings = []
method_specific_settings: dict[str, str | float] = {}
for method in args.methods:
    if method == REF_COL:
        continue
    if method == "CEH-v2":
        method_specific_settings = {
            "color": "#07529a",
            "marker": "o",
            "edgecolor": "black",
        }
    if method == "CEH-v1":
        method_specific_settings = {
            "color": "#909085",
            "marker": "D",
            "edgecolor": "black",
        }
    elif method == "EEQ":
        method_specific_settings = {
            "color": "#fcba00",
            "marker": "p",
            "edgecolor": "black",
        }
    elif method == "GFN1-xTB":
        method_specific_settings = {
            "color": "#A367C2",
            "marker": "s",
            "edgecolor": "black",
        }
    elif method == "GFN2-xTB":
        method_specific_settings = {
            "color": "#C73C5F",
            "marker": "^",
            "edgecolor": "black",
        }
    elif method == "EspalomaCharge":
        method_specific_settings = {
            "color": "#7AC284",
            "marker": "<",
            "edgecolor": "black",
        }
    list_method_specific_settings.append(method_specific_settings)
    print(f"Plotting '{method}' data...")
    # general settings
    settings = {
        "x": REF_COL,
        "y": method,
        "data": df,
        "s": POINTSIZE,
        "linewidth": 0.1,
        "alpha": ALPHA_DOT,
        # add PCC value to the legend
        # (use the pearsons values calculated before, always with respect to the reference column)
        "label": method,  #  + f" ({pearsons[method]:.3f})",
    }
    # combine the settings with the specific marker and edgecolor
    settings.update(method_specific_settings)
    p = sns.scatterplot(**settings)

    tmp_df = df_noclean.copy(deep=True)
    if method == "EspalomaCharge":
        tmp_df = tmp_df[tmp_df["EspalomaCharge"] > -9.0]
    # remove outliers of the specific method (get outliers from the dictionary outliers)
    tmp_df = tmp_df[~tmp_df["CID"].isin({item for item in outliers[method]})]
    # Add a linear regression line
    sns.regplot(
        x=REF_COL,
        y=method,
        data=tmp_df,
        color=method_specific_settings["color"],
        scatter=False,
        line_kws={"linestyle": "--", "linewidth": 0.5},
        ci=95,
        # elongate the line to the plot borders
        truncate=False,
    )


# Adding a 1:1 line
plt.plot(
    [df[REF_COL].min() - 0.1, df[REF_COL].max() + 0.1],
    [df[REF_COL].min() - 0.1, df[REF_COL].max() + 0.1],
    ls="-",
    lw=0.5,
    color="black",
    # create legend entry
    label="1:1 line",
)

# Customizing plot
xlab = "ωB97M-V/def2-TZVPPD Hirshfeld charge / a.u."
if args.actinoids:
    xlab = "ωB97M-V/ma-def-TZVP Hirshfeld charge / a.u."
plt.xlabel(
    xlabel=xlab,
    fontsize=8,
    color="black",
)
plt.ylabel("calc. atomic charge / a.u.", fontsize=8, color="black")
plt.xticks(fontsize=8, color="black")
plt.yticks(fontsize=8, color="black")

# get current legend's handles and labels
old_handles, old_labels = p.get_legend_handles_labels()
new_handles = []
for i, method_specific_settings in enumerate(list_method_specific_settings):
    new_handle = plt.Line2D(
        [],
        [],
        color=method_specific_settings["color"],  # type: ignore
        marker=method_specific_settings["marker"],  # type: ignore
        alpha=0.75,
        linestyle="",
        markeredgecolor="black",
        markersize=1.5,  # set markersize for better visibility
        markeredgewidth=0.2,
    )
    new_handles.append(new_handle)

# update the legend
legend = plt.legend(
    new_handles + old_handles[-1:],
    old_labels,
    loc="upper left",
    fontsize=8,
    frameon=True,
    title="",
    fancybox=True,
    markerscale=2,
    title_fontsize=8,
)
legend.get_frame().set_facecolor("white")
legend.get_frame().set_edgecolor("black")
# set linewidth of the legend frame
legend.get_frame().set_linewidth(0.5)

# Customizing plot margins
plt.subplots_adjust(top=0.9, right=0.9, bottom=0.1, left=0.1)

# Generate output file name from the method names
filename = "plot_" + "_".join(args.methods) + "_correlation" + "." + FFORMAT

# Saving the plot
plt.savefig(filename, format=FFORMAT, bbox_inches="tight", dpi=600)
