"""
Python script to plot the partial charge of an atom in a two-atom-molecule depending on the distance between the atoms.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # type: ignore

# Constants
CM_TO_INCH = 1 / 2.54

# Define a dictionary with the partial charges of the atoms for each method

charges: dict[str, dict[float, float]] = {}
charges.update(
    {
        "EEQ": {
            2.87: 4.8380136419754471e-03,
            4.0: -6.7526919789808737e-03,
            5.0: -9.7851229732465184e-03,
            6.0: -0.009594,
            7.0: -0.009451,
            8.0: -0.009346,
            9.0: -0.009266,
            10.0: -0.009203,
        },
        "CEH-v2": {
            2.87: -1.3419657058664991e-02,
            4.0: -4.0178820249867464e-02,
            5.0: -8.8992126111635678e-02,
            6.0: -0.161494,
            7.0: -0.193031,
            8.0: -0.198278,
            9.0: -0.198879,
            10.0: -0.198937,
        },
        "CEH-v1": {
            2.87: -1.5773404484514392e-01,
            4.0: -7.7287983177519359e-01,
            5.0: -9.9182032348117544e-01,
            6.0: -0.999850,
            7.0: -0.999982,
            8.0: -0.999983,
            9.0: -0.999983,
            10.0: -0.999983,
        },
        "ωB97M-V/def2-TZVPPD": {
            2.87: 0.003324,
            4.0: 0.015204,
            5.0: 0.013791,
            6.0: 0.010426,
            7.0: 0.008752,
            8.0: 0.006940,
            9.0: 0.004972,
            10.0: 0.003447,
        },
    }
)

gaps: dict[str, dict[float, float]] = {}
gaps.update(
    {
        "CEH-v2": {
            2.87: 0.370545,
            4.0: 0.240656,
            5.0: 0.112569,
            6.0: 0.044947,
            7.0: 0.018619,
            8.0: 0.011434,
            9.0: 0.010338,
            10.0: 0.010227,
        },
        "CEH-v1": {
            2.87: 0.081452,
            4.0: 0.027180,
            5.0: 0.022341,
            6.0: 0.022234,
            7.0: 0.022234,
            8.0: 0.022234,
            9.0: 0.022234,
            10.0: 0.022234,
        },
    }
)

AU2EV = 27.211386245988  # conversion factor from atomic units to electron volts


# Plot the partial charge of the atom in the two-atom-molecule depending on the distance between the atoms
def plot_charges(chrgs: dict[str, dict[float, float]]):
    """
    Actual plotting using seaborn.
    """
    sns.set(style="darkgrid")
    fig, ax = plt.subplots(figsize=(8 * CM_TO_INCH, 4 * CM_TO_INCH))

    # define Roboto Condensed as the default font
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = "Roboto Condensed"
    plt.rcParams["font.weight"] = "regular"

    for method, chrgdata in chrgs.items():
        x = list(chrgdata.keys())
        y = list(chrgdata.values())
        if method == "EEQ":
            color = "#fcba00"
            label = method
        elif method == "CEH-v2":
            color = "#07529a"
            label = method
        elif method == "CEH-v1":
            color = "#909085"
            label = method
        elif method == "ωB97M-V/def2-TZVPPD":
            color = "#C73C5F"
            label = "ωB97M-V"
        else:
            raise ValueError(f"Unknown method: {method}")

        ax.plot(x, y, marker="o", label=label, color=color, linewidth=1, markersize=2)
    # set y range for ax (partial charge)
    ax.set_ylim(-1.2, 0.05)

    ax2 = ax.twinx()
    ax2.set_ylabel("gap / eV", fontsize=8)
    for method, gapdata in gaps.items():
        x = list(gapdata.keys())
        # convert gaps from atomic units to electron volts
        y = [gap * AU2EV for gap in gapdata.values()]
        if method == "CEH-v2":
            color = "#07529a"
            label = "CEH-v2"
        elif method == "CEH-v1":
            color = "#909085"
            label = "CEH-v1"
        ax2.plot(
            x,
            y,
            marker="^",
            label=label,
            color=color,
            linewidth=0.75,
            markersize=1.5,
            linestyle="--",
        )
    ax2.set_ylim(0, 12.5)

    ax.set_xlabel("distance / Å", fontsize=8)
    ax.set_ylabel("partial charge on Na / $e^{-}$", fontsize=8)
    ax.xaxis.set_label_position("top")
    ax.xaxis.set_ticks_position("top")
    # set legend to middle right
    ax.legend(loc="center right", fontsize=6)
    ax.tick_params(axis="both", which="major", labelsize=6)
    ax.tick_params(left=False, top=False)
    # recuce distance between tick labels and axis
    ax.tick_params(axis="x", pad=-2)
    ax.tick_params(axis="y", pad=-2)
    ax2.tick_params(labelsize=6)
    ax2.tick_params(axis="y", pad=-2)
    ax2.tick_params(left=False, right=False)
    # reduce line width of tick symbols
    ax2.tick_params(width=0.5)
    ax2.grid(False)
    # save figure as svg
    fig.savefig("LiNa_chargeplot.svg", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    plot_charges(charges)
