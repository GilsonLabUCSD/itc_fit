import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm

# 01/31/19 Rejection of bad fits based on sum of squared errors
# Reference: http://www.isbg.fr/IMG/pdf/microcal-itc200-system-user-manual-malvern.pdf
# ITC settings
V0 = 0.202 / 1000  # Cell Volume (L)
M0 = 0.005  # Cell Concentration (M)
X0 = 0.07500  # Injectant Concentration (M)
syringe_error = 0.02  # percent
cell_error = 0.02  # percent
heat_error = 0.01  # percent
base_error = 0.00000015  # calories


def process(file, skip):
    """
    Read in the raw ITC data and return the heat and volume injected.
    Always skip the first row which is the header from the instrument.
    Always skip the final row which is a summary of total heat and total injected
    volume from the instrument.
    Optionally skip additional injections from the beginning of the experiment.
    """

    data = np.genfromtxt(file, skip_header=1, skip_footer=1)

    heat = data[:, 0] / 1e6  # microcalories to calories
    volume = data[:, 1] / 1e6  # microliters to liters

    return heat[skip:], volume[skip:]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--file", help="Raw ITC data", required=True)
    ap.add_argument(
        "-s", "--skip", help="Number of injections to skip (default: 0)", required=False
    )
    ap.add_argument(
        "-t",
        "--temperature",
        help="Temperature of the experiment (default: 300.15 K)",
        required=False,
    )
    ap.add_argument(
        "-M",
        "--M0",
        help="Cell concentration in Molar (default: 0.005 M)",
        required=False,
    )
    ap.add_argument(
        "-X",
        "--X0",
        help="Injectant concentration in Molar (default: 0.075 M)",
        required=False,
    )
    args = vars(ap.parse_args())

    if not args["temperature"]:
        temperature = 300.15
    else:
        temperature = float(args["temperature"])
    if not args["skip"]:
        skip = 0
    else:
        skip = int(args["skip"])

    if args["X0"]:
        X0 = float(args["X0"])
    if args["M0"]:
        M0 = float(args["M0"])

    dQ, dV = process(args["file"], skip)

    # ITC = dQ / (dV * X0)

    ##################################################################################

    sampled_syringe_concentration = X0
    sampled_cell_concentration = M0

    XM = np.zeros([2, len(dQ) + 1])
    XM[0, 0] = 0
    XM[1, 0] = 0

    cumulative_volume = np.cumsum(dV)

    print(f"Going from {cumulative_volume} to...")
    # cumulative_volume = cumulative_volume[skip:]
    print(f"{cumulative_volume[skip:]} cumulative volume...")

    # New Injectant Concentration
    XM[0, 1:] = (cumulative_volume * sampled_syringe_concentration / V0) * (
        1 / (1 + (cumulative_volume / (2 * V0)))
    )

    print(f"Going from {XM[0, 1:]} to...")
    print(f"{XM[0, 1 + skip:]} guest concentration...")

    # New Cell Molecule Concentration
    XM[1, 1:] = sampled_cell_concentration * (
        (1 - cumulative_volume / (2 * V0)) / (1 + cumulative_volume / (2 * V0))
    )

    print(f"Going from {XM[1, 1:]} to...")
    print(f"{XM[1, 1 + skip:]} host concentration...")

    # Add heat error
    parITC = [
        np.random.normal(
            injection,
            abs(np.sqrt(((injection * heat_error) ** 2) + ((base_error) ** 2))),
        )
        for injection in dQ
    ]
    # Scale nominal Wiseman plot by new bootstrapped syringe concentration
    varITC = [
        injection / (volume * sampled_syringe_concentration)
        for injection, volume in zip(parITC, dV)
    ]

    print(XM[0, 1:] / XM[1, 1:])

    print(XM[0, 1 + skip :] / XM[1, 1 + skip :])
