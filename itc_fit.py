import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm

# 01/31/19 Rejection of bad fits based on sum of squared errors
# Reference: http://www.isbg.fr/IMG/pdf/microcal-itc200-system-user-manual-malvern.pdf
# ITC settings
# XM is the variable combining the concentration of the guest (X) and the host (M)
# according to the manual referenced above
# 05/08/19 Changed how the code deals with skipping of the first n data points
# Renamed variables

M0 = 0.005  # Cell Concentration (M)
X0 = 0.07500  # Injectant Concentration (M)
syringe_error = 0.02  # fraction (out of 1.0)
cell_error = 0.02  # fraction (out of 1.0)
heat_error = 0.01  # fraction (out of 1.0)
base_error = 0.000_000_15  # calories


def process(file):
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

    return heat, volume


def plot(XM, ITC, vardQ, name, dH, K, N, dG):
    """
    Plot the heat as a function of molar ratio, as measured, and the final (?) curve fit.
    I think this takes the best fit from the last bootstrapped cycle, but I'm not sure.
    """
    fig, ax = plt.subplots(1, figsize=(6 * 1.2, 6))
    ax.scatter(XM[0, 1:] / XM[1, 1:], ITC, c="k", label="ITC data")
    for index, dQ in enumerate(vardQ.T):
        ax.plot(
            XM[0, 1 + skip :] / XM[1, 1 + skip :],
            dQ,
            c="r",
            label="Equation fit" if index == 0 else "",
            alpha=0.2,
            zorder=-1,
        )

    ax.set_xlabel("Molar ratio")
    ax.set_ylabel("Heat (cal/mol)")
    ax.set_title(
        r"dH = {:9.2f} cal/mol, K = {:9.2f} M$^{{-1}}$, N = {:5.3f}, dG = {:9.2f} kcal/mol".format(
            dH, K, N, dG
        )
    )
    ax.legend(loc=4)
    fig.savefig("{}.png".format(os.path.splitext(name)[0]), bbox_inches="tight")


def fit(XMa, dH, K, N):
    """
    Fit dQ given dH, K, N.
    """
    injections = len(XMa[0])
    Q = np.zeros(injections)  # Total Heat Array Initialized
    Q[0] = 0.0  # Total Heat before injections is zero
    dQ = np.zeros(injections)
    X = XMa[0, :]
    M = XMa[1, :]
    for i in range(1, injections):  # Loop over injections
        # Total Heat Equation 9
        Q[i] = (N * M[i] * dH * V0 / 2) * (
            1
            + (X[i] / (N * M[i]))
            + 1 / (N * K * M[i])
            - np.sqrt(
                ((1 + X[i] / (N * M[i]) + 1 / (N * K * M[i])) ** 2)
                - 4 * X[i] / (N * M[i])
            )
        )

        # Change in heat normalized by amount of injectant
        # Similar to equation 10 except for the normalization factor dV*X0.
        # Now the unit is cal/mol; same unit as for exp_dQ_normalized
        dQ[i] = (Q[i] + (dV[i - 1] / V0) * ((Q[i] + Q[i - 1]) / 2) - Q[i - 1]) / (
            dV[i - 1] * X0
        )
    # We are not fitting the first skip points and and heat release before injection 1 is 0.0
    return dQ[skip + 1 :]


def bootstrap(dQ, dV, temperature, cycles=100):
    """
    Bootstrap the fitting of dQ given re-sampled uncertainties.
    """
    # Increase the number of bootstrapping runs by 10 percent. We delete the bad fits afterwards and want to fill up
    # the array so that we always use the same number of bootstrapping cycles
    realcycles = cycles
    cycles = int(cycles * 1.3)

    # All variables which are different between different cycles of the bootstrapping process and stored  start with a
    # leading 'bootstrap_'
    # bootstrap_heat = np.zeros([len(dQ)-skip, cycles]) I think this is unnecessary
    bootstrap_dH = np.zeros([cycles])
    bootstrap_K = np.zeros([cycles])
    bootstrap_N = np.zeros([cycles])
    bootstrap_SS = np.zeros([cycles])
    # We are not fitting the first 'skip' elements. Therefore the array of fitted dQ is smaller.
    # As we are only fitting the differences we can ignore the first skip elements in the fitting easily
    # For every cycle we are storing the full set of calculate dQ values
    bootstrap_dQ = np.zeros([len(dQ) - skip, cycles])
    # print(V0)
    for cycle in tqdm(range(cycles)):
        # Initial Guesses
        initial_guess = np.zeros(3)  # Guess Array
        initial_guess[0] = -1000.0  # Guess dH
        initial_guess[1] = 1000  # Guess K
        initial_guess[2] = 1.000  # Guess N

        # Picking concentrations from gaussian distribution
        sampled_syringe_concentration = np.random.normal(X0, abs(syringe_error * X0))
        sampled_cell_concentration = np.random.normal(M0, abs(cell_error * M0))

        # Find concentrations after each injection
        XM = np.zeros([2, len(dQ) + 1])

        # Concentration beore 1st injection
        XM[0, 0] = 0  # Guest concentration in cell
        XM[1, 0] = M0  # Starting host concentration in cell

        cumulative_volume = np.cumsum(dV)

        # New Injectant Concentration Equation 4
        XM[0, 1:] = (cumulative_volume * sampled_syringe_concentration / V0) * (
            1 / (1 + (cumulative_volume / (2 * V0)))
        )

        # New Cell Molecule Concentration Equation 2
        XM[1, 1:] = sampled_cell_concentration * (
            (1 - cumulative_volume / (2 * V0)) / (1 + cumulative_volume / (2 * V0))
        )

        # Add heat error
        exp_dQ = [
            np.random.normal(
                injection,
                abs(np.sqrt(((injection * heat_error) ** 2) + ((base_error) ** 2))),
            )
            for injection in dQ
        ]
        # Scale nominal Wiseman plot by new bootstrapped syringe concentration
        exp_dQ_normalized = [
            injection / (volume * sampled_syringe_concentration)
            for injection, volume in zip(exp_dQ, dV)
        ]

        # Fit the data
        # We are only fitting the experimental heat realizes after skip. Therefore, ignore the first 'skip'
        # heat releases
        # XM still has all datapoints
        try:
            fitting_variables, _ = curve_fit(
                fit, XM, exp_dQ_normalized[skip:], initial_guess, maxfev=100
            )
        except RuntimeError:
            # print("Curve fit failure. Possibly weak binder.")
            fitting_variables, _ = curve_fit(
                fit, XM, exp_dQ_normalized[skip:], initial_guess, maxfev=10000
            )
            pass
        dH = fitting_variables[0]
        K = fitting_variables[1]
        N = fitting_variables[2]

        # (Print) Original Data, Fit, and find SumSqr
        fitdQ = fit(XM, dH, K, N)

        SumSqr = 0.0
        for i in range(skip, len(dQ)):
            SumSqr += (exp_dQ_normalized[i] - fitdQ[i - skip]) ** 2

        bootstrap_dH[cycle] = dH
        bootstrap_K[cycle] = K
        bootstrap_N[cycle] = N
        bootstrap_SS[cycle] = SumSqr
        bootstrap_dQ[:, cycle] = fitdQ

    # Reject all values which are over a threshold; definition is arbitrary.
    threshold = 1 * bootstrap_SS.mean() + 7.0 * np.sqrt(bootstrap_SS.std())
    for cycle in range(cycles - 1, -1, -1):
        if bootstrap_SS[cycle] > threshold:
            bootstrap_dH = np.delete(bootstrap_dH, cycle)
            bootstrap_K = np.delete(bootstrap_K, cycle)
            bootstrap_N = np.delete(bootstrap_N, cycle)
            bootstrap_SS = np.delete(bootstrap_SS, cycle)
            bootstrap_dQ = np.delete(bootstrap_dQ, cycle, 1)
    # Only use number of bootstrapping cycles defined in the main program.
    bootstrap_dH = bootstrap_dH[:realcycles]
    bootstrap_K = bootstrap_K[:realcycles]
    bootstrap_N = bootstrap_N[:realcycles]
    bootstrap_SS = bootstrap_SS[:realcycles]
    bootstrap_dQ = bootstrap_dQ[:realcycles]

    # Shortcut to get the exact uncertainty in delta G assuming both delta G and K are well-behaved Gaussians.
    # This could be resampled.
    R = 1.987_203_6 * 10 ** -3  # kcal K^-1 mol^-1
    dG_sem = R * temperature * np.std(bootstrap_K) / np.mean(bootstrap_K)
    dG = -R * temperature * np.log(np.mean(bootstrap_K))

    return bootstrap_dQ, XM, bootstrap_dH, bootstrap_K, bootstrap_N, dG, dG_sem


def report(
    V0,
    M0,
    X0,
    syringe_error,
    cell_error,
    heat_error,
    base_error,
    vardQ,
    dH,
    K,
    N,
    dG,
    dG_sem,
    temperature,
):
    print(f"{'Cell volume = ':<20} {V0:5.7f} L")
    print(f"{'Cell conc. = ':<20} {M0:5.7f} M")
    print(f"{'Injectant conc. = ':<20} {X0:5.7f} M")
    print(f"{'Syringe error = ':<20} {syringe_error * 100:4.2f} percent")
    print(f"{'Cell error = ':<20} {cell_error * 100 :4.2f} percent")
    print(f"{'Heat error = ':<20} {heat_error * 100:4.2f} percent")
    print(f"{'Base error = ':<20} {base_error * 1e6:4.2f} ucal")

    print(
        f"{'# dH = ':<20} {np.mean(dH) / 1000.:>10.2f} +/- {np.std(dH) / 1000.:>10.5f} kcal/mol"
    )
    print(
        f"{'# dH = ':<20} {np.mean(dH) /1000. * 4.184:>10.2f} +/- {np.std(dH)/ 1000. * 4.184:>10.5f} kJ/mol"
    )
    print(f"{'# K = ':<20} {np.mean(K):>10.2f} +/- {np.std(K):>10.5f} M^{{-1}}")
    print(f"{'# N = ':<20} {np.mean(N):>10.2f} +/- {np.std(N):>10.5f}")
    # print(f"{'# SS = ':<20} {np.mean(SS):>10.2f} +/- {np.std(SS):>10.5f}")

    print(f"{'# dG = ':<20} {dG:>10.2f} +/- {dG_sem:>10.5f} kcal/mol")
    print(f"{'# dG = ':<20} {dG * 4.184:>10.2f} +/- {dG_sem * 4.184:>10.5f} kJ/mol")

    dS = -1 * (dG - np.mean(dH) / 1000.0) / temperature
    dS_sem = np.sqrt(dG_sem ** 2 + (np.std(dH) / 1000.0) ** 2) / temperature

    print(f"{'# dS = ':<20} {dS * 1000:>10.5f} +/- {dS_sem * 1000:>10.5f} cal/mol/K")
    print(
        f"{'# dS = ':<20} {dS * 4.184 * 1000:>10.5f} +/- {dS_sem * 4.184 * 1000:>10.5f} J/mol/K"
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--file", help="Raw ITC data", required=True)
    ap.add_argument(
        "-s",
        "--skip",
        help="Number of injections to skip (default: 0)",
        default=0,
        required=False,
    )
    ap.add_argument(
        "-t",
        "--temperature",
        help="Temperature of the experiment (default: 300.15 K)",
        default=300.15,
        required=False,
    )
    ap.add_argument(
        "-M",
        "--M0",
        help="Cell concentration in Molar (default: 0.005 M)",
        default=0.005,
        required=False,
    )
    ap.add_argument(
        "-X",
        "--X0",
        help="Injectant concentration in Molar (default: 0.075 M)",
        default=0.075,
        required=False,
    )
    ap.add_argument(
        "-V",
        "--V0",
        help="Volume of the ITC cell in liters (default: 0.000202 L)",
        default=0.202 / 1000,
        required=False,
    )
    args = vars(ap.parse_args())

    temperature = float(args["temperature"])
    skip = int(args["skip"])
    X0 = float(args["X0"])
    M0 = float(args["M0"])
    V0 = float(args["V0"])

    # Load heat (dQ) and injection volumes (dV)
    dQ, dV = process(args["file"])

    # dQ and XM include the complete data from the ITC output file.
    # We refer to dQ as the heat release measured by the ITC and dH as the reaction enthalpy obtained from the fitting
    # process
    vardQ, XM, dH, K, N, dG, dG_sem = bootstrap(dQ, dV, temperature, cycles=1000)

    report(
        V0,
        M0,
        X0,
        syringe_error,
        cell_error,
        heat_error,
        base_error,
        vardQ,
        dH,
        K,
        N,
        dG,
        dG_sem,
        temperature,
    )

    ITC = dQ / (dV * X0)

    plot(XM, ITC, vardQ, args["file"], np.mean(dH), np.mean(K), np.mean(N), dG)

