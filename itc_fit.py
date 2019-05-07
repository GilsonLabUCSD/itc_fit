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

    return heat, volume


def plot(XM, ITC, skip, vardQ, name, dH, K, N, dG):
    """
    Plot the heat as a function of molar ratio, as measured, and the final (?) curve fit.
    I think this takes the best fit from the last bootstrapped cycle, but I'm not sure.
    """

    print(XM[0, 1:] / XM[1, 1:])

    fig, ax = plt.subplots(1, figsize=(6 * 1.2, 6))
    ax.scatter(XM[0, 1:] / XM[1, 1:], ITC[skip:], c="k", label="ITC data")
    for index, dQ in enumerate(vardQ.T):
        ax.errorbar(
            XM[0, 1:] / XM[1, 1:],
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
    dQ = np.zeros(injections - 1)
    X = XMa[0, :]
    M = XMa[1, :]
    for i in range(1, injections):  # Loop over injections
        # Total Heat
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
        dQ[i - 1] = (Q[i] + (dV[i - 1] / V0) * ((Q[i] + Q[i - 1]) / 2) - Q[i - 1]) / (
            dV[i - 1] * X0
        )

    return dQ


def bootstrap(dQ, dV, V0, skip, temperature, cycles=100):
    """
    Bootstrap the fitting of dQ given re-sampled uncertainties.
    """
    # Increase the number of bootstrapping runs by 10 percent. We delete the bad fits afterwards and want to fill up
    # the array so that we always use the same number of bootstrapping cycles
    realcycles = cycles
    cycles = int(cycles * 1.1)
    varITC = np.zeros([len(dQ) - skip])
    varheat = np.zeros([len(dQ) - skip, cycles])
    vardH = np.zeros([cycles])
    varK = np.zeros([cycles])
    varN = np.zeros([cycles])
    varSS = np.zeros([cycles])
    vardQ = np.zeros([len(dQ) - skip, cycles])
    for n in tqdm(range(cycles)):
        # Initial Guesses
        p0 = np.zeros(3)  # Guess Array
        p0[0] = -1000.0  # Guess dH
        p0[1] = 1000  # Guess K
        p0[2] = 1.000  # Guess N

        # Picking concentrations from gaussian distribution
        sampled_syringe_concentration = np.random.normal(X0, abs(syringe_error * X0))
        sampled_cell_concentration = np.random.normal(M0, abs(cell_error * M0))

        # Find concentrations after each injection
        XM = np.zeros([2, len(dQ) + 1 - skip])

        print(f"Number of injections to fit = {(len(dQ) + 1 - skip)}")

        XM[0, 0] = 0
        XM[1, 0] = 0

        cumulative_volume = np.cumsum(dV)
        cumulative_volume = cumulative_volume[skip:]

        print(np.shape(cumulative_volume))
        print(np.shape(XM[0, 1]))
        # New Injectant Concentration
        XM[0, 1:] = (cumulative_volume * sampled_syringe_concentration / V0) * (
            1 / (1 + (cumulative_volume / (2 * V0)))
        )

        # New Cell Molecule Concentration
        XM[1, 1:] = sampled_cell_concentration * (
            (1 - cumulative_volume / (2 * V0)) / (1 + cumulative_volume / (2 * V0))
        )

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

        print(varITC)
        print(varITC[skip:])

        # Fit the data
        try:
            popt, _ = curve_fit(fit, XM, varITC[skip:], p0, maxfev=100)
        except RuntimeError:
            print("Curve fit failure. Possibly weak binder.")
            popt, _ = curve_fit(fit, XM, varITC[skip:], p0, maxfev=10000)
            pass
        dH = popt[0]
        K = popt[1]
        N = popt[2]

        # (Print) Original Data, Fit, and find SumSqr

        fitdQ = fit(XM, dH, K, N)

        SumSqr = 0.0
        # for i in range(len(dQ)):
        #     print(np.shape(varITC))
        #     print(np.shape(fitdQ))
        #     print(np.shape(varheat))
        #     SumSqr += (varITC[i] - fitdQ[i]) ** 2
        #     varheat[i, n] = fitdQ[i]

        vardH[n] = dH
        varK[n] = K
        varN[n] = N
        # varSS[n] = SumSqr
        vardQ[:, n] = fitdQ

    # Reject all values which are over a threshold; definition is arbitrary.
    threshold = 1 * varSS.mean() + 5.0 * np.sqrt(varSS.std())
    for n in range(cycles - 1, -1, -1):
        if varSS[n] > threshold:
            # print('Rejecting cycle number {}'.format(n))
            vardH = np.delete(vardH, n)
            varK = np.delete(varK, n)
            varN = np.delete(varN, n)
            varSS = np.delete(varSS, n)
            vardQ = np.delete(vardQ, n, 1)
    # Only use number of bootstrapping cycles defined in the main program.
    # print(realcycles)
    vardH = vardH[:realcycles]
    varK = varK[:realcycles]
    varN = varN[:realcycles]
    varSS = varSS[:realcycles]
    vardQ = vardQ[:realcycles]

    # Shortcut to get the exact uncertainty in delta G assuming both delta G and K are well-behaved Gaussians.
    # This could be resampled.
    R = 1.9872036 * 10 ** -3  # kcal K^-1 mol^-1
    dG_sem = R * temperature * np.std(varK) / np.mean(varK)
    dG = -R * temperature * np.log(np.mean(varK))

    return vardQ, XM, vardH, varK, varN, dG, dG_sem


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

    vardQ, XM, dH, K, N, dG, dG_sem = bootstrap(dQ, dV, V0, skip, temperature, cycles=1)
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
    plot(XM, ITC, skip, vardQ, args["file"], np.mean(dH), np.mean(K), np.mean(N), dG)

