# =========================
# SamplingRVs_v2.py
#  - changes to overall plotting structure, lots of options added here
# =========================

def runsamples(
    runs, orbsamplerange, limit_on_samples, run_number,
    dispCDF=False, dispElse=False, elseerror=False, saveElse=False,
    SaveThisFig=False, PrintTestResults=False, savefig_name="Test",
    show=True, ax=None,
    # -------------------------
    # Plotting controls
    # -------------------------
    plot_obs=True,                 # master toggle: draw observation CDF/errorbars?
    plot_obs_full=True,            # show the "3+ measurements" observation set
    plot_obs_samplim=True,         # show the "{samplim}+ measurements" observation set
    plot_true=True,                # draw "actual" line(s)?
    plot_measured=True,            # draw "measured" line(s)?
    plot_fill=True,                # fill between actual and measured?
    sim_alpha=0.18,                # transparency for simulated overlays
    label_once=True,               # only label sim curves on run_number==1

    # -------------------------
    # label toggles
    # -------------------------
    label_sim=True,                # label simulated curves?
    label_obs=True,                # label observation series?

    # step vs line toggle for simulated curves
    sim_plot_style="step",         # "step" or "line"
    step_where="mid",              # where for step plots ("pre","post","mid")

    # -------------------------
    # observation grouping
    # -------------------------
    obs_group="all",               # "all", "ssg", "rs", "multi", "ssg_and_rs"

    # -------------------------
    # special observed stars overlay
    # -------------------------
    plot_obs_spec=False,           # toggle: overplot only selected stars?
    spec_stars=None,               # list of star names (strings) to overlay
    spec_label="Selected stars",   # legend label
    spec_style=None,               # optional dict to override style
    warn_missing_spec=True,        # print any requested names not found

    # -------------------------
    # return simulated distributions for aggregation
    # -------------------------
    return_sim_data=False,
    return_else_data=False,

    # -------------------------
    # external observation data input
    # -------------------------
    obs_data=None,                 # if provided, use this instead of the hard-coded 'data' list

    # -------------------------
    # label star names next to obs points (only meaningful for obs_group="all")
    # -------------------------
    label_obs_names=False,         # annotate star names next to their plotted obs point(s)
    obs_name_style=None,           # optional dict for text style (fontsize, alpha, etc.)
    obs_name_xytext=(4, 2),        # offset in points for label placement
    obs_name_subset=None,          # list of star names to label (None -> label all when label_obs_names=True)
):
    lower_samp = orbsamplerange[0]  # define upper and lower sampling values
    upper_samp = orbsamplerange[1]

    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from scipy.optimize import curve_fit

    import astropy
    import random
    from scipy.optimize import fsolve

    #####################################################
    # Generate Tidally synchronized distributions

    file = 'table1.dat'   # This is the RSCVn data table
    import glob
    folder_path = 'C:/Users/Jonah/Astro/Stars/DATA/' + file  # Change this to the path of the folder with table1.dat

    def is_float(string):
        # True if given string is float else False
        try:
            return float(string)
        except ValueError:
            return False

    stardata = []
    with open(folder_path, 'r') as f:
        d = f.readlines()
        for i in d:
            k = i.rstrip().split("     ")
            stardata.append([float(i) if is_float(i) else i for i in k])
    stardata = np.array(stardata, dtype='O')

    # Read a plain text .dat file
    data_txt = []
    with open(folder_path, 'r') as file:
        for line in file:
            columns = line.split()  # Splits by any whitespace
            if len(columns) == 18:
                combined_element = columns[1] + ' ' + columns[2]
                new_array = [columns[0], combined_element] + columns[3:]
                data_txt.append(new_array)
            else:
                data_txt.append(columns)

    periods = []
    RSsP = []
    SSGsP = []
    extra = []
    test = []
    extraind = []
    extraid = []

    TessRS = []
    TessSSGs = []
    Tessperiods = []
    for i in range(len(data_txt)):

        if data_txt[i][16] == 'RS':
            RSsP.append(data_txt[i][11])
            periods.append(data_txt[i][11])
            if float(data_txt[i][12]) != -999:
                TessRS.append(data_txt[i][12])
                Tessperiods.append(data_txt[i][12])
        elif data_txt[i][16] == 'SSG':
            SSGsP.append(data_txt[i][11])
            test.append(data_txt[i][16])
            periods.append(data_txt[i][11])
            if float(data_txt[i][12]) != -999:
                TessSSGs.append(data_txt[i][12])
                Tessperiods.append(data_txt[i][12])
        elif data_txt[i][16] != 'cut':
            extra.append(data_txt[i][16])
            extraind.append(i)
            extraid.append(data_txt[i][len(data_txt[i]) - 1])

    periods = [float(item) for item in periods]
    SSGsP = [float(item) for item in SSGsP]
    RSsP = [float(item) for item in RSsP]

    Tessperiods = [float(item) for item in Tessperiods]
    TessSSGs = [float(item) for item in TessSSGs]
    TessRS = [float(item) for item in TessRS]

    SSGperiods = SSGsP
    Tessperiods = TessSSGs

    #overwrite SSGperiods as RSsP here - may be able to adjust/change if we want SSG exclusive 
    #overwritten to provide a larger sample size for the "tidally locked" distribution
    #We aren't looking at just the SSGs here (as we reclassify plenty of them)
    #we want 'tidally locked' distribution of stars we assume are binaries, so including all RS in this tidally locked sample is valid.
    SSGperiods = RSsP

    randPvals = []
    for i in range(1000):
        hist_counts, bins = np.histogram(SSGperiods, bins=25, density=True)
        pdfdata = hist_counts / sum(hist_counts)
        cdfdata = np.cumsum(pdfdata)
        cdfdata = np.insert(cdfdata, 0, 0)
        uniform_random_value = np.random.rand()
        PorbRand = np.interp(uniform_random_value, cdfdata, bins)
        randPvals.append(PorbRand)

    #####################################################
    # Begin CDF Generation

    totcorr = []
    numsamp = []
    logPvals, logPfullvals, logP1000vals, logPtidal = [], [], [], []
    inclinations = []
    randPvals = []

    actual1, actual2, actual3 = [], [], []
    measured1, measured2, measured3 = [], [], []
    periods1, periods2, periods3 = [], [], []
    Kvals1, Kvals2, Kvals3 = [], [], []

    actual1e, actual2e, actual3e = [], [], []
    measured1e, measured2e, measured3e = [], [], []
    periods1e, periods2e, periods3e = [], [], []
    Kvals1e, Kvals2e, Kvals3e = [], [], []
    eccens, eccens_logPvals = [], []
    ec_all = []

    actualinfe, measuredinfe, periodsinfe, Kvalsinfe = [], [], [], []

    from tqdm import tqdm

    for i in tqdm(range(runs), desc=f"Run {run_number}"):

        Msol = 1.9891 * 10**30
        G = 6.67430e-11

        # Lognormal Distribution
        ux = 4.8
        ox = 2.3
        mode = 4.8

        u = np.log((ux**2) / np.sqrt(ux**2 + ox**2))
        o = np.sqrt(np.log(1 + (ox**2 / ux**2)))
        um = np.log(mode) + o**2

        m1 = random.uniform(0.9, 2)
        m2 = random.uniform(0.1, m1)

        def lognormal(mean, sigma, range):
            while True:
                logP = np.random.lognormal(mean, sigma)
                if logP > np.log10(2) and logP < range:
                    break
            return logP

        while True:
            logP = lognormal(um, o, 1.5)
            logvalue = 10**logP
            if logvalue <= 30 and logvalue >= 2:  # 30
                break

        while True:
            logP1000 = lognormal(um, o, 3)
            logvalue1000 = 10**logP1000
            if logvalue1000 <= 1000 and logvalue1000 >= 2:
                break

        while True:
            logPfull = lognormal(um, o, 10)
            logvaluefull = 10**logPfull
            if logvaluefull >= 2:
                break

        # Tidally Locked Distribution
        hist_counts, bins = np.histogram(SSGperiods, bins=25, density=True)
        pdfdata = hist_counts / sum(hist_counts)
        cdfdata = np.cumsum(pdfdata)
        cdfdata = np.insert(cdfdata, 0, 0)

        while True:
            uniform_random_value = np.random.rand()
            PorbRand = np.interp(uniform_random_value, cdfdata, bins)
            if PorbRand <= 30:
                break

        logtidal = np.log10(PorbRand)

        logPvals.append(logP)
        logPfullvals.append(logPfull)
        logP1000vals.append(logP1000)
        logPtidal.append(logtidal)

        # Simulate Mass Ratios, Inclinations
        M1 = m1 * Msol
        M2 = m2 * Msol
        v = random.uniform(0, 1)
        incl = np.arccos(2 * v - 1)
        inclinations.append(incl)

        # Decide on Period
        b1 = 10**logP
        binf = 10**logPfull
        b2 = 10**logP1000
        b3 = 10**logtidal

        # Generate Orbit time grid
        x = np.linspace(0, 365*2, 10000)

        def genorbit(b):
            Porb = b * 86400
            K = (((M2**3 * np.sin(incl)**3 * (2 * np.pi * G)) / (Porb * (M1 + M2)**2)))**(1/3)
            a = K/1000
            c = np.random.random() * 10
            equation = a * np.sin(((2*np.pi)/b) * x + c)
            return equation, a, b, 0

        def genorbit_e(Period):
            P = Period
            Porb = Period * 86400

            if np.log10(P) < 1:
                e = 0
            elif np.log10(P) >= 1 and np.log10(P) < 2:
                e = random.uniform(0, 0.4)
            elif np.log10(P) >= 2 and np.log10(P) < 3:
                e = random.uniform(0.2, 0.6)
                if elseerror:
                    e = e + random.uniform(-0.19, 0.19)
            elif np.log10(P) >= 3:
                e = random.uniform(0.2, 0.8)
                if elseerror:
                    e = e + random.uniform(-0.19, 0.19)

            eccens_logPvals.append(np.log10(Period))

            Ke = (((M2**3 * np.sin(incl)**3 * (2 * np.pi * G)) / (Porb * (M1 + M2)**2 * (1-e**2)**3/2)))**(1/3)

            omega = np.deg2rad(np.random.uniform(0, 360))
            t = x
            n = 2 * np.pi / P
            M = n * t

            E = M.copy()

            def newton_solve_kepler(E, M, e, tol=1e-6, max_iter=100):
                for _ in range(max_iter):
                    f_E = E - e * np.sin(E) - M
                    f_prime_E = 1 - e * np.cos(E)
                    delta_E = -f_E / f_prime_E
                    E += delta_E
                    if np.all(np.abs(delta_E) < tol):
                        break
                return E

            E = newton_solve_kepler(E, M, e)

            nu = 2 * np.arctan2(
                np.sqrt(1 + e) * np.sin(E / 2),
                np.sqrt(1 - e) * np.cos(E / 2)
            )

            v_r = Ke * (np.cos(omega + nu) + e * np.cos(omega))
            return v_r/1000, Ke/1000, P, e

        orbit3, a3, b3, e3 = genorbit(b3)
        orbit1e, a1e, b1e, e1e = genorbit_e(b1)
        orbit2e, a2e, b2e, e2e = genorbit_e(b2)
        orbitinfe, ainfe, binfe, einfe = genorbit_e(binf)

        ec_all.append([einfe, e2e, e1e, e3])

        samples = random.randint(lower_samp, upper_samp)

        def sample(orbit, samples):
            xval = []
            yval = []
            obsday = 0
            prevsamp = int((np.random.random()*10000) - 1)
            for k in range(samples):
                obsday = obsday + 1
                if obsday < int(np.random.random()*2) + 2:
                    samp = prevsamp + int(np.random.random()*75)
                    if samp > 9999:
                        samp = 9999
                    prevsamp = samp
                    yval.append(orbit[samp])
                    xval.append(x[samp])
                else:
                    samp = int((np.random.random()*10000))
                    if samp > 9999:
                        samp = 9999
                    prevsamp = samp
                    yval.append(orbit[samp])
                    xval.append(x[samp])
                    obsday = 0
            return xval, yval

        xval3, yval3 = sample(orbit3, samples)
        xval1e, yval1e = sample(orbit1e, samples)
        xval2e, yval2e = sample(orbit2e, samples)
        xvalinfe, yvalinfe = sample(orbitinfe, samples)

        def measure(yval, a):
            std = np.std(yval)
            Amp = (np.max(yval) - np.min(yval)) / 2
            return Amp, std

        Amp3, std3 = measure(yval3, a3)
        Amp1e, std1e = measure(yval1e, a1e)
        Amp2e, std2e = measure(yval2e, a2e)
        Ampinfe, stdinfe = measure(yvalinfe, ainfe)

        actual3.append(a3)
        actual1e.append(a1e)
        actual2e.append(a2e)
        actualinfe.append(ainfe)

        measured3.append(Amp3)
        measured1e.append(Amp1e)
        measured2e.append(Amp2e)
        measuredinfe.append(Ampinfe)

        periods3.append(b3)
        periods1e.append(b1e)
        periods2e.append(b2e)
        periodsinfe.append(binf)

        Kvals3.append(a3)
        Kvals1e.append(a1e)
        Kvals2e.append(a2e)
        Kvalsinfe.append(ainfe)

    # Now, generate CDF plots
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import rankdata

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # Helpers for step vs line
    sim_plot_style = str(sim_plot_style).lower().strip()
    if sim_plot_style not in ("step", "line"):
        raise ValueError(f"sim_plot_style must be 'step' or 'line' (got {sim_plot_style!r})")

    def draw_sim_curve(xv, yv, *, linestyle, color, alpha, label):
        if sim_plot_style == "step":
            ax.step(xv, yv, where=step_where, linestyle=linestyle, color=color, alpha=alpha, label=label)
        else:
            ax.plot(xv, yv, linestyle=linestyle, color=color, alpha=alpha, label=label)

    # --- SIMULATION overlays ---
    for actual, measured, label, color in zip(
        [actualinfe, actual2e, actual1e, actual3],
        [measuredinfe, measured2e, measured1e, measured3],
        [r'Field Binary Distribution',
         'Truncated Field Distribution',
         'Short Period, Log-normal Distribution',
         'Synchronized SSG Distribution'],
        ['blue', 'green', 'red', 'y']
    ):
        actual_sorted = np.sort(actual)
        actual_cdf = rankdata(actual_sorted, method='average') / len(actual_sorted)

        measured_sorted = np.sort(measured)
        measured_cdf = rankdata(measured_sorted, method='average') / len(measured_sorted)

        # force curves to start at (0,0)
        actual_sorted = np.insert(actual_sorted, 0, 0.0)
        actual_cdf    = np.insert(actual_cdf,    0, 0.0)

        measured_sorted = np.insert(measured_sorted, 0, 0.0)
        measured_cdf    = np.insert(measured_cdf,    0, 0.0)

        if not dispCDF:
            continue

        # NEW: label toggle applied here
        do_label = label_sim and ((not label_once) or (run_number == 1))
        true_lbl = f"{label} (true)" if do_label else None
        meas_lbl = f"{label} (measured)" if do_label else None
        fill_lbl = f"{label} (bias band)" if do_label else None

        if plot_true:
            draw_sim_curve(
                actual_sorted, actual_cdf * 100,
                linestyle='-', color=color, alpha=sim_alpha, label=true_lbl
            )

        if plot_measured:
            draw_sim_curve(
                measured_sorted, measured_cdf * 100,
                linestyle='-', color=color, alpha=sim_alpha, label=meas_lbl
            )

        if plot_fill:
            ax.fill_betweenx(
                actual_cdf * 100,
                actual_sorted,
                measured_sorted,
                color=color,
                alpha=sim_alpha * 0.6,
                label=fill_lbl
            )

    ##################
    nan = np.nan
    # Star, N, RV1_amp, rv1_hi, rv1_lo, RV2_amp, rv2_hi, rv2_lo
    if obs_data is not None:
        data = obs_data
    else:
        data = [
            ['29Dra', 5, 0.992, 0.151, 0.167, nan, nan, nan, 'SSG', 'SB1'],
            ['2MASSJ08171221+0736164', 4, 12.751, 0.193, 0.153, nan, nan, nan, 'RS', 'SB1'],
            #['2MASSJ08504952+1217158', 4, 18.098, 0.324, 0.380, 4.965, 0.369, 0.565, '--', 'SB2'],
            ['2MASSJ17270871+2700144', 3, 3.493, 0.051, 0.045, nan, nan, nan, 'SSG', 'SB1'],
            ['54Cam', 5, 49.191, 0.157, 0.155, 48.586, 0.248, 0.225, 'RS', 'SB2'],
            ['ASASJ073642+0354.3', 8, 0.576, 0.108, 0.091, nan, nan, nan, 'SSG', 'SB1'],
            ['ASASJ155307+2028.6', 4, 47.053, 0.211, 0.258, nan, nan, nan, 'SSG', 'SB1'],
            ['ASASJ162510+0514.9', 9, 55.507, 0.970, 1.043, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['ASASJ210323+1930.9', 12, 28.870, 0.344, 0.367, nan, nan, nan, 'SSG', 'SB1'],
            ['BD+142519', 7, 34.534, 0.264, 0.272, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['BD+201895', 3, 1.123, 0.197, 0.210, nan, nan, nan, 'RS', 'Asymmetric'],
            ['BMCVn', 7, 12.139, 0.136, 0.153, nan, nan, nan, 'SSG', 'SB1'],
            ['BPSCS22950-0027', 12, 26.011, 0.153, 0.150, nan, nan, nan, 'SSG', 'SB1'],
            ['CGTri', 12, 39.050, 0.654, 0.553, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['EPCVn', 12, 49.956, 0.106, 0.117, nan, nan, nan, 'SSG', 'SB1'],
            ['FMLyn', 6, 31.651, 0.287, 0.290, nan, nan, nan, 'RS', 'SB1'],
            ['FNCom', 13, 43.760, 0.303, 0.297, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['FRBoo', 9, 3.316, 0.161, 0.156, nan, nan, nan, 'RS', 'SB1'],
            #['FVCnc', 4, 27.034, 0.077, 0.069, 28.581, 0.118, 0.126, '--', 'SB2'],
            ['GSC00563-00384', 10, 14.422, 0.107, 0.125, nan, nan, nan, 'SSG', 'SB1'],
            ['GSC02769-00149', 11, 23.567, 0.323, 0.348, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['HD113714', 5, 74.895, 0.248, 0.246, nan, nan, nan, 'RS', 'SB1'],
            ['HD161570', 6, 21.266, 0.061, 0.062, nan, nan, nan, 'RS', 'SB1'],
            ['HD347929', 13, 0.853, 0.132, 0.142, nan, nan, nan, 'SSG', 'SB1'],
            ['HD53476', 4, 20.336, 0.152, 0.135, nan, nan, nan, 'RS', 'SB1'],
            ['HRPsc', 12, 21.484, 0.707, 0.667, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['HUVir', 12, 46.728, 0.691, 0.683, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['KLCnc', 5, 35.939, 0.845, 0.897, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['KNCnc', 9, 27.324, 0.153, 0.160, 6.229, 0.319, 0.273, 'SSG', 'SB2'],
            ['KYBoo', 10, 0.510, 0.184, 0.205, nan, nan, nan, 'SSG', 'SB1'],
            ['MPBoo', 5, 0.276, 0.079, 0.079, nan, nan, nan, 'RS', 'SB1'],
            ['MSBoo', 5, 18.595, 0.140, 0.199, nan, nan, nan, 'SSG', 'SB1'],
            ['MXBoo', 4, 0.262, 0.108, 0.111, nan, nan, nan, 'RS', 'SB1'],
            ['MXDra', 9, 39.695, 0.266, 0.305, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['NPCom', 6, 61.883, 0.431, 0.413, nan, nan, nan, 'SSG', 'SB1'],
            ['NRLib', 8, 50.539, 1.143, 1.246, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['OOBoo', 4, 30.136, 0.137, 0.117, 22.988, 0.377, 0.340, 'RS', 'SB2'],
            ['QZBoo', 12, 32.924, 0.146, 0.190, nan, nan, nan, 'SSG', 'SB1'],
            ['ROTSE1J135649.21+242927.1', 6, 23.842, 0.109, 0.092, 2.318, 0.113, 0.118, 'RS', 'SB2'],
            ['ROTSE1J160518.19+372623.5', 4, 46.531, 0.138, 0.132, 45.682, 0.348, 0.369, 'SSG', 'SB2'],
            ['ROTSE1J162157.22+381733.6', 4, 27.572, 0.103, 0.099, nan, nan, nan, 'RS', 'SB1'],
            ['ROTSE1J170537.88+335052.5', 4, 20.952, 0.274, 0.282, 16.809, 1.050, 0.991, 'RS', 'SB2'],
            ['ROTSE1J172404.95+182937.3', 3, 18.915, 0.155, 0.160, nan, nan, nan, 'SSG', 'SB1'],
            ['ROTSE1J183703.58+280017.1', 5, 12.077, 0.226, 0.246, nan, nan, nan, 'RS', 'Asymmetric'],
            ['ROTSE1J184706.07+434032.7', 5, 7.540, 0.127, 0.137, nan, nan, nan, 'SSG', 'SB1'],
            ['RXJ1112.4+3950', 15, 51.708, 0.774, 0.626, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['TYC1339-572-1', 4, 12.046, 0.169, 0.178, 8.673, 0.131, 0.121, 'SSG', 'SB2'],
            ['TYC1543-131-1', 3, 7.990, 0.086, 0.088, nan, nan, nan, 'RS', 'SB1'],
            ['TYC2268-394-1', 11, 58.238, 1.174, 1.262, nan, nan, nan, 'SSG', 'SB1'],
            ['TYC2645-156-1', 6, 25.825, 0.236, 0.250, nan, nan, nan, 'RS', 'Asymmetric'],
            #['TYC3401-730-1', 3, 62.355, 1.527, 1.685, nan, nan, nan, '--', 'SB1'],
            ['TYC3557-919-1', 10, 26.490, 0.129, 0.127, nan, nan, nan, 'SSG', 'SB1'],
            ['TYC3927-1130-1', 5, 35.127, 0.386, 0.328, nan, nan, nan, 'RS', 'SB1'],
            ['TYC4216-800-1', 5, 15.280, 0.145, 0.116, nan, nan, nan, 'RS', 'SB1'],
            ['TYC4223-399-1', 23, 26.459, 1.976, 2.145, nan, nan, nan, 'SSG', 'SB1'],
            ['TYC4415-328-1', 10, 12.916, 0.154, 0.212, nan, nan, nan, 'SSG', 'SB1'],
            ['TYC4422-1970-1', 18, 38.306, 0.552, 0.511, nan, nan, nan, 'RS', 'SB1'],
            ['TYC4428-785-1', 3, 9.378, 0.106, 0.091, nan, nan, nan, 'RS', 'SB1'],
            ['TYC4667-90-1', 11, 52.485, 0.344, 0.450, nan, nan, nan, 'SSG', 'Asymmetric'],
            ['TYC5003-138-1', 8, 46.438, 0.246, 0.288, 60.411, 0.286, 0.297, 'SSG', 'SB2'],
            ['TYC583-566-1', 11, 51.934, 0.516, 0.503, nan, nan, nan, 'SSG', 'SB1'],
            ['UCAC228842571', 35, 45.614, 0.509, 1.075, 46.522, 7.848, 17.727, 'SSG', 'SB2'],
            ['V424Gem', 4, 28.076, 0.091, 0.106, nan, nan, nan, 'SSG', 'SB1'],
            ['V457Vir', 6, 47.867, 0.446, 0.405, nan, nan, nan, 'RS', 'Asymmetric'],
            ['V498And', 11, 30.630, 0.280, 0.373, nan, nan, nan, 'SSG', 'SB1'],
            #['V503Hya', 5, 70.636, 0.647, 0.578, 54.455, 0.904, 0.902, '--', 'SB2'],     #
            #['V832Her', 5, 10.286, 0.059, 0.058, nan, nan, nan, '--', 'SB1'],            #
            #['V834Her', 7, 1.159, 0.165, 0.194, nan, nan, nan, '--', 'SB1'],             #
            #['V846Her', 5, 17.532, 0.206, 0.262, nan, nan, nan, '--', 'SB1'],            #
        ]

    # ----------------------------
    # NEW: obs_group filtering (NO OTHER FUNCTIONALITY CHANGED)
    # ----------------------------
    obs_group = str(obs_group).lower().strip()
    if obs_group not in ("all", "ssg", "rs", "multi", "ssg_and_rs"):
        raise ValueError(
            "obs_group must be 'all', 'ssg', 'rs', 'multi', or 'ssg_and_rs' "
            f"(got {obs_group!r})"
        )

    def filter_rows(rows, grp):
        if grp == "all":
            return rows
        if grp == "ssg":
            return [r for r in rows if str(r[8]).strip().upper() == "SSG"]
        if grp == "rs":
            return [r for r in rows if str(r[8]).strip().upper() == "RS"]
        if grp == "ssg_and_rs":
            return [r for r in rows if str(r[8]).strip().upper() in ("SSG", "RS")]
        return rows

    # styles to differentiate overlays
    group_style = {
        "all": dict(fmt=".", color="black", ecolor="black", markersize=18, markeredgewidth=0,
                    elinewidth=1.5, capsize=2, capthick=1.2, alpha=1, zorder=200),
        "ssg": dict(fmt=".", color="salmon", ecolor="black", markersize=22,  markeredgewidth=1.5, markeredgecolor="black",
                    elinewidth=1.5, capsize=2, capthick=1.2, alpha=0.75, zorder=210),
        "rs":  dict(fmt=".", color="cyan", ecolor="black", markersize=22,  markeredgewidth=1.5, markeredgecolor="black",
                    elinewidth=1.5, capsize=2, capthick=1.2, alpha=0.75, zorder=205),
    }
    group_name = {"all": "All", "ssg": "SSGs", "rs": "RS CVn"}

    def compute_obs_arrays(rows):
        import numpy as np

        sampnums = [row[1] for row in rows]
        amps = [row[2] for row in rows]
        amps_hi = [row[3] for row in rows]
        amps_lo = [row[4] for row in rows]

        amps_hi_all = amps_hi.copy()
        amps_lo_all = amps_lo.copy()
        amps_hi_all += [row[6] for row in rows if not np.isnan(row[5])]
        amps_lo_all += [row[7] for row in rows if not np.isnan(row[5])]

        samplelim = limit_on_samples
        ampnew, hi_new, lo_new = [], [], []
        for s, a, hi, lo in zip(sampnums, amps, amps_hi_all, amps_lo_all):
            if s >= samplelim or a > 50:
                ampnew.append(a)
                hi_new.append(hi)
                lo_new.append(lo)

        vals1 = np.array(sorted(amps))
        lo1 = np.array([lo for _, lo in sorted(zip(amps, amps_lo_all))])
        hi1 = np.array([hi for _, hi in sorted(zip(amps, amps_hi_all))])
        p1 = np.arange(len(vals1)) / (len(vals1) - 1) if len(vals1) > 1 else np.array([1.0])

        vals2 = np.array(sorted(ampnew))
        lo2 = np.array([lo for _, lo in sorted(zip(ampnew, lo_new))]) if len(vals2) else np.array([])
        hi2 = np.array([hi for _, hi in sorted(zip(ampnew, hi_new))]) if len(vals2) else np.array([])
        p2 = (np.arange(len(vals2)) / (len(vals2) - 1)) if len(vals2) > 1 else (np.array([1.0]) if len(vals2) == 1 else np.array([]))

        return vals1, lo1, hi1, p1, vals2, lo2, hi2, p2, ampnew

    def y_from_cdf(x, xs, ps):
        import numpy as np
        xs = np.asarray(xs, dtype=float)
        ps = np.asarray(ps, dtype=float)
        if xs.size == 0:
            return np.nan
        x = float(x)
        mask = (xs == x)
        if np.any(mask):
            return float(np.mean(ps[mask]) * 100.0)
        if xs.size == 1:
            return float(ps[0] * 100.0)
        return float(np.interp(x, xs, ps) * 100.0)

    # ----------------------------
    # NEW: obs name annotation style defaults
    # ----------------------------
    default_name_style = dict(
        fontsize=7,
        alpha=0.85,
        ha="left",
        va="bottom",
        zorder=1000,
    )
    if isinstance(obs_name_style, dict):
        default_name_style.update(obs_name_style)

    # --- Observations (per-layer toggles) ---
    if dispCDF == True and plot_obs:

        if obs_group == "multi":
            groups_to_plot = ["all", "ssg", "rs"]
        elif obs_group == "ssg_and_rs":
            groups_to_plot = ["ssg", "rs"]
        else:
            groups_to_plot = [obs_group]

        # Track whether we should annotate names this run:
        annotate_names = bool(label_obs_names)

        # Build optional name filter (case-insensitive)
        subset_lc = None
        if obs_name_subset is not None:
            subset_lc = set([str(s).strip().lower() for s in obs_name_subset if str(s).strip()])

        def want_label(starname: str) -> bool:
            if subset_lc is None:
                return True
            return str(starname).strip().lower() in subset_lc

        for grp in groups_to_plot:
            data_g = filter_rows(data, grp)
            if len(data_g) == 0:
                continue

            vals1, lo1, hi1, p1, vals2, lo2, hi2, p2, ampnew = compute_obs_arrays(data_g)
            sty = group_style[grp]
            gname = group_name[grp]
            samplelim = limit_on_samples

            if plot_obs_full:
                ax.errorbar(
                    vals1, p1 * 100,
                    xerr=[lo1, hi1],
                    label=((f"{gname} observations") if label_obs else None),
                    **sty
                )

            if plot_obs_samplim and len(vals2) > 0:
                sty2 = dict(sty)
                sty2["alpha"] = 0.85
                ax.errorbar(
                    vals2, p2 * 100,
                    xerr=[lo2, hi2],
                    label=((f"{gname} observations with {samplelim}+ Measurements") if label_obs else None),
                    **sty2
                )

            # ----------------------------
            # CHANGED: annotate star names for the CURRENT group
            # (moved OUT of the samplim block so FULL labels work too)
            # ----------------------------
            if annotate_names:

                # FULL set labels
                if plot_obs_full:
                    for row in data_g:
                        star = str(row[0]).strip()
                        if not want_label(star):
                            continue
                        amp = float(row[2])
                        y = y_from_cdf(amp, vals1, p1)
                        if np.isfinite(y):
                            ax.annotate(
                                star,
                                xy=(amp, y),
                                xytext=obs_name_xytext,
                                textcoords="offset points",
                                **default_name_style
                            )

                # SAMPLIM set labels (only those that survive the ampnew criterion in THIS group)
                if plot_obs_samplim and len(vals2) > 0:
                    ampnew_set = set([float(a) for a in ampnew])
                    for row in data_g:
                        star = str(row[0]).strip()
                        if not want_label(star):
                            continue
                        amp = float(row[2])
                        if amp in ampnew_set:
                            y2 = y_from_cdf(amp, vals2, p2)
                            if np.isfinite(y2):
                                ax.annotate(
                                    star,
                                    xy=(amp, y2),
                                    xytext=obs_name_xytext,
                                    textcoords="offset points",
                                    **default_name_style
                                )

        # ----------------------------
        # Special stars overlay at their TRUE CDF positions
        # ----------------------------
        if plot_obs_spec and spec_stars:

            if obs_group in ("all", "ssg", "rs", "ssg_and_rs"):
                ref_grp = obs_group
            else:
                ref_grp = "all"

            ref_rows = filter_rows(data, ref_grp)

            if len(ref_rows) > 0:

                vals1_ref, lo1_ref, hi1_ref, p1_ref, vals2_ref, lo2_ref, hi2_ref, p2_ref, _ = compute_obs_arrays(ref_rows)

                wanted = [str(s).strip() for s in spec_stars if str(s).strip()]
                wanted_lc = [w.lower() for w in wanted]

                available = {str(r[0]).strip().lower(): r for r in data}

                data_spec = []
                found_lc = set()
                for nm_lc in wanted_lc:
                    if nm_lc in available:
                        data_spec.append(available[nm_lc])
                        found_lc.add(nm_lc)

                if warn_missing_spec:
                    missing = [w for w, wlc in zip(wanted, wanted_lc) if wlc not in found_lc]
                    if missing:
                        print("[runsamples] WARNING: spec_stars not found in data array:")
                        for m in missing:
                            print("  -", m)

                if len(data_spec) > 0:

                    default_spec_style = dict(
                        fmt="*",
                        color="magenta",
                        ecolor="black",
                        markersize=16,
                        markeredgewidth=1,
                        markeredgecolor="black",
                        elinewidth=1,
                        capsize=3,
                        capthick=1.4,
                        alpha=0.95,
                        zorder=400,
                    )
                    if isinstance(spec_style, dict):
                        default_spec_style.update(spec_style)

                    used_label = False
                    for row in data_spec:
                        amp = float(row[2])
                        hi = float(row[3])
                        lo = float(row[4])

                        if plot_obs_full:
                            y = y_from_cdf(amp, vals1_ref, p1_ref)
                            ax.errorbar(
                                [amp], [y],
                                xerr=[[lo], [hi]],
                                label=(spec_label if (label_obs and not used_label) else None),
                                **default_spec_style
                            )
                            used_label = True

                        if plot_obs_samplim and len(vals2_ref) > 0:
                            y2 = y_from_cdf(amp, vals2_ref, p2_ref)
                            sty2 = dict(default_spec_style)
                            sty2["alpha"] = min(1.0, default_spec_style.get("alpha", 1.0) * 0.9)
                            ax.errorbar(
                                [amp], [y2],
                                xerr=[[lo], [hi]],
                                label=((f"{spec_label} ({limit_on_samples}+)")
                                       if (label_obs and not used_label) else None),
                                **sty2
                            )
                            used_label = True

        ax.set_xlabel('K (km/s)')
        ax.set_ylabel('Probability (%)')
        ax.set_xlim(-5, 105)
        ax.set_ylim(-5, 105)
        ax.legend(loc='lower right', fontsize=10)
        plt.grid(alpha=0.5, linestyle="--")
        plt.tight_layout()

        if SaveThisFig:
            plt.savefig(savefig_name + ".png", dpi=300)

    # ----------------------------
    # KS testing (unchanged)
    # ----------------------------
    from scipy.stats import ks_2samp

    aKSval = []
    mKSval = []
    apval = []
    mpval = []

    def KS(amps_in, start):
        for actual, measured in zip(
            [actualinfe, actual2e, actual1e, actual3],
            [measuredinfe, measured2e, measured1e, measured3]
        ):
            ks_statistic, p_value = ks_2samp(amps_in, measured)
            mKSval.append(ks_statistic)
            mpval.append(p_value)

            ks_statistic, p_value = ks_2samp(amps_in, actual)
            aKSval.append(ks_statistic)
            apval.append(p_value)

    sampnums_all = [row[1] for row in data]
    amps_all = [row[2] for row in data]
    amps_hi_all0 = [row[3] for row in data]
    amps_lo_all0 = [row[4] for row in data]

    amps_hi_all2 = amps_hi_all0.copy()
    amps_lo_all2 = amps_lo_all0.copy()
    amps_hi_all2 += [row[6] for row in data if not np.isnan(row[5])]
    amps_lo_all2 += [row[7] for row in data if not np.isnan(row[5])]

    samplelim0 = limit_on_samples
    ampnew_all, hi_new_all, lo_new_all = [], [], []
    for s, a, hi, lo in zip(sampnums_all, amps_all, amps_hi_all2, amps_lo_all2):
        if s >= samplelim0 or a > 50:
            ampnew_all.append(a)
            hi_new_all.append(hi)
            lo_new_all.append(lo)

    KS(amps_all, "Full Distribution")
    KS(ampnew_all, f"Observations with {samplelim0}+ Measurements")

    if PrintTestResults == True:
        print("KS, measured sets:", [f"{val:.4f}" for val in mKSval])
        print("KS, actual sets:  ", [f"{val:.4f}" for val in aKSval])
        print()
        print("p-val, measured sets:", [f"{val:.4f}" for val in mpval])
        print("p-val, actual sets:  ", [f"{val:.4f}" for val in apval])

    # ----------------------------
    # NEW: package "else" data for pooling across runs (NO plotting here)
    # ----------------------------
    else_data = None
    if return_else_data:
        else_data = {
            "logPfullvals":  np.array(logPfullvals, dtype=float),
            "logP1000vals":  np.array(logP1000vals, dtype=float),
            "logPvals":      np.array(logPvals, dtype=float),
            "logPtidal":     np.array(logPtidal, dtype=float),
            "ec_all":        np.array(ec_all, dtype=float),  # shape (runs, 4)
        }

    sim_data = None
    if return_sim_data:
        sim_data = {
            "ecc_lognorm_2_1e10": {
                "actual":   np.array(actualinfe, dtype=float),
                "measured": np.array(measuredinfe, dtype=float),
            },
            "ecc_lognorm_2_1000": {
                "actual":   np.array(actual2e, dtype=float),
                "measured": np.array(measured2e, dtype=float),
            },
            "ecc_lognorm_2_30": {
                "actual":   np.array(actual1e, dtype=float),
                "measured": np.array(measured1e, dtype=float),
            },
            "tidal_locked": {
                "actual":   np.array(actual3, dtype=float),
                "measured": np.array(measured3, dtype=float),
            },
        }

        if return_else_data:
            return mKSval, aKSval, mpval, apval, sim_data, else_data

        return mKSval, aKSval, mpval, apval, sim_data

    if return_else_data:
        return mKSval, aKSval, mpval, apval, else_data

    return mKSval, aKSval, mpval, apval