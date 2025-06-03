

def runsamples(runs, orbsamplerange, limit_on_samples, run_number, dispCDF = False, SaveThisFig = False, PrintTestResults = False, savefig_name = "Test"):

    lower_samp = orbsamplerange[0]  #define upper and lower sampling values
    upper_samp = orbsamplerange[1]
    #print(lower_samp, upper_samp)
    
    #package imports
    #If using Jupyter locally, make sure you have pip installed saphires and any dependencies.
    
    #print(f"Run {run_number}")   #If doing many runs, uncomment this line

    
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from scipy.optimize import curve_fit   
    


    import astropy
    import random
    from scipy.optimize import fsolve
    
    #x = [0,1,2,3,4,5,6,7,8,9]
    #y = x
    #plt.plot(x, y)    #Testing out matplotlib inline to make sure this runs
    #print()
    #plt.show()

    #####################################################
    #Generate Tidally synchronized distributions
    
    file = 'table1.dat'   #This is the RSCVn data table
    import os
    import glob
    folder_path = 'C:/Users/Jonah/Astro/Stars/DATA/' + file  # Change this to the path of the folder with table1.dat
    
    #print(folder_path)
    
    #.dat file extraction, checks .p
    
    
    def is_float(string):
        #True if given string is float else False
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
    data = []
    with open(folder_path, 'r') as file:
        for line in file:
            columns = line.split()  # Splits by any whitespace
            #print(len(columns))
            if len(columns) == 18:
                combined_element = columns[1] + ' ' + columns[2]
                new_array = [columns[0], combined_element] + columns[3:]
                #print(new_array)
                #print(len(new_array))
                data.append(new_array)
                #print()
            else:
                #print(columns)
                data.append(columns)
                #print()


    
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
    for i in range(len(data)):
        
        if data[i][16] == 'RS':
            RSsP.append(data[i][11])
            periods.append(data[i][11])
            if float(data[i][12]) != -999:
                TessRS.append(data[i][12])
                Tessperiods.append(data[i][12])
        elif data[i][16] == 'SSG':
            SSGsP.append(data[i][11])
            test.append(data[i][16])
            periods.append(data[i][11])
            if float(data[i][12]) != -999:
                TessSSGs.append(data[i][12])
                Tessperiods.append(data[i][12])
        elif data[i][16] != 'cut':
            extra.append(data[i][16])
            extraind.append(i)
            extraid.append(data[i][len(data[i])-1])
            #periods.append(data[i][11])
    '''print()
    print("periods length: ", len(periods))
    print("RS length: ", len(RSsP))
    print("SSG length: ", len(SSGsP))
    print("Extra length: ", len(extra))
    print()
    print("TESS RS length: ", len(TessRS))
    print("TESS SSGs length: ", len(TessSSGs))
    print("TESS Total length: ", len(Tessperiods))
    print()
    print(extra)
    print(extraind)
    print(extraid)'''
    periods = [float(item) for item in periods]
    SSGsP = [float(item) for item in SSGsP]
    RSsP = [float(item) for item in RSsP]
    #print(len(periods))
    #print(periods)
    #print(periods)
    
    Tessperiods = [float(item) for item in Tessperiods]
    TessSSGs = [float(item) for item in TessSSGs]
    TessRS = [float(item) for item in TessRS]
    
    
    
    '''
    fig = plt.figure(figsize = (10,4))
    ax1 = fig.add_subplot(121)
    
    
    
    ax1.hist(periods, 50, color = 'r', label = "Total Periods")
    ax1.hist(RSsP, 50, color = 'orange', label = "RS")#, density = True)
    ax1.hist(SSGsP, 50, color = 'g', label = "SSGs")#, density = True)
    ax1.set_title("Leiner 2022 Samples (VSX)")
    ax1.set_xlabel("Periods (Days)")
    #ax1.set_xlim(right = 50)
    plt.legend()
    
    
    ax2 = fig.add_subplot(122)
    ax2.hist(Tessperiods, 50, color = 'r', label = "Total Periods")
    ax2.hist(TessRS, 50, color = 'orange', label = "RS")#, density = True)
    ax2.hist(TessSSGs, 50, color = 'g', label = "SSGs")#, density = True)
    ax2.set_title("Leiner 2022 Samples (TESS)")
    ax2.set_xlabel("Periods (Days)")
    #ax2.set_xlim(right = 50)
    plt.legend()
    plt.show()'''
    
    
    
    
    SSGperiods = SSGsP
    Tessperiods = TessSSGs
    
    
    SSGperiods = RSsP
    
    randPvals = []
    for i in range(1000):
        ###
        hist_counts, bins = np.histogram(SSGperiods, bins=25, density=True)
            
        # Compute the PDF
        pdfdata = hist_counts / sum(hist_counts)
            
        # Compute the CDF
        cdfdata = np.cumsum(pdfdata)
        cdfdata = np.insert(cdfdata, 0, 0)  # Insert 0 at the beginning for the first bin edge
            
        # Generate a uniform random number
        uniform_random_value = np.random.rand()
            
        # Map the uniform random number to the histogram distribution
        PorbRand = np.interp(uniform_random_value, cdfdata, bins)
        randPvals.append(PorbRand)
        #print(f"Random value drawn from the histogram distribution: {PorbRand}")
        
        ###
    '''
    plt.hist(randPvals, bins = 25)
    plt.xlabel("Periods (Days)")
    plt.show()'''

    
    #####################################################
    #####################################################
    #Begin CDF Generation


    
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
    
    
    from tqdm import tqdm
    
    # Adding progress bar to the loop with tqdm
    # Generate CDFs for each run.
    for i in tqdm(range(runs), desc=f"Run {run_number}"):
        
        Msol = 1.9891 * 10**30 #kg
        G = 6.67430e-11  # m^3 kg^-1 s^-2
        ###############
        #Create Period options
        
        
        
        #Lognormal Distribution
        ux = 4.8
        ox = 2.3
        mode = 4.8
        #Raghavan numbers: ux = 5.03, ox = 2.28, mode = 5.03
        u = np.log( (ux**2) / np.sqrt(ux**2 + ox**2))
        o = np.sqrt( np.log(1 + (ox**2 / ux**2)))
        um = np.log(mode) + o**2
            
        #print("Ux = ", ux)
        #print("Ox = ", ox)
        #print("Sigma input: ", o)
        #print("Mu input:", um)
        #print("------------------------------")
        #print()
        
        m1 = random.uniform(0.9, 2)
        m2 = random.uniform(0.1, m1)
        def lognormal(mean, sigma, range):
            while True:
                logP = np.random.lognormal(mean, sigma)
                if logP>0 and logP < range:
                    break
            return logP
        
        
        while True:
            logP = lognormal(um, o, 1.5) #10^1.5 = 31 days
            logvalue = 10**logP
        
            if logvalue <= 30 and logvalue >= 2:
                break  # Accept the value if it's within the desired range
    
        while True:
            logP1000 = lognormal(um, o, 3)  #10^3 = 1000 days
            logvalue1000 = 10**logP1000
    
            if logvalue1000 <= 1000 and logvalue1000 >=2:
                break  # Accept the value if it's within the desired range

        
        logPfull = lognormal(um, o, 30)  # 10^30  days, don't end up using this, this is just the full distribution.
        
        #Tidally Locked Distribution
        hist_counts, bins = np.histogram(SSGperiods, bins=25, density=True)
        # Compute the PDF
        pdfdata = hist_counts / sum(hist_counts)  
        # Compute the CDF
        cdfdata = np.cumsum(pdfdata)
        cdfdata = np.insert(cdfdata, 0, 0)  # Insert 0 at the beginning for the first bin edge
    
        
        while True:
            uniform_random_value = np.random.rand()
            PorbRand = np.interp(uniform_random_value, cdfdata, bins)
        
            if PorbRand <= 35:
                break  # Accept the value if it's within the desired range
    
        
        logtidal = np.log10(PorbRand)
        
        logPvals.append(logP)
        logPfullvals.append(logPfull)
        logP1000vals.append(logP1000)
        logPtidal.append(logtidal)
        ###########################
        ###########################
    
        
    
        
        ###########################
        ###########################     
        
        #Simulate Mass Ratios, Inclinations
        M1 = m1 * Msol
        M2 = m2 * Msol
        v = random.uniform(0,1)
        incl = np.arccos(2 * v - 1)
        inclinations.append(incl)
        
        ###########################
        ###########################
        #Decide on Period
        
        b1 = 10**logP                   #Lognormal, 0-30 days
        #b = 10**logPfull               #Lognormal, full range 0-inf days
        b2 = 10**logP1000               #Lognormal, 0-1000 days
        #b = np.random.uniform(1,30)    #Uniform, 1-30 days
        #b = np.random.uniform(1,1000)  #Uniform, 1-1000 days
        b3 = 10 ** logtidal                   #Tidally Locked sample
        
        ###########################
        ###########################
        #Generate Orbit
        x = np.linspace(0, 365*2, 10000) 
        def genorbit(b):
            Porb = b * 86400
            #Period range, 0-2 years 
                
            #Generate amplitude from orbit, mass ratios
            K = (  (M2**3 * np.sin(incl)**3 * (2 * np.pi * G)) / (Porb * (M1 + M2)**2)  )**(1/3)
            #Ke = (  (M2**3 * np.sin(incl)**3 * (2 * np.pi * G)) / (Porb * (M1 + M2)**2  * (1-e**2)**3/2)  )**(1/3)
            #Kvals.append(K/1000)
            #Kvalse.append(Ke/1000)
            
            a = K/1000                      #Amplitude
            c = np.random.random() * 10     #Shift
            #d = np.random.random() * np.random.random()*10           #Shift up/down
            
            #print("M1:", M1, "Msol")
            #print("M2:", M2, "Msol")
            #print("Inclination:", np.rad2deg(incl), "Degrees")
                
            #print("Amplitude (km/s):", a)
            #print("Period (Days):", b)
            #print("Shift:", c)
                
            #print(a, "* sin(2pi/", b, "*x +", c, ")")# + d")
            #print()    
            equation = a* np.sin(((2*np.pi)/b)*x + c)
            return equation, a, b
    
        ###########################
        ###########################
        #Generate Eccentric Orbit
        
        def genorbit_e(Period):   
            
            P = Period
            Porb = Period * 86400
            
            #print(f"P: {Period} days")
                
            if np.log10(P) < 1:
                e = 0
            elif np.log10(P) >= 1 and np.log10(P) < 2:
                e = random.uniform(0, 0.4)
            elif np.log10(P) >= 2:
                e = random.uniform(0.2, 0.6)
            eccens.append(e)
            eccens_logPvals.append(np.log10(Period))
        
                
            Ke = (  (M2**3 * np.sin(incl)**3 * (2 * np.pi * G)) / (Porb * (M1 + M2)**2  * (1-e**2)**3/2)  )**(1/3)
            #print("Ke:", Ke/1000)
            
            # Orbital Parameters
            omega = np.deg2rad(np.random.uniform(0, 360))  # Argument of periapsis in radians
        
            t = x
            n = 2 * np.pi / P  # Mean motion (rad/day)
            
            # Mean anomaly M(t)
            M = n * t
            
            # Solve Kepler's equation using Newton-Raphson method
            #This improves speed of calculation, as opposed to scipy fsolve
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
            
            # Compute True Anomaly nu(t)
            nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
            
            # Compute Radial Velocity v_r(t)
            v_r = Ke * (np.cos(omega + nu) + e * np.cos(omega))
            
            # Plot results
            return v_r/1000, Ke/1000, P
    
    
        ###########################
        ###########################
        #Create different orbits
        #Will only use three (Eccentric Long period, Eccentric Shrot period, Circular Tidally synchronized)
        
        orbit1, a1, b1 = genorbit(b1)
        orbit2, a2, b2 = genorbit(b2)
        orbit3, a3, b3 = genorbit(b3)
    
    
        orbit1e, a1e, b1e = genorbit_e(b1)
        orbit2e, a2e, b2e = genorbit_e(b2)
        orbit3e, a3e, b3e = genorbit_e(b3)
    
        #orbit = orbit1
        ###########################
        ###########################
        
        #Generate random samples

        
        samples = random.randint(lower_samp, upper_samp)
        #print("Number of Samples:", samples)
        
        ###########################
        ###########################
        # Generate Oribt Sampling
        
        #Sample from random orbit
        def sample(orbit, samples):
            xval = []
            yval = []
            obsday = 0
            prevsamp = int((np.random.random()*10000) - 1)
            for k in range(samples):
                obsday = obsday + 1
                #print(obsday)
                if obsday < int(np.random.random()*2)+2:
                    #print("<rand")
                    #print(prevsamp)
                    samp = prevsamp + int(np.random.random()*75)
                    if samp > 9999:
                        samp = 9999
                        #print("Limit Hit")
                    prevsamp = samp
                    #print("Km/s:", orbit[samp])
                    #print("Day:", x[samp])
                    yval.append(orbit[samp])
                    xval.append(x[samp])
                    #print(samp)
                    
                else:
                    #print("else")
                    samp = int((np.random.random()*10000))
                    #print(prevsamp)
                    if samp > 9999:
                        samp = 9999
                        #print("Limit Hit")
                    prevsamp = int((np.random.random()*10000) - 1)
                    prevsamp = samp
                    #print("Km/s:", orbit[samp])
                    #print("Day:", x[samp])
                    
                    yval.append(orbit[samp])
                    xval.append(x[samp])
                    obsday = 0
                    #print(samp)
                    #print(obsday)
                #print()
            return xval, yval
    
        xval1, yval1 = sample(orbit1, samples)
        xval2, yval2 = sample(orbit2, samples)
        xval3, yval3 = sample(orbit3, samples)
    
    
        xval1e, yval1e = sample(orbit1e, samples)
        xval2e, yval2e = sample(orbit2e, samples)
        xval3e, yval3e = sample(orbit3e, samples)
    
    
        
        def plotorbs(xval, yval, orbit):
                
            
                    #print()
            fig, axs = plt.subplots(1, 2, figsize=(14, 4)) 
            # First plot
            axs[0].scatter(x, orbit, marker=".", s=1)
            axs[0].axhline(0, color="red")
            axs[0].set_title("Simulated Orbit")
            axs[0].set_xlabel("Period (days)")
            axs[0].set_ylabel("Velocity (km/s)")
            axs[0].set_xlim(-10,(365*2)+10)
                
            # Second plot
            axs[1].scatter(xval, yval, marker=".")
            axs[1].plot(x, orbit, alpha = 0.1)
            axs[1].axhline(0, color="red", alpha=0.5)
            axs[1].set_title("Simulated Observations")
            axs[1].set_xlabel("Period (days)")
            axs[1].set_ylabel("Velocity (km/s)")
            axs[1].set_xlim(-10,(365*2)+10)
                
            plt.tight_layout()  # Adjust spacing between plots
            plt.show()

        # Plot orbits if desired
        #plotorbs(xval1, yval1, orbit1)
        #plotorbs(xval2, yval2, orbit2)
        #plotorbs(xval3, yval3, orbit3)
        
        #plotorbs(xval1e, yval1e, orbit1e)
        #plotorbs(xval2e, yval2e, orbit2e)
        #plotorbs(xval3e, yval3e, orbit3e)
        
        ###########################
        ###########################
        #Calculate Amplitude of this sample (measured K)
        
        def measure(yval, a):
            std = np.std(yval)
            #print("Standard Deviation:", std)
            Amp = (np.max(yval) - np.min(yval)) / 2
            
            
            #print("Measured Amplitude (K Value):", Amp) 
            #print("Actual Amplitude (K Value):", a)
            #print()
            return Amp, std
    
        Amp1, std1 = measure(yval1, a1)
        Amp2, std2 = measure(yval2, a2)
        Amp3, std3 = measure(yval3, a3)
    
        Amp1e, std1e = measure(yval1e, a1e)
        Amp2e, std2e = measure(yval2e, a2e)
        Amp3e, std3e = measure(yval3e, a3e)
        
        ###########################
        ###########################
        # Record P, actual K, measured K for analysis
        
        #Do correlation:
        #corr = Amp/a
        #totcorr.append(corr * 100)
        #print("Corr:", corr*100)
    
        actual1.append(a1)
        actual2.append(a2)
        actual3.append(a3)
        actual1e.append(a1e)
        actual2e.append(a2e)
        actual3e.append(a3e)
    
        measured1.append(Amp1)
        measured2.append(Amp2)
        measured3.append(Amp3)
        measured1e.append(Amp1e)
        measured2e.append(Amp2e)
        measured3e.append(Amp3e)
    
        periods1.append(b1)
        periods2.append(b2)
        periods3.append(b3)
        periods1e.append(b1e)
        periods2e.append(b2e)
        periods3e.append(b3e)
    
        Kvals1.append(a1)
        Kvals2.append(a2)
        Kvals3.append(a3)
        Kvals1e.append(a1e)
        Kvals2e.append(a2e)
        Kvals3e.append(a3e)
    
        #print()
        #print(Amp1, a1, b1)
        #print(Amp2, a2, b2)
        #print(Amp3, a3, b3)
    
        #################

    #Now, generate CDF plots
    #Inlcuded a few different options for CDF setups
    #Choose [actual2e, actual1e, actual3], [measured2e, measured1e, measured3] for our desired plots.
    #(Eccentric Long period, Eccentric Shrot period, Circular Tidally synchronized)

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import rankdata
    
    plt.figure(figsize=(10, 6), dpi = 300)
    '''
    for actual, measured, label, color in zip([actual1, actual2, actual3], [measured1, measured2, measured3], 
                                              ['Lognormal (0-30 days)', 'Lognormal (0-1000 days)', 'Tidally Locked (Leiner SSGs Distribution)'],
                                              ['red', 'green', 'blue']):
        actual_sorted = np.sort(actual)
        actual_cdf = rankdata(actual_sorted, method='average') / len(actual_sorted)
        measured_sorted = np.sort(measured)
        measured_cdf = rankdata(measured_sorted, method='average') / len(measured_sorted)
        
        plt.plot(actual_sorted, actual_cdf * 100, linestyle='-', color=color, alpha = 0.5)
        plt.plot(measured_sorted, measured_cdf * 100, linestyle='--', color=color, alpha = 0.5)
        plt.fill_betweenx(actual_cdf*100, actual_sorted, measured_sorted, color=color, alpha=0.25, label = f"{label}")'''
    
    
    ##################
    ##################
    '''
    for actual, measured, label, color in zip([actual1e, actual2e, actual3e], [measured1e, measured2e, measured3e], 
                                              ['Eccentric, Lognormal (0-30 days)', 'Eccentric, Lognormal (0-1000 days)', 'Eccentric, Tidally Locked (Leiner SSGs Distribution)'],
                                              ['orangered', 'lightgreen', 'cyan']):
        actual_sorted = np.sort(actual)
        actual_cdf = rankdata(actual_sorted, method='average') / len(actual_sorted)
        measured_sorted = np.sort(measured)
        measured_cdf = rankdata(measured_sorted, method='average') / len(measured_sorted)
        
        plt.plot(actual_sorted, actual_cdf * 100, linestyle='-', color=color, alpha = 0.5)
        plt.plot(measured_sorted, measured_cdf * 100, linestyle='--', color=color, alpha = 0.5)
        plt.fill_betweenx(actual_cdf*100, actual_sorted, measured_sorted, color=color, alpha=0.25, label = f"{label}")'''
    
    ##################
    ##################
    '''
    for actual, measured, label, color in zip([actual1e, actual2e], [measured1e, measured2e], 
                                              ['Eccentric, Lognormal (0-30 days)', 'Eccentric, Lognormal (0-1000 days)'],
                                              ['magenta', 'lightgreen']):
        actual_sorted = np.sort(actual)
        actual_cdf = rankdata(actual_sorted, method='average') / len(actual_sorted)
        measured_sorted = np.sort(measured)
        measured_cdf = rankdata(measured_sorted, method='average') / len(measured_sorted)
        
        plt.plot(actual_sorted, actual_cdf * 100, linestyle='-', color=color)
        plt.plot(measured_sorted, measured_cdf * 100, linestyle='--', color=color)
        plt.fill_betweenx(actual_cdf*100, actual_sorted, measured_sorted, color=color, alpha=0.3, label = f"{label}")
    '''
    
    ##################
    ##################
    
    for actual, measured, label, color in zip([actual2e, actual1e, actual3], [measured2e, measured1e, measured3], 
                                              ['Eccentric, Lognormal (P < 1000 days)', 'Eccentric, Lognormal (P < 30 days)', 'Tidally Locked (Leiner SSGs Distribution)'],
                                              ['green', 'red', 'blue']):
        actual_sorted = np.sort(actual)
        actual_cdf = rankdata(actual_sorted, method='average') / len(actual_sorted)
        measured_sorted = np.sort(measured)
        measured_cdf = rankdata(measured_sorted, method='average') / len(measured_sorted)
        if dispCDF == True:
            plt.plot(actual_sorted, actual_cdf * 100, linestyle='-', color=color)
            plt.plot(measured_sorted, measured_cdf * 100, linestyle='--', color=color)
            plt.fill_betweenx(actual_cdf*100, actual_sorted, measured_sorted, color=color, alpha=0.3, label = f"{label}")
    
    
    
    ##################
    ##################
    #Original data values, copied from txt files
    
    sampnums = [5, 4, 3, 8, 4, 9, 11, 7, 3, 7, 12, 12, 12, 5, 13, 9, 9, 11, 5, 6, 13,
                 4, 11, 12, 9, 10, 5, 5, 4, 9, 6, 12, 4, 3, 5, 5, 15, 3, 11, 6, 3, 10, 5,
                 5, 20, 10, 17, 3, 11, 8, 11, 4, 6, 11, 5, 6, 4]
    amps = [0.9850000000000003, 12.76, 3.4899999999999984, 0.5749999999999993, 47.065, 55.515,
             28.86, 34.535, 1.1100000000000012, 12.155000000000001, 26.020000000000003, 39.035000000000004, 
             49.955, 31.655, 43.735, 3.335000000000001, 14.43, 23.58, 74.89, 21.265, 0.86, 20.34, 21.5, 
             46.72, 27.32, 0.504999999999999, 0.27500000000000036, 18.595, 0.27000000000000046, 39.67, 61.86, 
             32.93, 27.565, 18.915, 12.08, 7.54, 51.705, 7.985000000000001, 58.185, 25.81, 62.315, 
             26.490000000000002, 35.105, 15.275, 26.19, 12.915, 38.31, 9.4, 52.47, 46.435, 51.935, 28.075, 47.86, 
             30.625, 10.285, 1.0350000000000001, 16.474999999999998]
    
    ampnew = []
    sampnew = []
    
    #print(len(sampnums), len(amps))

    #Create a 'good' dataset, where only values with a certain number of samples are included
    #The number is set by the user when running the function
    #We generally assume 7+ samples is "good"
    #We also include values with amplitdues above 50 km/s, as these are clearly extreme variable/fast rotating stars
    
    samplelim = limit_on_samples
    for i in range(len(sampnums)):
        if sampnums[i] >= samplelim:
            ampnew.append(amps[i])
            sampnew.append(sampnums[i])
        if sampnums[i] < samplelim and amps[i] > 50:
            ampnew.append(amps[i])
            sampnew.append(sampnums[i])
    
    
    #print(len(ampnew))
    #print(len(amps))
    STDs = amps
    
    
    
    #Begin creating CDFs;
    #If dispCDF = True, will plot the cdf.
    
    vals1 = np.sort(amps)
    p1 = np.arange(len(vals1)) / (len(vals1) - 1)
    coefficients = np.polyfit(vals1, p1, 3)
    best_fit_line1 = coefficients[0]*vals1**3 + coefficients[1]*vals1**2  + coefficients[2]*vals1 + coefficients[3] 
    
    #plt.scatter(vals1, p1*100, color = 'g', marker = '.', label = "Coude Sample Deviations", zorder = 2, s = 25)
    #plt.plot(vals1, best_fit_line1*100, color = 'g', linewidth = 6, alpha = 0.35, zorder = 3)#, Intersect: " + str(np.round(best_fit_line[slot3-1] * 100, 2)) + "%", zorder = 4)
    
    
    
    
    vals2 = np.sort(ampnew)
    p2 = np.arange(len(vals2)) / (len(vals2) - 1)
    coefficients = np.polyfit(vals2, p2, 3)
    
    best_fit_line2 = coefficients[0]*vals2**3 + coefficients[1]*vals2**2  + coefficients[2]*vals2 + coefficients[3] 
    #plt.scatter(vals2, p2*100, color = 'g', marker = '.', label = "Coude Sample Deviations", zorder = 2, s = 25)
    #plt.plot(vals2, best_fit_line2*100, color = 'g', linewidth = 6, alpha = 0.35, zorder = 3)#, Intersect: " + str(np.round(best_fit_line[slot3-1] * 100, 2)) + "%", zorder = 4)
    
    # Ensure we use the same x-values for both fits
    common_x = np.linspace(min(vals1.min(), vals2.min()), max(vals1.max(), vals2.max()), 100)
    
    # Evaluate both polynomial fits at these x-values
    best_fit_line1_common = np.polyval(np.polyfit(vals1, p1, 3), common_x)
    best_fit_line2_common = np.polyval(np.polyfit(vals2, p2, 3), common_x)



    
    if dispCDF == True:
        # Fill the area between the two curves
        #plt.fill_between(common_x, best_fit_line1_common * 100, best_fit_line2_common * 100, color='black', alpha=0.25, label='Our Observations', zorder=5)
        # Plot original data and fits
        plt.scatter(vals1, p1 * 100, color='black', marker='x', label=f"Observations with 3+ Measurements", zorder=6, s=50)
        #plt.plot(common_x, best_fit_line1_common * 100, color='g', linewidth=6, alpha=0.35, zorder=3)
        plt.scatter(vals2, p2 * 100, color='black', marker='.', label=f"Observations with {samplelim}+ Measurements", zorder=6, s=50)
        #plt.plot(common_x, best_fit_line2_common * 100, color='g', linewidth=6, alpha=0.35, zorder=5)
        
        
        
        
        
        ##################
        ##################
        
        
        
        
        
        
        plt.xlabel('K (km/s)')
        plt.ylabel('Probability (%)')
        #plt.title(f'Run {run_number}: CDF of Simulated and Actual Observation Methods')
        plt.legend(loc='lower right', fontsize=9)
        plt.xlim(-5, 105)
        plt.ylim(-5, 105)
        #plt.xlim(105, -5)
        #plt.ylim(105, -5)
        plt.legend()
        plt.grid(alpha = 0.5, linestyle = "--")
        plt.tight_layout()
        #plt.savefig("Sim_Full_Eccens_No_Fill.png", dpi = 500)
        #plt.savefig("CDF1.png", dpi = 300)
        
        if SaveThisFig == True:
            plt.savefig(savefig_name + ".png", dpi = 300)
        plt.show()

        
    ############################
    ############################
    #Begin Kolmogorov–Smirnov testing


    
    from scipy.stats import ks_2samp
    '''print(f"""
    KS Test - Comparing two CDF values to identify similarities in CDF and underlying distribution:
        KS statistic: Max difference between both CDFs, value of 0-1 (0-100%, y axis)
            Large KS = Large difference = LESS similarity
            Small KS = Small difference = MORE similarity
        
        p-value: are these data sets drawn from the same underlying distribution?
            Large p (p > 0.05) = sets are drawn from SIMILAR (or same) distribution
            Small p (p < 0.05) = sets are drawn from DIFFERENT distributions
    
    The following code looks at the difference between our sampled data
        Full Distribution
        {samplelim}+ measurements
    It compares these datasets to the "Measured" and "Actual" datasets created from the three underlying distributions we want to compare:
        Eccentric, Lognormal (0-1000 days)
        Eccentric, Lognormal (0-30 days)
        Tidally Locked (Leiner SSGs Distribution)
    
    So, we're looking for the results with the LOWEST KS value, and HIGHEST p values:
    
    """)'''
    #print("##"*25)
    #print()
    aKSval = []
    mKSval = []
    apval = []
    mpval = []



    #KS test for comparison.
    #Want to compare each actual and measured set to the amplitudes we measured in our real datset
    #Uses ks_2samp to draw ks statistics (KS-value) and p values from comparing distributions
    #First, 'Measured' vs Amps
    #Second, 'Actual' vs Amps
    
    def KS(amps, start):
        #print(start)
        for actual, measured, label, color in zip([actual2e, actual1e, actual3], [measured2e, measured1e, measured3], 
                                                  ['Eccentric, Lognormal (0-1000 days)', 'Eccentric, Lognormal (0-30 days)', 'Tidally Locked (Leiner SSGs Distribution)'],
                                                  ['green', 'red', 'blue']):
            #print(label)
            #"Measured" values compared to the real amplitudes
            ks_statistic, p_value = ks_2samp(amps, measured)
            #print(f"  Measured: KS Statistic = {ks_statistic:.4f}, P-value = {p_value:.4f}")
            mKSval.append(ks_statistic)
            mpval.append(p_value)
            # Interpretation
            alpha = 0.05  # Significance level
            '''if p_value < alpha:
                print("p < 0.05; The distributions are different.")
            else:
                print("p > 0.05; The distributions are similar.")'''
    
    
            #"Actual" value compared to the real amplitudes
            ks_statistic, p_value = ks_2samp(amps, actual)
            #print(f"  Actual:   KS Statistic = {ks_statistic:.4f}, P-value = {p_value:.4f}")
            aKSval.append(ks_statistic)
            apval.append(p_value)
            '''
            # Interpretation
            alpha = 0.05  # Significance level
            if p_value < alpha:
                print("p < 0.05; The distributions are different.")
            else:
                print("p > 0.05; The distributions are similar.")'''
            #print()
        #print("##"*25)
        #print()


    #run KS Test
    #Two sets of data: Full dataset (amps array) and limited dataset (ampsnew array)
    # ampsnew is limited by the number of samples set in the function
    # e.g, 7+ samples in the dataset is considered "accurate" or "good" amplitude values from the real dataset.
    KS(amps, "Full Distribution")
    KS(ampnew, f"Observations with {samplelim}+ Measurments")
    
    
    if PrintTestResults == True:
        print("KS, measured sets:", [f"{val:.4f}" for val in mKSval])
        print("KS, actual sets:  ", [f"{val:.4f}" for val in aKSval])
        print()
        print("p-val, measured sets:", [f"{val:.4f}" for val in mpval])
        print("p-val, actual sets:  ", [f"{val:.4f}" for val in apval])

    plt.close()


    #Return KS and p values
    return mKSval, aKSval, mpval, apval
