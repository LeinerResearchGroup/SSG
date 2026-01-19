

def runsamples(runs, orbsamplerange, limit_on_samples, run_number, 
               dispCDF = False, dispElse = False, elseerror = False, saveElse = False,
               SaveThisFig = False, PrintTestResults = False, savefig_name = "Test", 
               show = True, ax = None):

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
    ec_all = []

    actualinfe, measuredinfe, periodsinfe, Kvalsinfe = [], [], [], []
    
    from tqdm import tqdm
    
    # Adding progress bar to the loop with tqdm
    # Generate CDFs for each run.
    for i in tqdm(range(runs), desc=f"Run {run_number}"):
        
        Msol = 1.9891 * 10**30 #kg
        G = 6.67430e-11  # m^3 kg^-1 s^-2
        ###############
        #Create Period options
        
        
        
        #Lognormal Distribution
        ux = 4.8   #5.03 / 4.8
        ox = 2.3   #2.28 / 2.3
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
                if logP>np.log10(2) and logP < range:
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

        while True:
            logPfull = lognormal(um, o, 10)  # 10^10  
            logvaluefull = 10**logPfull
            if logvaluefull >= 2:
                break   #Accept if the logpfullval is greater than 2 days. upper limit is already 10^10

        #IF we want to see very short period performance, here's a quick rough replacement of the long-period samples:
        #while True:
            #logPfull = np.random.uniform(-0.0, 0.69897) #0.2 days to 5 days, rough uniform distribution to see estimate
            #logPfull = np.random.uniform(0.69897, 1) #5 to 10 days
            #logPfull = np.random.uniform(1, 1.30103) #10 to 20 days
            #logPfull = np.random.uniform(1.30103, 1.5) #20 to 31 days
        #    logvaluefull = 10**logPfull
        #    if logvaluefull <= 30: #limits to 30 days max
        #        break   #Accept if the logpfullval is less than 5 days


        
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
        binf = 10**logPfull               #Lognormal, full range 0-inf days
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
            return equation, a, b, 0
    
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
            elif np.log10(P) >= 2 and np.log10(P) < 3:
                e = random.uniform(0.2, 0.6)
                if elseerror:
                    e = e + random.uniform(-0.19, 0.19)
            elif np.log10(P) >=3:
                e = random.uniform(0.2, 0.8)
                if elseerror:
                    e = e + random.uniform(-0.19, 0.19)
            
            ecc = e
            #eccens.append(e)
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
            return v_r/1000, Ke/1000, P, ecc
    
    
        ###########################
        ###########################
        #Create different orbits
        #Will only use three (Eccentric Long period, Eccentric Shrot period, Circular Tidally synchronized)
        
        #orbit1, a1, b1 = genorbit(b1)
        #orbit2, a2, b2 = genorbit(b2)
        orbit3, a3, b3, e3 = genorbit(b3) #tidal
    
    
        orbit1e, a1e, b1e, e1e = genorbit_e(b1)  #30 days
        orbit2e, a2e, b2e, e2e = genorbit_e(b2)  #1000 days
        #orbit3e, a3e, b3e = genorbit_e(b3)

        orbitinfe, ainfe, binfe, einfe = genorbit_e(binf)  #10^10
        #orbit = orbit1



        ec_all.append([einfe, e2e, e1e, e3])
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
    
        #xval1, yval1 = sample(orbit1, samples)
        #xval2, yval2 = sample(orbit2, samples)
        xval3, yval3 = sample(orbit3, samples)
    
    
        xval1e, yval1e = sample(orbit1e, samples)
        xval2e, yval2e = sample(orbit2e, samples)
        #xval3e, yval3e = sample(orbit3e, samples)
    
        xvalinfe, yvalinfe = sample(orbitinfe, samples)
        
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
    
        #Amp1, std1 = measure(yval1, a1)
        #Amp2, std2 = measure(yval2, a2)
        Amp3, std3 = measure(yval3, a3)
    
        Amp1e, std1e = measure(yval1e, a1e)
        Amp2e, std2e = measure(yval2e, a2e)
        #Amp3e, std3e = measure(yval3e, a3e)

        Ampinfe, stdinfe = measure(yvalinfe, ainfe)
        
        ###########################
        ###########################
        # Record P, actual K, measured K for analysis
        
        #Do correlation:
        #corr = Amp/a
        #totcorr.append(corr * 100)
        #print("Corr:", corr*100)
    
        #actual1.append(a1)
        #actual2.append(a2)
        actual3.append(a3)
        actual1e.append(a1e)
        actual2e.append(a2e)
        #actual3e.append(a3e)
        actualinfe.append(ainfe)
    
        #measured1.append(Amp1)
        #measured2.append(Amp2)
        measured3.append(Amp3)
        measured1e.append(Amp1e)
        measured2e.append(Amp2e)
        #measured3e.append(Amp3e)
        measuredinfe.append(Ampinfe)
    
        #periods1.append(b1)
        #periods2.append(b2)
        periods3.append(b3)
        periods1e.append(b1e)
        periods2e.append(b2e)
        #periods3e.append(b3e)
        periodsinfe.append(binf)
    
        #Kvals1.append(a1)
        #Kvals2.append(a2)
        Kvals3.append(a3)
        Kvals1e.append(a1e)
        Kvals2e.append(a2e)
        #Kvals3e.append(a3e)
        Kvalsinfe.append(ainfe)
    
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

    logPvals = np.array(logPvals)
    logPfullvals = np.array(logPfullvals)
    logP1000vals = np.array(logP1000vals)
    logPtidal = np.array(logPtidal)

    if dispElse:
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))  # 1 row, 4 columns
        
        axes[0].hist(logPfullvals, 15, histtype = 'stepfilled', 
                     facecolor = 'b', edgecolor = 'b', alpha = 0.6, linewidth = 2.5)
        axes[0].set_title(r"Porb Values, (P: 2-$10^{10}$ days)")
        axes[0].set_xlabel("Log10(P)")
        
        axes[1].hist(logP1000vals, 15, histtype = 'stepfilled', 
                     facecolor='g', edgecolor = 'g', alpha = 0.6, linewidth = 2.5)
        axes[1].set_title("Porb Values, (P: 2-1000 days)")
        axes[1].set_xlabel("Log10(P)")
        
        axes[2].hist(logPvals, 15,histtype = 'stepfilled',
                     facecolor='r', edgecolor = 'r', alpha = 0.5, linewidth = 2.5)
        axes[2].set_title("Porb Values, (P: 2-30 days)")
        axes[2].set_xlabel("Log10(P)")
        
        axes[3].hist(logPtidal, 15, histtype = 'stepfilled', 
                     facecolor='y', edgecolor = 'y', alpha = 0.6, linewidth = 2.5)
        axes[3].set_title("Porb Values, (Leiner SSGs Distribution)")
        axes[3].set_xlabel("Log10(P)")
        
        plt.tight_layout()
        if saveElse:
            plt.savefig("Distributions.png", dpi = 300)
        plt.show()
        plt.close()


    
    #print(ec_all)
    ec_all = np.array(ec_all)
    if dispElse:
            
        fig, axes = plt.subplots(figsize=(8,4))
        axes.scatter(logPfullvals, ec_all[:,0], label=r"$10^{10}$ d",   color='b', alpha = 0.5, s = 5)
        axes.scatter(logP1000vals, ec_all[:,1], label='1000 d',color='g', alpha = 0.5, s= 5)
        axes.scatter(logPvals, ec_all[:,2], label='30 d',  color='r', alpha = 0.5, s = 5)
        axes.scatter(logPtidal, ec_all[:,3], label='Leiner 2022', color='y', alpha = 0.5, s = 5)
        axes.legend()
        axes.set_xlabel("Log P")
        axes.set_ylabel("Eccentricities")
        plt.tight_layout()
        if saveElse:
            plt.savefig("Eccentricities.png", dpi = 300)
        plt.show()
        plt.close()



    
    
    #plt.figure(figsize=(10, 6), dpi = 300)
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,6), dpi = 300)
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
    
    for actual, measured, label, color in zip([actualinfe, actual2e, actual1e, actual3], [measuredinfe, measured2e, measured1e, measured3], 
                                              [r'Eccentric, Lognormal (P: 2-$10^{10}$ days)', 'Eccentric, Lognormal (P: 2-1000 days)', 
                                               'Eccentric, Lognormal (P: 2-30 days)', 'Tidally Locked (Leiner SSGs Distribution)'],
                                              ['blue', 'green', 'red', 'y']):
        actual_sorted = np.sort(actual)
        actual_cdf = rankdata(actual_sorted, method='average') / len(actual_sorted)
        measured_sorted = np.sort(measured)
        measured_cdf = rankdata(measured_sorted, method='average') / len(measured_sorted)
        if dispCDF == True:
            ax.step(actual_sorted, actual_cdf * 100, linestyle='-', color=color, where = 'mid')
            ax.step(measured_sorted, measured_cdf * 100, linestyle='--', color=color, where = 'mid')
            ax.fill_betweenx(actual_cdf*100, actual_sorted, measured_sorted, color=color, alpha=0.3, label = f"{label}")
    
    
    
    ##################
    ##################
    nan = np.nan
    #Star, N, RV1_amp, rv1_hi, rv1_lo, RV2_amp, rv2_hi, rv2_lo
    data = [
    ['29Dra', 5, 0.992, 0.151, 0.167, nan, nan, nan, 'SSG', 'SB1'],
    ['2MASSJ08171221+0736164', 4, 12.751, 0.193, 0.153, nan, nan, nan, 'RS', 'SB1'],
    ['2MASSJ08504952+1217158', 4, 18.098, 0.324, 0.380, 4.965, 0.369, 0.565, '--', 'SB2'],
    ['2MASSJ17270871+2700144', 3, 3.493, 0.051, 0.045, nan, nan, nan, 'SSG', 'SB1'],
    ['54Cam', 5, 49.191, 0.157, 0.155, 48.586, 0.248, 0.225, 'RS', 'SB2'],
    ['ASASJ073642+0354.3', 8, 0.576, 0.108, 0.091, nan, nan, nan, 'SSG', 'SB1'],
    ['ASASJ155307+2028.6', 4, 47.053, 0.211, 0.258, nan, nan, nan, 'SSG', 'SB1'],
    ['ASASJ162510+0514.9', 9, 55.507, 0.970, 1.043, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['ASASJ210323+1930.9', 12, 28.870, 0.344, 0.367, nan, nan, nan, 'SSG', 'SB1'],
    ['BD+142519', 7, 34.534, 0.264, 0.272, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['BD+201895', 3, 1.123, 0.197, 0.210, nan, nan, nan, 'RS', 'Asymmetric'],
    ['BMCVn', 7, 12.139, 0.136, 0.153, nan, nan, nan, 'RS', 'SB1'],
    ['BPSCS22950-0027', 12, 26.011, 0.153, 0.150, nan, nan, nan, 'SSG', 'SB1'],
    ['CGTri', 12, 39.050, 0.654, 0.553, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['EPCVn', 12, 49.956, 0.106, 0.117, nan, nan, nan, 'Cut', 'SB1'],
    ['FMLyn', 6, 31.651, 0.287, 0.290, nan, nan, nan, 'RS', 'SB1'],
    ['FNCom', 13, 43.760, 0.303, 0.297, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['FRBoo', 9, 3.316, 0.161, 0.156, nan, nan, nan, 'RS', 'SB1'],
    ['FVCnc', 4, 27.034, 0.077, 0.069, 28.581, 0.118, 0.126, '--', 'SB2'],
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
    ['KYBoo', 10, 0.510, 0.184, 0.205, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['MPBoo', 5, 0.276, 0.079, 0.079, nan, nan, nan, 'RS', 'SB1'],
    ['MSBoo', 5, 18.595, 0.140, 0.199, nan, nan, nan, 'SSG', 'SB1'],
    ['MXBoo', 4, 0.262, 0.108, 0.111, nan, nan, nan, 'RS', 'SB1'],
    ['MXDra', 9, 39.695, 0.266, 0.305, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['NPCom', 6, 61.883, 0.431, 0.413, nan, nan, nan, 'SSG', 'SB1'],
    ['NRLib', 8, 50.539, 1.143, 1.246, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['OOBoo', 4, 30.136, 0.137, 0.117, 22.988, 0.377, 0.340, 'SSG', 'SB2'],
    ['QZBoo', 12, 32.924, 0.146, 0.190, nan, nan, nan, 'SSG', 'SB1'],
    ['ROTSE1J135649.21+242927.1', 6, 23.842, 0.109, 0.092, 2.318, 0.113, 0.118, 'Cut', 'SB2'],
    ['ROTSE1J160518.19+372623.5', 4, 46.531, 0.138, 0.132, 45.682, 0.348, 0.369, 'RS', 'SB2'],
    ['ROTSE1J162157.22+381733.6', 4, 27.572, 0.103, 0.099, nan, nan, nan, 'RS', 'SB1'],
    ['ROTSE1J170537.88+335052.5', 4, 20.952, 0.274, 0.282, 16.809, 1.050, 0.991, 'RS', 'SB2'],
    ['ROTSE1J172404.95+182937.3', 3, 18.915, 0.155, 0.160, nan, nan, nan, 'SSG', 'SB1'],
    ['ROTSE1J183703.58+280017.1', 5, 12.077, 0.226, 0.246, nan, nan, nan, 'RS', 'Asymmetric'],
    ['ROTSE1J184706.07+434032.7', 5, 7.540, 0.127, 0.137, nan, nan, nan, 'RS', 'SB1'],
    ['RXJ1112.4+3950', 15, 51.708, 0.774, 0.626, nan, nan, nan, 'SSG', 'Asymmetric'],
    ['TYC1339-572-1', 4, 12.046, 0.169, 0.178, 8.673, 0.131, 0.121, 'Cut', 'SB2'],
    ['TYC1543-131-1', 3, 7.990, 0.086, 0.088, nan, nan, nan, 'RS', 'SB1'],
    ['TYC2268-394-1', 11, 58.238, 1.174, 1.262, nan, nan, nan, 'SSG', 'SB1'],
    ['TYC2645-156-1', 6, 25.825, 0.236, 0.250, nan, nan, nan, 'RS', 'Asymmetric'],
    ['TYC3401-730-1', 3, 62.355, 1.527, 1.685, nan, nan, nan, '--', 'SB1'],
    ['TYC3557-919-1', 10, 26.490, 0.129, 0.127, nan, nan, nan, 'SSG', 'SB1'],
    ['TYC3927-1130-1', 5, 35.127, 0.386, 0.328, nan, nan, nan, 'RS', 'SB1'],
    ['TYC4216-800-1', 5, 15.280, 0.145, 0.116, nan, nan, nan, 'RS', 'SB1'],
    ['TYC4223-399-1', 23, 26.459, 1.976, 2.145, nan, nan, nan, 'SSG', 'SB1'],
    ['TYC4415-328-1', 10, 12.916, 0.154, 0.212, nan, nan, nan, 'SSG', 'SB1'],
    ['TYC4422-1970-1', 18, 38.306, 0.552, 0.511, nan, nan, nan, 'RS', 'SB1'],
    ['TYC4428-785-1', 3, 9.378, 0.106, 0.091, nan, nan, nan, 'RS', 'SB1'],
    ['TYC4667-90-1', 11, 52.485, 0.344, 0.450, nan, nan, nan, 'Cut', 'Asymmetric'],
    ['TYC5003-138-1', 8, 46.438, 0.246, 0.288, 60.411, 0.286, 0.297, 'SSG', 'SB2'],
    ['TYC583-566-1', 11, 51.934, 0.516, 0.503, nan, nan, nan, 'SSG', 'SB1'],
    ['UCAC228842571', 35, 45.614, 0.509, 1.075, 46.522, 7.848, 17.727, 'RS', 'SB2'],
    ['V424Gem', 4, 28.076, 0.091, 0.106, nan, nan, nan, 'SSG', 'SB1'],
    ['V457Vir', 6, 47.867, 0.446, 0.405, nan, nan, nan, 'RS', 'Asymmetric'],
    ['V498And', 11, 30.630, 0.280, 0.373, nan, nan, nan, 'SSG', 'SB1'],
    ['V503Hya', 5, 70.636, 0.647, 0.578, 54.455, 0.904, 0.902, '--', 'SB2'],
    ['V832Her', 5, 10.286, 0.059, 0.058, nan, nan, nan, '--', 'SB1'],
    ['V834Her', 7, 1.159, 0.165, 0.194, nan, nan, nan, '--', 'SB1'],
    ['V846Her', 5, 17.532, 0.206, 0.262, nan, nan, nan, '--', 'SB1']
    ]
    ##################
    ##################



    # Build primary sample & amp lists (SB1 + SB2)
    sampnums     = [row[1] for row in data]
    amps         = [row[2] for row in data]
    amps_hi      = [row[3] for row in data]
    amps_lo      = [row[4] for row in data]
    
    # Identify SB2 systems (those with RV2_amp not NaN)
    sampnums_sb2 = [row[1] for row in data if not np.isnan(row[5])]
    amps_sb2     = [row[2] for row in data if not np.isnan(row[5])]
    
    # Combine into full lists
    sampnums_all = sampnums
    amps_all     = amps

    if run_number == 1:
        print(f'Total number of Star observations: {len(amps)}')
        print(f'Single peak stars (SB1): {len(amps) - len(amps_sb2)}')
        print(f'Double peak stars (SB2): {len(amps_sb2)}\n')

    # And correspondingly extend hi/lo for the "all" set (assumes SB2 amps use RV2_hi/lo if present)
    # First copy the SB1 hi/lo
    amps_hi_all = amps_hi.copy()
    amps_lo_all = amps_lo.copy()
    # Then append RV2 errors for SB2 entries:
    amps_hi_all += [row[6] for row in data if not np.isnan(row[5])]
    amps_lo_all += [row[7] for row in data if not np.isnan(row[5])]
    
    # Filtered set (by sample limit or amp>50)
    samplelim = limit_on_samples  # user-defined
    ampnew, sampnew, hi_new, lo_new = [], [], [], []
    for s,a,hi,lo in zip(sampnums_all, amps_all, amps_hi_all, amps_lo_all):
        if s >= samplelim or a > 50:
            ampnew.append(a)
            sampnew.append(s)
            hi_new.append(hi)
            lo_new.append(lo)
    
    # Compute CDFs
    # All observations
    vals1 = np.array(sorted(amps_all))
    lo1   = np.array([lo for _,lo in sorted(zip(amps_all, amps_lo_all))])
    hi1   = np.array([hi for _,hi in sorted(zip(amps_all, amps_hi_all))])
    p1    = np.arange(len(vals1)) / (len(vals1) - 1)
    
    # Filtered observations
    vals2 = np.array(sorted(ampnew))
    lo2   = np.array([lo for _,lo in sorted(zip(ampnew, lo_new))])
    hi2   = np.array([hi for _,hi in sorted(zip(ampnew, hi_new))])
    p2    = np.arange(len(vals2)) / (len(vals2) - 1)
    



    if dispCDF == True:
        # CDF curves
        #plt.step(vals1, p1*100, where='post', label='All Obs. CDF', linewidth=2)
        #plt.step(vals2, p2*100, where='post', label=f'Filtered (≥{samplelim} samples or amp>50)', linewidth=2)
        
        #plt.step(vals1, p1*100, where='post', linewidth=2, linestyle = '--', color = 'black', alpha = 0.5)
        #plt.step(vals2, p2*100, where='post', linewidth=2, color = 'black', alpha = 0.5)
        
        # Errorbars on all observations
        ax.errorbar(vals1, p1*100,
                    xerr=[lo1, hi1],
                    color = 'black', fmt='x',
                    ecolor='grey', alpha=1, label=f"Observations with 3+ Measurements",
                    markeredgewidth=2)
                
                # Errorbars on filtered observations
        ax.errorbar(vals2, p2*100,
                    xerr=[lo2, hi2],
                    color = 'black', fmt='.', 
                    ecolor='grey', alpha=1, label=f"Observations with {samplelim}+ Measurements",
                    markersize = 8)
                
        
        
        
        # Highlight SB2
        '''sb2_mask1 = np.isin(vals1, amps_sb2)
        ax.scatter(vals1[sb2_mask1], p1[sb2_mask1]*100,
                    s=20, facecolors='orange', label='SB2 (All)',
                    zorder = 10, marker = 'x')
        
        sb2_mask2 = np.isin(vals2, amps_sb2)
        ax.scatter(vals2[sb2_mask2], p2[sb2_mask2]*100, zorder = 10, marker = '.',
                    s=25, facecolors='orange')'''
            
        ##################
        ##################
        
        
        
        
        
        
        ax.set_xlabel('K (km/s)')
        ax.set_ylabel('Probability (%)')
        #plt.title(f'Run {run_number}: CDF of Simulated and Actual Observation Methods')
        
        ax.set_xlim(-5, 105)
        ax.set_ylim(-5, 105)
        #plt.xlim(105, -5)
        #plt.ylim(105, -5)
        ax.legend(loc='lower right', fontsize=10)
        #ax.legend()
        plt.grid(alpha = 0.5, linestyle = "--")
        plt.tight_layout()
        #plt.savefig("Sim_Full_Eccens_No_Fill.png", dpi = 500)
        #plt.savefig("CDF1.png", dpi = 300)
        
        if SaveThisFig:
            plt.savefig(savefig_name + ".png", dpi = 300)

        #if show:
        #    plt.show()

        
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
        for actual, measured, label, color in zip([actualinfe, actual2e, actual1e, actual3], [measuredinfe, measured2e, measured1e, measured3], 
                                                  [r'Eccentric, Lognormal (2-$10^{10}$ days)', 'Eccentric, Lognormal (2-1000 days)', 
                                                   'Eccentric, Lognormal (2-30 days)', 'Tidally Locked (Leiner SSGs Distribution)'],
                                                  ['green', 'red', 'blue', 'y']):
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

    #plt.close()


    #Return KS and p values
    return mKSval, aKSval, mpval, apval
