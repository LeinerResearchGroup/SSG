The set of files contained in this folder are all the final versions of the major files I used to create and analyze the data for this paper. This ReadMe outlines what each file contains.

---------------------------------------
---------------------------------------
RV_and_Vsini_Code_combined.ipynb
    This is where I started the analysis; creating broadening functions for each star using the saphires package developed for this data. It extracts Radial Velocity (RV) and Vsin(i) values for each observation date based on each star, though it requires each star be looked at individually - each file will be automated, but only one star at a time can be assessed. This is due to the broadening function auto-estimation, which was a bit finnicky at times and needed to be adjusted. Overall, this file extracts RV and Vsini data and plots the variation over time, with options for SB1 and SB2 fitting options. This file saves PNGs of each RV plot and RV/Vsini/Overall summary txt files, and a summary text file, saved to the PNGs and starstxt folders.
    Requires:
        Fit_SB1.py
        Fit_SB2.py
        Saphires
        Star Data
            .dat files
            .p files

PNGs, starstxt, and summary folders
    Contains the RV plots (.png) and all extracted data (.txt) for each star.
---------------------------------------
---------------------------------------
RV_and_Vsini_Code_combined-nres.ipynb
    Same format as RV_and_Vsini_Code_combined, for the samples from nres (seperate from the Coude sample). Same requirements, different save folders.
    Requires:
        Fit_SB1.py
        Fit_SB2.py
        Saphires
        Star Data
            .dat files
            .p files

PNGs_2, starstxt_nres, and summary_nres folders
    Contains the RV plots (.png) and all extracted data (.txt) for each star.
---------------------------------------
---------------------------------------
coude_nres_combine.ipynb
    This file combines the Coude and nres datasets (where applicable) and plots them, saving the plots together. To be used for future work.
---------------------------------------
---------------------------------------
Fit_SB1.py and Fit_SB2.py
    These are the two fitting functions for SB1 and SB2 sources, called in RV_and_Vsini_Code_combined.ipynb. This does the Broadening Function fitting and Vsini extraction based on the parameters fed into the function. 
---------------------------------------
---------------------------------------
RV_Vsini_statistics.ipynb
    This file takes the data generated from the RV and Vsini extraction (RV_and_Vsini_Code_combined.ipynb) and generates statistics for the dataset. Primarily used to confirm which stars have been processed and as confirmation of the Table used in our paper.
    Requries
        starstxt folder (generated data for each star from RV_and_Vsini_Code_combined.ipynb)
---------------------------------------
---------------------------------------
RV_K_investigation.ipynb
    This file is the final version of the full CDF creation. It includes simulated data generation and CDF analysis in one file separated into blocks for *hopefully* a easier step-by-step understanding of the process. This was used as the template for the final SamplingRVs.py file. RV_K_investigation.ipynb goes through the process of generating each distribution, generating orbits, sampling/observing from these simulated orbits, and creating the CDF from the simulations. It then overlays the data obtained from RV_and_Vsini_Code. Finally, the KS test is performed between each simulation and the real dataset. It also includes a few PDFs we had considered for implementation as well. Slight chagnes have been made for the SamplingRVs.py file, but this .ipynb should give a good overview of the main steps for our process.
    Requires:
        Leiner 2022 table (.dat file, included in this folder as table1.dat)
        RV_and_Vsini_Code data (this is actually included in the file, though it might change with further star observation.)
---------------------------------------
---------------------------------------
Sampling_Plots.ipynb
    This file is a test/demonstration of the orbit creation and sampling. It includes a snippet of the code used to create the distributions and then sample from the generated orbit.
    Requires:
        Leiner 2022 table (.dat file, included in this folder as table1)
---------------------------------------
---------------------------------------
SamplingRVs.py
    This .py file uses the same code as RV_K_investigation.ipynb, but formatted as a single function for simplicity in the Run_CDF_Plot.ipynb file.
    Requires:
        Leiner 2022 table (.dat file, included in this folder as table1)
        RV_and_Vsini_Code data (Included in the file, though it might change with further star observation.)
---------------------------------------
---------------------------------------
Run_CDF_Plot.ipynb
    This file is the main run file of our data analysis. It includes the primary function from SamplingRVs.py with options for function execution, including CDF generation, KS test runs, plot saving, sampling and number of run changes, etc. Run this set to get all final values obtained for our analysis.
    Requires
        SamplingRVs.py file
            Leiner 2022 table (.dat file, included in this folder as table1)
            RV_and_Vsini_Code data/array (Included in the SamplingRVs.py file.)
---------------------------------------
---------------------------------------
Run_CDF_Mult_Plot.ipynb
    This file is setup to run/generate the multiple low-sample CDF sets. It is currently set up to run 10 iterations of 100 runs - so during iteration 1, it runs the RV Monte Carlo generation 100 seperate times, each with a low sample count (set to 69, the number of stars in our sample), and runs a KS test on that iteration, saving the resutls. It repeats for the selected number of iterations (10), and in the end takes every iteration's KS test result and compares across sets, as described in the paper. The resulting estimated sigma deviations are a combination of variance within AND across iterations, in the form: Var(X) = E[Var(X|I)] + Var(E[X|I]). So the KS test result from this file accounts for both deviation within each set of 100 runs, as well as across all 10 iterations, resulting in a more accurate estimation of the variance in results of the real observations. Think of it as the difference between observing 10000-100000 random sources straight-up, and observing 100 different sets of 69 random sources 10 seperate times.
    Requires
        SamplingRVs.py file
            Leiner 2022 table (.dat file, included in this folder as table1)
            RV_and_Vsini_Code data/array (Included in the SamplingRVs.py file.)
---------------------------------------
---------------------------------------
Each file should have its own notes in the margins for each block with some explanation as to what each part does. You will have to change most of the file/directory paths before running the files once you have everything downloaded. 
Enjoy!

    Jonah Wilkes
    Illinois Institute of Technology

========== End of Line ==========