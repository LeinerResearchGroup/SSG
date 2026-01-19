The set of files contained in this folder are all the final versions of the major files I used to create and analyze the data for this paper. This ReadMe outlines what each file contains.

---------------------------------------
---------------------------------------
RV_and_Vsini_Code.ipynb
    This is where I started the analysis; creating broadening functions for each star using the saphires package developed for this data. It extracts Radial Velocity (RV) and Vsin(i) values for each observation date based on each star, though it requires each star be looked at individually - each file will be automated, but only one star at a time can be assessed. This is due to the broadening function auto-estimation, which was a bit finnicky at times and needed to be adjusted. Overall, this file extracts RV and Vsini data and plots the variation over time. The outputs for each star are contained in the RV-Plots and starstxt folders.
    Requires:
        Saphires
        Star Data
            .dat files
            .p files


RV-Plots and starstxt folders
    Contains the RV plots (.png) and all extracted data (.txt) for each star.
---------------------------------------
---------------------------------------
RV_K_investigation.ipynb
    This file is the final version of the full CDF creation. It includes simulated data generation and CDF analysis in one file separated into blocks for *hopefully* a easier step-by-step understanding of the process. This was used as the template for the final SamplingRVs.py file. RV_K_investigation.ipynb goes through the process of generating each distribution, generating orbits, sampling/observing from these simulated orbits, and creating the CDF from the simulations. It then overlays the data obtained from RV_and_Vsini_Code. Finally, the KS test is performed between each simulation and the real dataset. It also includes a few PDFs we had considered for implementation as well
    Requires:
        Leiner 2022 table (.dat file, included in this folder as table1)
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
            RV_and_Vsini_Code data (Included in the file, though it might change with further star observation.)
---------------------------------------
---------------------------------------

Each file should have its own notes in the margins for each block with some explanation as to what each part does. You will have to change most of the file/directory paths before running the files once you have everything downloaded. 
Enjoy!

    Jonah Wilkes
    Illinois Institute of Technology