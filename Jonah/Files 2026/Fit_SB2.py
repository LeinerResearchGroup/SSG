
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter

import saphires as saph


def Run_SB2(tar,tar_spec, vel,bf,bf_sterr,bf_wstd, filename, star, date_title, BVCORR,  p0):

  R_coude = 60000.0     #spectral resolution of the Coude.
  c = (2.9979245*10**5)


 
        
  n_boot = 500 #number of bootstrap interations to perform
  rv_boot_dist = np.zeros(n_boot) #a blank array of RV values to fill
  bf_boot_dist = np.zeros([n_boot,tar_spec[tar[0]]['vel'].size]) #an array to save the BFs we create

  rv_boot_dist_primary = np.zeros(n_boot)
  rv_boot_dist_secondary = np.zeros(n_boot)
    
  rv_amp_primary = np.zeros(n_boot)
  rv_amp_secondary = np.zeros(n_boot)
 
  for i in range(n_boot):
      rindwr = np.random.randint(0,tar.size,tar.size) #randomly sample the contributing orders with replacement.
      tar_i = tar[rindwr] #pick out the selected orders
      vel_i,bf_i,bf_ste_i,bf_wstd_i = saph.bf.weight_combine(tar_i,tar_spec,vel_gt_lt=(+150,-150)) #create a version of the combined BF
      fit_i = (curve_fit(saph.utils.d_gaussian_off,vel_i,bf_i,p0=p0, maxfev=10000))[0] #fit it

          # Identify which is the primary peak (stronger amplitude)

      rv_boot_dist_primary[i] = fit_i[1]  # mu1
      rv_boot_dist_secondary[i] = fit_i[4]  # mu2

      
      #bf_boot_dist_primary[i,:] = bf_i
      #bf_boot_dist_secondary[i,:] = bf_i
      rv_amp_primary[i] = fit_i[0]  # A1
      rv_amp_secondary[i] = fit_i[3]  # A2
    
    
  rv_boot_16_primary,rv_boot_med_primary,rv_boot_84_primary = np.percentile(rv_boot_dist_primary,[16,50,84])
  rv_booterr_lo_primary = rv_boot_med_primary - rv_boot_16_primary
  rv_booterr_hi_primary = rv_boot_84_primary - rv_boot_med_primary


  rv_boot_16_secondary,rv_boot_med_secondary,rv_boot_84_secondary = np.percentile(rv_boot_dist_secondary,[16,50,84])
  rv_booterr_lo_secondary = rv_boot_med_secondary - rv_boot_16_secondary
  rv_booterr_hi_secondary = rv_boot_84_secondary - rv_boot_med_secondary

  ### Last step, convert to barcentric RV. This has already been computed and is in the starname_coude_head.dat file.
  brv_boot_primary = (rv_boot_med_primary) + BVCORR + ((rv_boot_med_primary) * BVCORR / c)
  brv_boot_secondary = (rv_boot_med_secondary) + BVCORR + ((rv_boot_med_secondary) * BVCORR / c)
    
  ### Last step, convert to barcentric RV. This has already been computed and is in the starname_coude_head.dat file.
  A1_mean = np.mean(rv_amp_primary)
  A2_mean = np.mean(rv_amp_secondary)   
  ######
  #adds RV uncorrceted and corrected values to appropriate arrays
  #RVuncor[k].extend([np.round(rv_boot_med,2), np.round(rv_booterr_hi,2), -np.round(rv_booterr_lo,2)])
  #RVcor[k].extend([np.round(brv_boot,2), np.round(rv_booterr_hi,2), -np.round(rv_booterr_lo,2)])
    
    ######
    
    
      #this function creates the rotationally broadened line profile function for a given set of additonal broadening sources.
      
  def make_double_rot_pro(R, micro_turb, macro_turb, R_smooth, R_syn):
      
    '''
    Returns a function that models the sum of two rotational broadening profiles
    with the same convolution kernel as make_rot_pro did
    '''
    if R > 0:
        FWHM = (2.997924e5)/R
        sig_R = FWHM/(2.0*np.sqrt(2.0*np.log(2.0)))
    else:
        sig_R = 0

    if R_smooth > 0:
        FWHM_Rs = (2.997924e5)/R_smooth
        sig_Rs = FWHM_Rs/(2.0*np.sqrt(2.0*np.log(2.0)))
    else:
        sig_Rs = 0

    if R_syn > 0:
        FWHM_Rsyn = (2.997924e5)/R_syn
        sig_Rsyn = FWHM_Rsyn/(2.0*np.sqrt(2.0*np.log(2.0)))
    else:
        sig_Rsyn = 0

    sig = np.sqrt(sig_R**2 + sig_Rs**2 + micro_turb**2 + macro_turb**2 - sig_Rsyn**2)

    def double_rot_pro(x, A1, rv1, rvw1, A2, rv2, rvw2, o):
        '''
        Two-component rotational broadening model with an offset.
        '''
        e = 0.75  # Limb-darkening coefficient

        def single_profile(A, rv, rvw):
            c1 = (2*(1-e)) / (np.pi*rvw*(1-e/3.0))
            c2 = e / (2*rvw*(1-e/3.0))
            core = 1 - ((x - rv)/rvw)**2
            prof = A * (c1 * np.sqrt(core) + c2 * core)
            prof[np.isnan(prof)] = 0.0
            return prof

        prof1 = single_profile(A1, rv1, rvw1)
        prof2 = single_profile(A2, rv2, rvw2)
        total = prof1 + prof2 + o

        v_spacing = x[1] - x[0]
        smooth_sigma = sig / v_spacing

        return gaussian_filter(total, sigma=smooth_sigma)

    return double_rot_pro
    
  #then you will use curve_fit to fit a rotational broadening profile to the broadening function.
  #You should not need to change anything in this cell for different observations
  rot_prof = make_double_rot_pro(R_coude, 0, 0, R_coude, 0)

  p0[1] = rv_boot_med_primary
  p0[4] = rv_boot_med_secondary
  print(f'Vsini p0: {p0}')
  fit_rp,covar_rp = curve_fit(rot_prof,vel,bf,p0=p0, maxfev = 10000)

  '''
  rot_prof = make_double_rot_pro(R_coude, 0, 0, R_coude, 0)
  # --- restrict fit to +/- 100 km/s around each RV component ---
  window_half_width = 50  # km/s
  # center windows on the bootstrap medians for primary & secondary
  mask = (
      (np.abs(vel - rv_boot_med_primary)   <= window_half_width) |
      (np.abs(vel - rv_boot_med_secondary) <= window_half_width)
  )
  vel_fit = vel[mask]
  bf_fit  = bf[mask]
  # Fit only within the union of those windows
  fit_rp, covar_rp = curve_fit(rot_prof, vel_fit, bf_fit, p0=p0, maxfev=10000)
  '''


    
  #print(fit_rp)
  '''
  ax.plot(vel,rot_prof(vel,*fit_rp), '--', lw= 2, label= "Model")
  print('RV (uncorrected):',np.round(rv_boot_med,2),'+',np.round(rv_booterr_hi,2),'-',np.round(rv_booterr_lo,2))
  print('RV (barycenter corrected):',np.round(brv_boot,2),'+',np.round(rv_booterr_hi,2),'-',np.round(rv_booterr_lo,2))
  #print('RV:',np.round(rv_boot_med,2),'+',np.round(rv_booterr_hi,2),'-',np.round(rv_booterr_lo,2), 'km/s')
  print('vsini:',np.round(vsini1_lsfit,2),'+/-',np.round(vsini1err_lsfit,2), 'km/s')
  '''
  # Plot full fitted model

  #fig,ax = plt.subplots(1)
  #for i in range(tar.size):
  #    ax.plot(tar_spec[tar[i]]['vel'],tar_spec[tar[i]]['bf_smooth'], color = "grey", alpha = 0.25)

  fig,ax = plt.subplots(1)
  for key, spec in tar_spec.items():
      
      # Skip if missing BF products
      if 'vel' not in spec or 'bf_smooth' not in spec:
          print(f"Skipping {key}: missing BF")
          continue
      ax.plot(spec['vel'], spec['bf_smooth'], color="grey", alpha=0.5)
  ax.plot(vel,bf, lw=3)#, label= "Observed")
  ax.plot(vel, rot_prof(vel, *fit_rp), '--', lw=2)#, label="Model")

  ax.axvline(rv_boot_med_primary, ls = "--", color = 'red', label = f'RV: {rv_boot_med_primary:.2f} km/s | BF Amp: {A1_mean:.3f}')
  ax.axvline(rv_boot_med_secondary, ls = "--", color = 'green', label = f'RV: {rv_boot_med_secondary:.2f} km/s | BF Amp: {A2_mean:.3f}')

  plt.legend()
  ax.set_xlabel('RV (km/s)')
  ax.set_ylabel('Broadening Function')

  ax.set_title(f"{star} {date_title}")
  plt.xlim(-250, 250)
  plt.ylim(-0.015, 0.08)
  #plt.show()




    
  
  #RV Median (Bootstrap):        {rv_boot_med:.2f} +{rv_booterr_hi:.2f} / -{rv_booterr_lo:.2f} km/s
  #RV (barycenter corrected):    {brv_boot:.2f} +{rv_booterr_hi:.2f} / -{rv_booterr_lo:.2f} km/s''')


  #print(fit_rp)
  A1, rv1, vsini1, A2, rv2, vsini2, offset = fit_rp

  # Ensure stronger peak is labeled "primary"
  if A2 > A1:
      A1, rv1, vsini1, A2, rv2, vsini2 = A2, rv2, vsini2, A1, rv1, vsini1
      
  if A2 > A1:
      vsini1_err, vsini2_err = np.sqrt(covar_rp[5, 5]), np.sqrt(covar_rp[2, 2])
  else:
      vsini1_err, vsini2_err = np.sqrt(covar_rp[2, 2]), np.sqrt(covar_rp[5, 5])

  print(f"""  BVCORR: {BVCORR} km/s
  
       Results
  Primary RV (uncorrected): {np.round(rv_boot_med_primary, 2)} + {np.round(rv_booterr_hi_primary, 2)} - {np.round(rv_booterr_lo_primary, 2)}
  Primary RV (Barycenter Corrected): {np.round(brv_boot_primary, 2)} + {np.round(rv_booterr_hi_primary, 2)} - {np.round(rv_booterr_lo_primary, 2)}

  Secondary RV (uncorrected): {np.round(rv_boot_med_secondary, 2)} + {np.round(rv_booterr_hi_secondary, 2)} - {np.round(rv_booterr_lo_secondary, 2)}
  Secondary RV (Barycenter Corrected): {np.round(brv_boot_secondary, 2)} + {np.round(rv_booterr_hi_secondary, 2)} - {np.round(rv_booterr_lo_secondary, 2)}

  Primary vsini:    {vsini1:.2f} +/- {vsini1_err:.2f} km/s
  Secondary vsini:  {vsini2:.2f} +/- {vsini2_err:.2f} km/s
  """)


  return brv_boot_primary, rv_booterr_hi_primary, rv_booterr_lo_primary, vsini1, vsini1_err, brv_boot_secondary, rv_booterr_hi_secondary, rv_booterr_lo_secondary, vsini2, vsini2_err