import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter

import saphires as saph


def Run_SB1(tar,tar_spec, vel,bf,bf_sterr,bf_wstd, filename, star, date_title, BVCORR,  p0):
  R_coude = 60000.0     #spectral resolution of the Coude.
  c = (2.9979245*10**5)

  n_boot = 500 #number of bootstrap interations to perform
  rv_boot_dist = np.zeros(n_boot) #a blank array of RV values to fill
  bf_boot_dist = np.zeros([n_boot,tar_spec[tar[0]]['vel'].size]) #an array to save the BFs we create

  rv_amp_primary = np.zeros(n_boot)
  #rv_amp_secondary = np.zeros(n_boot)
  for i in range(n_boot):
      rindwr = np.random.randint(0,tar.size,tar.size) #randomly sample the contributing orders with replacement.
      tar_i = tar[rindwr] #pick out the selected orders
      vel_i,bf_i,bf_ste_i,bf_wstd_i = saph.bf.weight_combine(tar_i,tar_spec,vel_gt_lt=(+150,-150)) #create a version of the combined BF
      fit_i = (curve_fit(saph.utils.gaussian_off,vel_i,bf_i,p0=p0, maxfev=10000))[0] #fit it
      rv_boot_dist[i] = fit_i[1] #save the RV
      bf_boot_dist[i,:] = bf_i


    # Store amplitudes
      rv_amp_primary[i] = fit_i[0]  # A1
      #rv_amp_secondary[i] = fit_i[3]  # A2


  rv_boot_16,rv_boot_med,rv_boot_84 = np.percentile(rv_boot_dist,[16,50,84])
  rv_booterr_lo = rv_boot_med - rv_boot_16
  rv_booterr_hi = rv_boot_84 - rv_boot_med

  ### Last step, convert to barcentric RV. This has already been computed and is in the starname_coude_head.dat file.
  brv_boot = (rv_boot_med) + BVCORR + ((rv_boot_med) * BVCORR / c)

  ### Last step, convert to barcentric RV. This has already been computed and is in the starname_coude_head.dat file.
  A1_mean = np.mean(rv_amp_primary)
  #A2_mean = np.mean(rv_amp_secondary)

  # -------------------------------------------------
  # RV single-Gaussian fit (for visualization)
  # -------------------------------------------------
  # Fit gaussian_off to the FINAL combined BF (vel, bf),
  # seeded using the bootstrap median RV so it locks onto the correct peak.
  p0_rv = np.array(p0, dtype=float).copy()
  p0_rv[1] = rv_boot_med  # x0 seed

  fit_g, cov_g = curve_fit(
      saph.utils.gaussian_off,
      vel,
      bf,
      p0=p0_rv,
      maxfev=10000
  )
  bf_gauss_model = saph.utils.gaussian_off(vel, *fit_g)

######
    #adds RV uncorrceted and corrected values to appropriate arrays
  #RVuncor[k].extend([np.round(rv_boot_med,2), np.round(rv_booterr_hi,2), -np.round(rv_booterr_lo,2)])
  #RVcor[k].extend([np.round(brv_boot,2), np.round(rv_booterr_hi,2), -np.round(rv_booterr_lo,2)])

######
    
  #this function creates the rotationally broadened line profile function for a given set of additonal broadening sources.
  print()
  ###
  #You don't need to change anything here
  def make_rot_pro(R,micro_turb,macro_turb,R_smooth,R_syn):
      '''
      A function to make a specific rotationally broadened fitting
      function for a given spectral resolution.
      '''

      if R >0:
          FWHM = (2.997924*10**5)/R
          sig_R = FWHM/(2.0*np.sqrt(2.0*np.log(2.0)))
      else:
          sig_R = 0

      if R_smooth > 0:
          FWHM_Rs = (2.997924*10**5)/R_smooth
          sig_Rs = FWHM_Rs/(2.0*np.sqrt(2.0*np.log(2.0)))
      else:
          sig_Rs = 0

      if R_syn > 0:
          FWHM_Rsyn = (2.997924*10**5)/R_syn
          sig_Rsyn = FWHM_Rsyn/(2.0*np.sqrt(2.0*np.log(2.0)))
      else:
          sig_Rsyn = 0

      sig = np.sqrt(sig_R**2 + sig_Rs**2 + micro_turb**2 + macro_turb**2 - sig_Rsyn**2)

      def rot_pro_ip(x,A,rv,rvw,o):
          '''
          Rotational line broadening function.

          To produce an actual line profile, you have to convolve this function
          with an acutal spectrum.

          In this form it can be fit directly to a the Broadening Fucntion.

          This is in velocity so uf you're going to convolve this with a spectrum make sure
          to take the appropriate cautions, whatever they may be.
          '''
          e = 0.75 # Limb-Darkening Coefficient

          c1 = (2*(1-e))/(np.pi*rvw*(1-e/3.0))
          c2 = e/(2*rvw*(1-e/3.0))

          prof=A*(c1*np.sqrt(1-((x-rv)/rvw)**2)+c2*(1-((x-rv)/rvw)**2))+o

          prof[np.isnan(prof)] = o

          v_spacing = x[1]-x[0]

          smooth_sigma = sig/v_spacing

          prof_conv=gaussian_filter(prof,sigma=smooth_sigma)

          return prof_conv

      return rot_pro_ip

  #then you will use curve_fit to fit a rotational broadening profile to the broadening function.
  #You should not need to change anything in this cell for different observations


    
  rot_prof = make_rot_pro(R_coude,0,0,R_coude,0)
  p0[1] = rv_boot_med
  print(f'Vsini p0: {p0}')

  fit_rp,covar_rp = curve_fit(rot_prof,vel,bf,p0=p0, maxfev = 10000)

  '''
  rot_prof = make_rot_pro(R_coude, 0, 0, R_coude, 0)
  # --- restrict fit region to +/- 100 km/s around the RV peak ---
  window_half_width = 4*p0[2]  # km/s
  print(f'\nvsini fitting width = {window_half_width} km/s\n')
  center = rv_boot_med            # center on the measured RV
  mask = (vel >= center - window_half_width) & (vel <= center + window_half_width)
  vel_fit = vel[mask]
  bf_fit  = bf[mask]
  # Now fit only in that restricted window
  fit_rp, covar_rp = curve_fit(rot_prof, vel_fit, bf_fit, p0=p0, maxfev=10000)
  '''
    
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

  ax.plot(vel,bf, lw=4, label="Observed BF")
  # NEW: RV Gaussian model overlay
  ax.plot(vel, bf_gauss_model, '--', lw=2)#, label="RV Gaussian fit")
  ax.axvline(rv_boot_med, ls = "--", label = f'RV: {rv_boot_med:.2f} km/s | BF Amp: {A1_mean:.3f}')


  #ax.plot(vel,rot_prof(vel,*fit_rp), '--', lw= 2, label="Rotational model")

  plt.legend()
  ax.set_xlabel('RV (km/s)')
  ax.set_ylabel('Broadening Function')

  ax.set_title(f"{star} {date_title}")
  plt.xlim(-250, 250)
  plt.ylim(-0.015, 0.08)
  plt.ylim(-0.015, 0.04)
  #plt.show()
  #print the results of the above fit
  vsini1_lsfit = fit_rp[2]
  vsini1err_lsfit = np.sqrt(covar_rp[2,2]) #note we are taking the error on the vsini to be the sqrt of the diagonal term in the covariance matrix



  '''print('RV (uncorrected):',np.round(rv_boot_med,2),'+',np.round(rv_booterr_hi,2),'-',np.round(rv_booterr_lo,2))
  print('RV (barycenter corrected):',np.round(brv_boot,2),'+',np.round(rv_booterr_hi,2),'-',np.round(rv_booterr_lo,2))
  #print('RV:',np.round(rv_boot_med,2),'+',np.round(rv_booterr_hi,2),'-',np.round(rv_booterr_lo,2), 'km/s')
  print('vsini:',np.round(vsini1_lsfit,2),'+/-',np.round(vsini1err_lsfit,2), 'km/s')'''


  print(f"""  BVCORR: {BVCORR} km/s

      Results
  Primary RV (uncorrected): {np.round(rv_boot_med, 2)} + {np.round(rv_booterr_hi, 2)} - {np.round(rv_booterr_lo, 2)}
  Primary RV (Barycenter Corrected): {np.round(brv_boot, 2)} + {np.round(rv_booterr_hi, 2)} - {np.round(rv_booterr_lo, 2)}

  Primary vsini: {np.round(vsini1_lsfit,2)} +/- {np.round(vsini1err_lsfit,2)} km/s
  """)



  return brv_boot, rv_booterr_hi, rv_booterr_lo, vsini1_lsfit, vsini1err_lsfit, np.nan, np.nan, np.nan, np.nan, np.nan
