PRO lognormal_photochem_fractal, noplots=noplots, nosavecloud=nosavecloud
;gna Dec 2013 (Fractals first try: gna Feb 4 2014)

;Purpose: Take in calculated values from the photochemical model
;and use reasonable assumptions to turn the monomodal photochem
;distribution into a lognormal distribution for spherical particles
;or keep as monomodal for fractal particles to input into SMART.
;We want the number density of a mixture of 2 modes on either size of
;the radius size bin such that we conserve the mass density
;that is found in the layer by the photochemical model.  i.e. Mass density
;per layer is conserved between the photochem output and this
;program's output.
;
;Relevant output are the arrays called numarray and tauarray, which
;contain the number density and optical depth for each particle mode
;in each layer.
;
;Currently, this program requires mie_single.pro, which is a freely
;available IDL routine written by an Earth science group at Oxford.
; Mie_single reads in the particle size parameter (2*pi*r)/lambda and
;the complex index of refraction (currently using Titan tholins at
;1um).  The output we need is Qext.  (for fractal particles, scatting
;parameters we use are derived from files sent by Eric Wolf to Shawn D.G.)
;
;Set noplots flag to suppress plots.  Set nosavecloud flag if you
;don't want to save a new cloud file.
;
;In the future, we may need to figure out how to do lognormal for fractals

;---------------------------------------------------------
;how many and which aerosol files are we going to do?

  readdir = '/astro/users/giada/archean_inputs/sulfur_biosig/Jan13/'
  readcol, readdir+'list_files.lst', fname, format='(a)'

  for yy=0, n_elements(fname) -1 DO BEGIN

;--------------------------------------------------------

     readcol, readdir+fname[yy], zz, aerosol, rpar, wfall, taused, tauedd, tauc, conver

     numlayers = n_elements(zz)
     tauarray = fltarr(numlayers, 5) ;to save optical depths in for the particle modes
     numarray = fltarr(numlayers, 5) ;to save particle densities in for the particle modes
     width = fltarr(numlayers)       ;array for width of layers for calculating optical depth


;----define values & arrays to be used:

     um_per_cm = 1.e4
     g_per_kg = 1000.
     cm_per_m = 100.                           
     ug_per_g = 1.e6
     cm3_per_m3 = 1.e6

     s= 1.5                     ;dimensionless geometric standard deviation.
                                ;Literature says 1.5 is a reasonable
                                ;value for a realistic aerosol
                                ;population.   We might ask Dave what
                                ;he thinks he should use.  What I
                                ;read suggests it should be between
                                ;1.5 and 2.5
     sigma = alog(s)
     


     aerosol = aerosol / ug_per_g / cm3_per_m3 ;I'm still not sure what the units of "aerosol" are, 
                                ;but assuming they're in ug/m^3
                                ;is giving me values that seem reasonable
     



     nindex = 1.65              ;tholin real index of refraction at 1 um (Khare and Sagan 1984) (1um = arbitrary reference wavelength)
     kindex = -1e-3             ;tholin imaginary index of refraction at 1 um (mie_single wants <0 for absorption)
     index = complex(nindex,kindex) 

     rho = 0.84                 ;density of diacetylene (tholins) [g/cm^3]
                                ;Note: Trainer et al (2006) has
                                ;suggested 0.64 g/cm^3 for early earth
                                ;hazes, but the photochemical model
                                ;hazes are diacetylene particles

     rmodes = [0.001, 0.01, 0.1, 1.0, 10.0]/um_per_cm ;log spacing radii (Wolf and Toon call it 'rpar')
     
     rmon = 0.05                ;um (monomer radius as defined by Wolf and Toon 2010 code))
     rmon = rmon / um_per_cm

;these are from the files Shawn got from Eric Wolf:
     filenum = [10, 19, 28, 37, 46]     ;For reference: RT file that corresponds to our chosen rmodes ('rpar')
     df = [3., 3., 1.5, 2.4, 2.4]       ;fractal parameter for each rmode (3 = spherical, 1 = linear chain)
     nmon = [1., 1., 8., 8000., 8.e6]   ;number of monomers in each rmode fractal particle
     Qext_ = [0.000013, 0.000138, 0.032373, 4.747904, 39.599220] ;Qext for fractals of rmode sizes (at 1um)
     w_ = [0.000042, 0.040516, 0.956332, 0.997023, 0.996430]     ;single scatt albedo for fractals of rmode sizes (at 1um)
     g_ = [0.000008, 0.000846, 0.248075, 0.909008, 0.976528]     ;asymm param for fractals of rmode sizes (at 1um)

;for fractal projected area:
     ka = 1.15      ;fractal area constant (Koylu et al 1995)
     alpha = 1.10   ;fractal area exponent (Koylu et al 1995)

;------------------------------------------------------
                                ;width of the layers:
     for q=0, numlayers -2 do width[q] = zz[q+1] - zz[q]
     width[numlayers-1] = zz[numlayers-1] - zz[numlayers -2] ;in cm

     for xx=0, numlayers -1 do begin ;cycle through layers and do calculations
        M = aerosol[xx]              ;g/cm**3 in layer
        

        rphot = rpar[xx]        ;already in cm
        z = width[xx]           ;already in cm

;----------------------------------------------------------------
;DEAL WITH BINS ON EITHER SIDE OF RPHOT:

;first thing is to determine which rmode bins rphot falls between
        binloc_lower = max(where(rmodes LT rphot))
        binloc_upper = min(where(rmodes GT rphot))
;smaller radius : 
        r1 = rmodes[binloc_lower]
        r_lim1_1 = r1/10. ;lower limit of integration for lognormal (I'm also using these limits in miescat.  please don't alter unless you also redo the .mie and .mom files that smart uses)
        r_lim2_1 = r1*10  ;upper limit of integration for lognormal 
        spacing1 = r1/100.      ;this can be played with. (stepsize for integration)
        radii1 = arrgen(r_lim1_1, r_lim2_1, spacing1)

;larger radius :
        r2 = rmodes[binloc_upper]
        r_lim1_2 = r2/10.       ;lower limit of integration for lognormal 
        r_lim2_2 = r2*10        ;upper limit of integration for lognormal 
        spacing2 = r2/100.      ;this can be played with. (stepsize for integration)
        radii2 = arrgen(r_lim1_2, r_lim2_2, spacing2)

;find distances in log space between nearest neighbor bins:
;note: this is how we want to do this, right?
        dist_to_top = alog10(r2) - alog10(rphot)
        dist_to_lower = alog10(rphot) - alog10(r1)
;here is how we will relate N1 and N2 (where N is the number density
;of particles #/cm^3)
        Nfactor = dist_to_lower/dist_to_top ; N2 = Nfactor*N1 
;e.g. if r1 = 1 and r2 = 10 and rphot =3 then in log space,
;we're about halfway between r1 and r2.  So then Nfactor ~1
;and we have roughly equal amounts of N2 and N1.

;saving these variables to make them available to functions I will use
        save, M, rphot, rmodes, sigma, rho, r1, r2, Z, Nfactor, index, rmon, df, nmon, filename='~/idl/photochem_input.sav'

;---------------------------------------------------------------
;FIND N1 AND N2

;now we want to determine the total number of particles with
;characteristic radii r1 and r2  required 
;to get the same mass here as in the photochemical 
;The Wolf & Toon (2010) RT files we were given are almost certainly
;monomodal distributions.  Therefore, we do not use a lognormal
;distribution as we did for the spherical particle case EXCEPT where
;df = 3 (then the particle IS spherical)

;Calculate number density for r1 and r2 distributions (N_r1 and N_r2)


        tot_inte1 = 0.
        tot_inte2 = 0.

;what we are doing:
; M = N1 * [integral(n(r1)) + Nfactor*integral(n(r2))]
        
;NOTE: for fractal particles, total mass is calculated by mass of the
;number of monomers in a particle.

;integration loop (Done really simplistically; could be improved if
;needed):

;to fix a weird issue:
num = [n_elements(radii1), n_elements(radii2)]
minnum = min(num)

        IF r2 LE 0.01 THEN BEGIN ;if both bins are in spherical regime
        for i=0, minnum-2 do begin
           
                                ;integrate
                                ;functions (e.g. mass_r1_n1) found below
           
           inte1 = spacing1 * (mass_r1_n1(radii1[i]) + mass_r1_n1(radii1[i+1]))/2.
           inte2 = spacing2 * (mass_r2_n1(radii2[i]) + mass_r2_n1(radii2[i+1]))/2.

           tot_inte1 = tot_inte1 + inte1
           tot_inte2 = tot_inte2 + inte2
           

        endfor
     ENDIF

        IF r2 EQ 0.1 THEN BEGIN ;if smaller bin is in spherical regime and larger in fractal
        for i=0, minnum-2 do begin
           
                                ;integrate
                                ;functions (e.g. mass_r1_n1) found below
           
           inte1 = spacing1 * (mass_r1_n1(radii1[i]) + mass_r1_n1(radii1[i+1]))/2.

           tot_inte1 = tot_inte1 + inte1
        endfor
        tot_inte2 = (4./3.)*!pi*rho*Nfactor*nmon[binloc_upper]*rmon^3.   ;total mass for r2 (monomodal)
     ENDIF

        IF r1 GE 0.1 THEN BEGIN ;if both bins are in fractal regime
        tot_inte1 = (4./3.)*!pi*rho* nmon[binloc_lower]*rmon^3. ;total mass for r1 (monomodal)
        tot_inte2 = (4./3.)*!pi*rho*Nfactor*nmon[binloc_upper]*rmon^3.   ;total mass for r2 (monomodal)
     ENDIF


        N_r1 = M/(tot_inte1+tot_inte2) ; number density of r1 particles
        N_r2 = Nfactor * N_r1          ; number density of r2 particles

        save, n_r1, n_r2,filename='~/idl/photochem_number_densities.sav'



;check with IDL's built in integrator
                                ; print, 'mass is: '
                                ; print, qsimp('mass_r1_n1', r_lim1, r_lim2) + qsimp('mass_r2_n2', r_lim1, r_lim2)
                                ;(yes, the mass out = the mass in)

;----------------------------------------------------------
; OPTICAL DEPTH:

;what we ultimately want for the cloud files is the total OPTICAL
;DEPTH in the atmospheric layer from each mode: tau = (cross sectional
;area) * (number density) * (thickness of layer) *
;(wavelength-dependent extinction efficiency)
;particles uniformly distributed in each layer.
;For the fractal particles, we multiply the optical depth by a factor
;involving the asymmetry param and single-scattering albedo to produce
;the effective optical depth (Wolf and Toon 2010 supplemental
;material).
        
        tot_tau1 = 0.
        tot_tau2 = 0.
;find the cross sectional area contribution from each radius in
;the mode's lognormal distribution of radii:
;note that the area_ functions used below require mie_single.pro
      IF r2 LE 0.01 THEN BEGIN ;if both bins are in spherical regime
        for i=0, minnum-2 do begin
           
                                ;integrate
           inte1 = spacing1 * (area_r1_n1(radii1[i]) + area_r1_n1(radii1[i+1]))/2.
           inte2 = spacing2 * (area_r2_n2(radii2[i]) + area_r2_n2(radii2[i+1]))/2.
           
           tot_tau1 = tot_tau1 + inte1
           tot_tau2 = tot_tau2 + inte2

        endfor
     ENDIF

   IF r2 EQ 0.1 THEN BEGIN ;if smaller bin is in spherical regime and larger in fractal
        for i=0, minnum-2 do begin
           
                                ;integrate
           inte1 = spacing1 * (area_r1_n1(radii1[i]) + area_r1_n1(radii1[i+1]))/2.
                    
           tot_tau1 = tot_tau1 + inte1
          
        endfor
        tot_tau2 = [1-w_[binloc_upper]/2. - (g_[binloc_upper]*w_[binloc_upper])/2.] * Qext_[binloc_upper] * n_r2 * rpar^2 * !pi * (nmon[binloc_upper]/ka)^(1/alpha) 
     ENDIF

   IF r1 GE 0.1 THEN BEGIN ;if both bins are in fractal regime
        tot_tau1 = [1-w_[binloc_lower]/2. - (g_[binloc_lower]*w_[binloc_lower])/2.] * Qext_[binloc_lower] * n_r1 * rpar^2 * !pi * (nmon[binloc_lower]/ka)^(1/alpha) 
        tot_tau2 = [1-w_[binloc_upper]/2. - (g_[binloc_upper]*w_[binloc_upper])/2.] * Qext_[binloc_upper] * n_r2 * rpar^2 * !pi * (nmon[binloc_upper]/ka)^(1/alpha) 
     ENDIF


        tau1 = tot_tau1 * Z
        tau2 = tot_tau2 * Z


;now save to our tauarray and numarray:
        tauarray[xx, binloc_lower] = tau1
        tauarray[xx, binloc_upper] = tau2

        numarray[xx, binloc_lower] = N_r1
        numarray[xx, binloc_upper] = N_r2

        print, 'end of layer number ' + strtrim(xx,1)
     endfor                     ;cycle through layers


;------------------------------------------------------
;CLOUD FILE

;last step is to write a cloud file:
     heights = zz / 1e5
     save, tauarray, heights, numarray, aerosol, rpar, fname, filename="~/idl/cldtauinput.sav"
     IF not keyword_set(nosavecloud) THEN write_cloud_file, tauarray, heights, fname[yy], fractal=1


;--------------------------------------------------------
;MAKE PLOTS
     savepathplots = '/astro/users/giada/archean_earth/'+fname[yy]

     IF not keyword_set(noplots) THEN BEGIN
        restore, filename="~/idl/cldtauinput.sav"
        set_colors
        !p.multi=0
        set_plot, 'ps'
        device, file=savepathplots+'fractal_number_density.eps', /encapsulated, /color
        plot, heights, numarray[*,0], color=0, /ylog, yrange=[1e-3, 1e12], xtitle='height [km]', ytitle='particles/cm**3', title='calculated number densities'
        oplot, heights, numarray[*,1], color=50
        oplot, heights, numarray[*,2], color=100
        oplot, heights, numarray[*,3], color=200
        oplot, heights, numarray[*,4], color=230
        legend, ['0.001', '0.01', '0.1', '1.0', '10.0'], color=[0,50,100,200,230], textcolor=[0,0,0,0,0], linestyle=[0,0,0,0,0]
        device, /close
        set_plot, 'x'

        set_plot, 'ps'
        device, file=savepathplots+'fractal_optical_depth.eps', /encapsulated, /color
        plot, heights, tauarray[*,0], color=0, /ylog, yrange=[1e-9, 1e5], xtitle='height [km]', ytitle='tau', title='calculated optical depths'
        oplot, heights, tauarray[*,1], color=50
        oplot, heights, tauarray[*,2], color=100
        oplot, heights, tauarray[*,3], color=200
        oplot, heights, tauarray[*,4], color=230
        legend, ['0.001', '0.01', '0.1', '1.0', '10.0'], color=[0,50,100,200,230], textcolor=[0,0,0,0,0], linestyle=[0,0,0,0,0]
        device, /close
        set_plot, 'x'

        !p.multi = [0,2,1]
        set_plot, 'ps'
        device, file=savepathplots+'fractal_aerosol_rpar.eps', /encapsulated, /color
        plot, heights, aerosol, color=0, xtitle='height [km]', ytitle='g/cm**3', title='photochem mass density', /ylog
        plot, heights, rpar*1e4, color=0, xtitle='height [km]', ytitle='radius [um]', title='particle radii',/ylog
        device, /close
        set_plot, 'x'

        rarr = arrgen(0.001, 10, 0.001)
        dx = (2*!pi*rarr)/(1)   ;1e-4 = um per cm (1 um reference wavelength)
        mie_single, dx, index, Dqext

        !p.multi = 0
        set_plot, 'ps'
        device, file=savepathplots+'fractal_qext.eps', /encapsulated, /color
        plot, rarr, dqext, color=0, /xlog, /ylog, xtitle='particle radius [um]', ytitle='Qext', title='Qext mie scattering'
        device, /close
        set_plot, 'x'

     ENDIF
  ENDFOR                        ;yy (number of aerosol files)
  stop
END








FUNCTION mass_r2_n1, r

  restore, '~/idl/photochem_input.sav'
  mass =  rho * (4./3.)* !pi *r^3.
  aa =  1/sqrt(2*!pi) * 1/sigma 
  bb = 2*sigma^2

  mass =   Nfactor * 1/r * aa * mass * exp(-((alog(r) - alog(r2))^2) / bb)

  return, mass

END


FUNCTION mass_r1_n1, r

  restore, '~/idl/photochem_input.sav'
  mass =  rho * (4./3.)* !pi *r^3.
  aa =  1/sqrt(2*!pi) * 1/sigma 
  bb = 2*sigma^2

  mass =  1/r * aa *  mass * exp(-((alog(r) - alog(r1))^2) / bb)

  return, mass

END


FUNCTION area_r1_n1, r
;this calculates the area * number density = B at our reference
;wavelength.  Relies on mie_single, which gives us the Qext,
;extinction efficiency based on our wavelength and the Tholin
;complex refractive index

  restore, '~/idl/photochem_input.sav'
  area = !pi * r^2.
  aa =  1/sqrt(2*!pi) * 1/sigma 
  bb = 2*sigma^2

  restore, '~/idl/photochem_number_densities.sav'
  
  dx = (2*!pi*r)/(1/1e4);1e-4 = um per cm (1 um reference wavelength)
  mie_single, dx, index, Dqext

  B =  Dqext * N_r1 * 1/r * aa * area * exp(-((alog(r) - alog(r1))^2) / bb)

  return, B

END


FUNCTION area_r2_n2, r
;this calculates the area * number density = B at our reference
;wavelength.  Relies on mie_single, which gives us the Qext,
;extinction efficiency based on our wavelength and the Tholin
;complex refractive index

  restore, '~/idl/photochem_input.sav'
  area = !pi * r^2.
  aa =  1/sqrt(2*!pi) * 1/sigma 
  bb = 2*sigma^2

  restore, '~/idl/photochem_number_densities.sav'

  dx = (2*!pi*r)/(1/1e4) ;1e-4 = um per cm (1 um reference wavelength)
  mie_single, dx, index, Dqext
  B =  Dqext * N_r2 * 1/r * aa * area * exp(-((alog(r) - alog(r2))^2 / bb))

  return, B

END


FUNCTION lognormal_function_n1, r
;never actually used this
  restore, '~/idl/photochem_input.sav'

  aa = 1/sqrt(2*!pi) * 1/sigma 
  bb = 2*sigma^2

  restore, '~/idl/photochem_number_densities.sav'

  num =   N_r1 * 1/r * aa * exp(-((alog(r) - alog(r1))^2) / bb)

  return, num

END

FUNCTION lognormal_function_n2, r
;never actually used this
  restore, '~/idl/photochem_input.sav'

  aa = 1/sqrt(2*!pi) * 1/sigma 
  bb = 2*sigma^2

  restore, '~/idl/photochem_number_densities.sav'

  num =   N_r2 * 1/r * aa * exp(-((alog(r) - alog(r2))^2) / bb)

  return, num

END
