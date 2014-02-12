PRO write_cloud_file, tauarray, heightarray, fname, fractal=fractal
;gna Dec 2013
  
;PURPOSE: read in an array of optical depths for each mode and also an
;array of altitudes.  Outputs a cloud file for each mode.

;Needs 'arrgen.pro', which is easy to find on google.  (want the version
;hosted at the ifa.hawaii.edu site, not the one hosted by
;boulder.swri.edu)

;fractal flag needs to be set to write to fractal directory and not
;overwrite mie clouds in main cloud directory.
;lognormal_photochem_fractal will take care of this if running from there.

;Set number of atmos levels to 55 (so we don't even hit smart's limit):
  heightarray_smaller = arrgen(min(heightarray), max(heightarray), nstep=55)
;Set tau to correspond to that many levels:
  s = size(tauarray)
  tauarr_smaller = fltarr(n_elements(heightarray_smaller),s[2])
  for i=0, s[2] -1 do tauarr_smaller[*,i] = spline(heightarray, tauarray[*,i], heightarray_smaller)

  tauarray_ = tauarr_smaller
  heightarray = heightarray_smaller

  s = size(tauarray_)
  nmodes = s[2]
  modenumber = indgen(nmodes)
  nheights = s[1]

;create mode arrays:
  for mm=0, nmodes-1 do begin
     temparr = tauarray_[*,mm]
     tootiny = where(temparr LT 1e-7) ;smart seems to ignore very small tau
     temparr[tootiny] = 0.0
     IF execute('mode'+strtrim(mm,1)+  $
                '=temparr') EQ 0 then $
                   print, 'bad'
     tauarray_[*,mm] = temparr 
  endfor

  

;define a header for the cloud file:
  header = strarr(1)
  header='alt(km)     mode'

  header0 = strarr(8)

  header0[0] = ' '                                  
  header0[1] = ' Model Atmosphere Aerosol Distribution'
  header0[2] = ' '                                   
  header0[3] = ' Number of aerosol layers = '+strtrim(nheights,1)  
  header0[4] = ' Aerosol reference wavenumber = 10000'
  header0[5] = ''
  header0[6] = ' Homogeneous layer aerosol optical depths'
  header0[7] = ''


  footer = strarr(1)
  footer[0] = ''

  IF keyword_set(fractal) THEN filepath = '/astro/users/giada/archean_inputs/sulfur_biosig/cloud/fractal/'
  IF NOT keyword_set(fractal) THEN filepath = '/astro/users/giada/archean_inputs/sulfur_biosig/cloud/'


  for qq=0, nmodes-1 do begin
     header0[1] = ' Model Atmosphere Aerosol Distribution: Mode ' +strtrim(modenumber[qq],1)

     filename = filepath+fname+'__MODE_'+strtrim(modenumber[qq],1)+'.cld'
     
     thismode = tauarray_[*,qq]

                                ; nonzeros = where(thismode NE 0) ;if all zeros
                                ; IF nonzeros[0] LT 0 THEN BEGIN ;if it's ALL zeros then create arrays anyway with just the first 50 entries
                                ;    usefulhts = heightarray[0:50]
                                ;    usefulmode = thismode[0:50]
                                ; ENDIF 

     usefulmode = [0.0, thismode]
     usefulhts = [0.0, heightarray]

     openw, lun, filename, /get_lun
     for j =0, n_elements(header0) -1 do begin 
        size = strlen(header0[j])
        printf, lun, header0[j], format='(a'+strtrim(size,1)+')'
     endfor

     free_lun, lun


     openw, 1, filename, /append
     printf, 1, header, format='(a'+strtrim(strlen(header),1)+')'
     close,1

     openw, 1, filename, /append
     for j=0, n_elements(usefulhts)-1 do printf, 1, usefulhts[j], usefulmode[j], format='(F4.1, 5X, E8.2)'
     close,1

     openw, 1, filename, /append
     printf, 1, footer, format='(a6)'
     close,1

  endfor                        ;qq (modes)

END
