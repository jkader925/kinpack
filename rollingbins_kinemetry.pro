PRO RollingBins_Kinemetry, data, q, bin, twobins=twobins, fullprofile=fullprofile, diagnostics=diagnostics, cosmo=cosmo, reff, sige, $
                            innerbinwidth,innerbinQ,innerbinPA,Nbootstraps_in,Nbootstraps_out,innerbinINC,innerbinINCR,profplot,F1,obs=obs, $
                            mininc, inc, incr, rrange, outx, outy, l, peakvQ, peakplot=peakplot
      
      
;Diagnostic Plots
if diagnostics GT 0 then window,0,xsize=500,ysize=500
if diagnostics GT 0 then plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3
data.vel = -data.vel
;...................................................

;Basic Bin Parameters
;rrange=30d  ;Elliptical radius (R/R_e) of furthest-out kinemetry bin
;mininc = 150d ; Minimum Nstars per rolling bin
;inc = 250d < n_elements(data); (Nominal) rolling bin size
;incr = 250d ; Nstars to increment by


;if cosmo eq 1 then begin
;  mininc=500d
;  inc=1000d
;  incr=10000d
;endif

if obs eq 1 then begin
  mininc=5
  inc=10
  incr=1
endif


;Bin Geometry
kaxisratio = q ; Shape of rolling bins
kpa = 0.*!dtor ; Position angle of semi-major axis (for eliptical bins)

;Generate elliptical coordinates
dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort) ;Install new a and b values for each object. Also sort by elliptical radius.
ind = lindgen(n_elements(data))             ;indices of the data

;...................................................



IF fullprofile EQ 1 THEN BEGIN


  j=0
  WHILE 1 DO BEGIN
    ;manually specified initial parameter guesses
    kpa = 0d*!dtor;30d*!dtor
    vsys = 0d
    kaxisratio = q
    
    outer=1
    inner=0
    vfit=1
    pafit=0
    qfit=0
    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile
    
;...................................................
      if ((incr*j)) GE n_elements(ind) then break
      if ((incr*j+inc)) GE n_elements(ind) then break
      use = ind[incr*j:min([incr*j+inc,n_elements(data)])-1]            ;isolate the data for this elliptical bin
      if (n_elements(use) lt mininc) then break
      if max(data[use].er) GT rrange*reff then break
      ;.......................plot elliptical bins
         if diagnostics eq 1 then begin
      if j ne 0 then oplot,data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].x, $
        data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].y,psym=3
      oplot,data[use].x,data[use].y,psym=3,color=fsc_color('purple')
      endif
      ;.......................
      hold = kinemetry_max(data[use].pa, data[use].vel, data[use].errvel, pi, start, lam=lam, verbose=verbose)
;      if size(hold,/N_DIMENSIONS) EQ 0 then break
      hold = struct_addtags(hold, {er:mean(data[use].er)})
      bin = struct_append(bin, hold) ; holds outer bin fitted vrot
;...................................................

    j++
  ENDWHILE
  
        binfullprof = bin
        stellar = 0 ;just for plot title
        discrete = 1
        twobin = 0
        bin1=bin
  
        Kinemetry_profile_adbin_sims, reff, sige, q, data, rrange, bin44, bin55, bin1, bin2, twobin=twobin, binfullprof, fullprofile=fullprofile, name, stellar, discrete, inplotdat, inmaxr, outx, $
              outy, obs, cosmo, l, peakvQ, peakplot=peakplot, maxr, bin
          
          
IF peakplot EQ 1 THEN BEGIN
  plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3
          j=0
      ind1 = [where(abs(reff/2d - data.er) LT 0.1d)] ;Define 1Re annulus
      if diagnostics eq 1 then  oplot,data[ind1].x,data[ind1].y,psym=3,color=fsc_color('red')
      ind=ind1
          
    kpa = 30d*!dtor
    vsys = 0d
    kaxisratio = q
    
    
    outer=1
    inner=0
    vfit=0
    pafit=0
    qfit=1
    fullprofile=1
    
    
    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile

;...................................................
      ll=0
      while 1 do begin
        oplot,data[ind1].x,data[ind1].y,psym=3,color=fsc_color('orange')
        holdq = kinemetry_max2(data[ind1].pa, data[ind1].vel, data[ind1].errvel, pi, start, lam=lam, verbose=verbose)
        holdq = struct_addtags(holdq, {er:mean(data[use].er)})
        IF holdq.axisratio LT 1d AND holdq.axisratio GE 0.1 THEN BEGIN
          binq = struct_append(binq, holdq) ; holds outer bin fitted vrot
          BREAK
        ENDIF
        ll=ll+1
        ind1 = [where(abs(reff/2d - data.er) LT 0.1d + 0.1d*ll)] ;Define 1Re annulus
      endwhile
;...................................................

     if peakplot eq 1 then begin
      if l eq 0 then begin
        openw,l+1,'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\qkin_cosmo_08_11_2016.dat'
        out=dblarr(3,1)
        out[0,*]=maxr
        out[1,*]=binq.axisratio
        out[2,*]=l
        printf,l+1,FORMAT='(3F)',out
        close,l+1
      endif


      if l gt 0 then begin
        openu,l+1,'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\qkin_cosmo_08_11_2016.dat',/append
        out=dblarr(3,1)
        out[0,*]=maxr
        out[1,*]=binq.axisratio
        out[2,*]=l
        printf,l+1,FORMAT='(3F)',out
        close,l+1
      endif
     endif
 
ENDIF
          
;        Kinemetry_profile, reff, sige, q, data, rrange, bin44, bin55, bin, binfullprof, twobin=twobin, fullprofile=fullprofile, name, stellar, discrete, inplotdat, inmaxr, outx, outy
ENDIF
;...................................................
;...................................................
;...................................................
IF twobins EQ 1 THEN BEGIN
obs=0 ;since this rollingbins code is only used on simulations anyway

      ;**********************1 R_e bin********************
      ;***************************************************

      ;...................................................
      ;Basic Bin Parameters
      inc = innerbinINC < n_elements(data); (Nominal) rolling bin size
      mininc = ROUND(inc/3d) ; Minimum Nstars per rolling bin
      incr = innerbinINCR ; Nstars to increment by

      ;Bin geometry
      innerbinwidth=innerbinwidth
      kaxisratio = innerbinQ ; Shape of annulus
      kpa = innerbinPA ; Orientation of annulus
      dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort) ;Install new a and b values for each object. Also sort by elliptical radius.
      ind1 = [where(abs(reff - data.er) LT innerbinwidth)] ;Define 1Re annulus
      if diagnostics eq 1 then  oplot,data[ind1].x,data[ind1].y,psym=3,color=fsc_color('red')
      ind=ind1
      ;...................................................
      ;...................................................
      
;          simfile = 'C:/Users/Justin/DATA_v14.fits'
;    obj = mrdfits(simfile,1,/silent)
;    info_loc = mrdfits(simfile,2,/silent)
;    dummy = ellippar(obj.x,obj.y,0.6d,0d,struc=obj,/sort)
;    an = obj[where(obj.er GT 0.9*info_loc.reff AND obj.er LT 1.1*info_loc.reff)]
;    an = struct_addtags(an, arr_struct({pa:pacalc(an.x, an.y, 0d)*!dtor}))
;    
;    raw = data[where(data.er GT 0.9*reff AND data.er LT 1.1*reff)]
;    
;    p=plot(data.x,data.y,linestyle=6,symbol='dot',xrange=[-30,30],yrange=[-30,30],title='raw data',xtitle='x',ytitle='y')
;    p=plot(raw.x,raw.y,linestyle=6,symbold='dot',color='red',/overplot)
;    
;    p=plot(obj.x,obj.y,linestyle=6,symbol='dot',xrange=[-30,30],yrange=[-30,30],title='voronoi cells',xtitle='x',ytitle='y')
;    p=plot(an.x,an.y,linestyle=6,symbol='dot',color='lime green',/overplot)
;    
;    p=plot(raw.pa/!dtor,raw.vel/sige,linestyle=6,symbol='o',sym_size=0.1,sym_filled=1,xrange=[0,360],yrange=[-1.5,1.5],title='combined rolling bin',xtitle='PA [deg]',ytitle='$v_{rad}/\sigma_e$')
;    p=plot(an.pa/!dtor,an.vel/sige,linestyle=6,symbol='o',sym_size=0.1,sym_filled=1,color='lime green',/overplot)
;      stop


      ;Bootstrapping 
      bootstrap1 = replicate({innerv:0d, outerv:0d, innerfixpa:0d, outerfixpa:0d, ivpa_mean:0d, ovpa_mean:0d, delta:0d, $
              deltafixpa:0d, deltapa_mean:0d, innerpa:0d, outerpa:0d, ipa_mean:0d, opa_mean:0d, gamma2bin:0d, gamma:0d},1)

      i=0
      FOR i=0,(Nbootstraps_in - 1) DO BEGIN
          
        randind = min(ind1) + ROUND((max(ind1)-min(ind1))*RANDOMU(seed,(max(ind1)-min(ind1))))
        data2=data[randind]
        data2=data2[sort(data2.er)]
        indresamp=[where(abs(reff - data2.er) LT innerbinwidth)]
        ind=randind[sort(randind)]
        
        
        
;        ;manually specified initial guesses
        kpa = 30d*!dtor
        vsys = 0d
        kaxisratio = q

        vfit=0
        pafit=1
        inner=1
        outer=0
        obs=0
;        rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, inner=inner, outer=outer, obs
;        while 1 do begin
;        holdin = kinemetry_max(data[ind].pa, data[ind].vel, data[ind].errvel, pi, start, lam=lam, verbose=verbose)
;        if size(holdin,/N_DIMENSIONS) EQ 1 then break
;        endwhile
;        holdin = struct_addtags(holdin, {er:mean(data[ind].er)})
;        binin = struct_append(binin, holdin) ; holds outer bin fitted vrot
        
      ;...................................................
      ;Start rolling bin algorithm in the inner annulus
      j=0
      WHILE 1 DO BEGIN
        ;manually specified initial guesses
        kpa = 10d*!dtor
        vsys = 0d
        kaxisratio = q

        vfit=0
        pafit=1
        qfit=0
        inner=1
        outer=0 ; this must be switched ON for simulation inner annulus, since we are using discrete velocity modeling
        obs=0
        rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile
;...................................................
    
        if ((incr*j)) GE n_elements(ind) then break
        if ((incr*j+inc)) GE n_elements(ind) then break
        use = ind[incr*j:min([incr*j+inc,n_elements(data)])-1]            ;isolate the data for this elliptical bin
        if (n_elements(use) lt mininc) then break
        if max(data[use].er) GT rrange*reff then break
        ;.......................plot elliptical bins
        if diagnostics eq 1 then begin
        if j ne 0 then oplot,data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].x, $
        data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].y,psym=3
        oplot,data[use].x,data[use].y,psym=3,color=fsc_color('purple')
        endif
        ;.......................
        datr = data[use].er
        holdin = kinemetry_max(data[use].pa, data[use].vel, data[use].errvel, pi, start, lam=lam, verbose=verbose, sige, datr, reff)
        if size(holdin,/N_DIMENSIONS) EQ 0 then break
        holdin = struct_addtags(holdin, {er:mean(data[use].er)})
        binin = struct_append(binin, holdin) ; holds outer bin fitted vrot

;...................................................

    j++
    ENDWHILE

    bin4=binin
    undefine,start
;     bin4 = bin4[where(bin4.pa/!dtor GT 0d)]
;     bin4 = bin4[where(bin4.pa/!dtor LT 180d)]
      bin4 = bin4[where(bin4.quan GT 0d)]
     IF n_elements(bin4.quan) EQ 1 THEN bootstrap1.innerv = bin4.quan/sige ELSE bootstrap1.innerv = mean(bin4.quan/sige)
;     IF n_elements(bin4.PA) EQ 1 THEN bootstrap1.ipa_mean = bin4.PA/!dtor ELSE bootstrap1.ipa_mean = abs(mean(bin4.pa/!dtor)) ;abs(median(bin4.PA/!dtor))
;     IF n_elements(bin4.PA) EQ 1 THEN bootstrap1.innerpa = bin4.PA/!dtor ELSE bootstrap1.innerpa = abs(median(bin4.pa/!dtor)) ;abs(median(bin4.PA/!dtor)) 
     IF n_elements(bin4.PA) EQ 1 THEN bootstrap1.ipa_mean = bin4.PA/!dtor ELSE bootstrap1.ipa_mean = mean(bin4[where(bin4.pa/!dtor GT -180d AND bin4.pa/!dtor NE 0d AND bin4.pa/!dtor LT 180d)].pa/!dtor) ;abs(median(bin4.PA/!dtor))
     IF n_elements(bin4.PA) EQ 1 THEN bootstrap1.innerpa = bin4.PA/!dtor ELSE bootstrap1.innerpa = median(bin4[where(bin4.pa/!dtor GT -180d AND bin4.pa/!dtor NE 0d AND bin4.pa/!dtor LT 180d)].pa/!dtor) ;abs(median(bin4.PA/!dtor)) 

  IF i eq 0 THEN BEGIN 
    bootstrap_in = bootstrap1
  ENDIF ELSE BEGIN
  bootstrap_in = struct_append(bootstrap_in,bootstrap1)
  ENDELSE


bin44=bin4
undefine, bin4
ENDFOR



undefine, ind











      ;**********************4 R_e bin********************
      ;***************************************************

      ;...................................................

      ;Bin geometry
      outerbinwidth=0.5*reff
      kaxisratio = 1d ; Shape of annulus
      kpa = 0.*!dtor ; Orientation of annulus
      dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort) ;Install new a and b values for each object. Also sort by elliptical radius.
      ind4 = [where(abs(4d*reff - data.er) LT outerbinwidth)] ;Define 4Re annulus
      if diagnostics eq 1 then  oplot,data[ind4].x,data[ind4].y,psym=3,color=fsc_color('red')
      ind=ind4
      ;...................................................
      
      ;...................................................
  
      ;Bootstrapping 
      bootstrap2 = replicate({innerv:0d, outerv:0d, innerfixpa:0d, outerfixpa:0d, ivpa_mean:0d, ovpa_mean:0d, delta:0d, $
              deltafixpa:0d, deltapa_mean:0d, innerpa:0d, outerpa:0d, ipa_mean:0d, opa_mean:0d, gamma2bin:0d, gamma:0d},1)

      i=0
      FOR i=0,(Nbootstraps_out - 1) DO BEGIN
;        print,i
          
        randind = min(ind4) + ROUND((max(ind4)-min(ind4))*RANDOMU(seed,(max(ind4)-min(ind4))))
        data2=data[randind]
        data2=data2[sort(data2.er)]
        indresamp=[where(abs(4d*reff - data2.er) LT outerbinwidth)]
        ind=randind[sort(randind)]
        
        
;        vfit=0
;        pafit=1
;        inner=0
;        outer=1
;        obs=0
;        rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, inner=inner, outer=outer, obs
;        while 1 do begin
;        holdout = kinemetry_max(data[ind].pa, data[ind].vel, data[ind].errvel, pi, start, lam=lam, verbose=verbose)
;        if size(holdout,/N_DIMENSIONS) EQ 1 then break
;        endwhile
;        holdout = struct_addtags(holdout, {er:mean(data[ind].er)})
;        binout = struct_append(binout, holdout) ; holds outer bin fitted vrot
        
      ;...................................................
      
      ;Start rolling bin algorithm in the inner annulus
      j=0
      WHILE 1 DO BEGIN

        ;manually specified initial guesses
        kpa = 30d*!dtor
        vsys = 0d
        kaxisratio = q
;
        vfit=0
        pafit=1
        qfit=0
        inner=0
        outer=1
        obs=0
        rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile
;...................................................
    

        if ((incr*j)) GE n_elements(ind) then break
        if ((incr*j+inc)) GE n_elements(ind) then break
        use = ind[incr*j:min([incr*j+inc,n_elements(data)])-1]            ;isolate the data for this elliptical bin
        if (n_elements(use) lt mininc) then break
        if max(data[use].er) GT rrange*reff then break
        ;.......................plot elliptical bins
        if diagnostics eq 1 then begin
        if j ne 0 then oplot,data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].x, $
        data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].y,psym=3
        oplot,data[use].x,data[use].y,psym=3,color=fsc_color('purple')
        endif
        ;.......................
        datr = data[use].er
        holdout = kinemetry_max(data[use].pa, data[use].vel, data[use].errvel, pi, start, lam=lam, verbose=verbose, sige, datr, reff)
        if size(holdout,/N_DIMENSIONS) EQ 0 then break
        holdout = struct_addtags(holdout, {er:mean(data[use].er)})
        binout = struct_append(binout, holdout) ; holds outer bin fitted vrot
;
;;...................................................
;
        undefine,start
;
;
        bin5=binout
;        ;print,strcompress('found PA:'+string(bin5.pa/!dtor))
        j++
     ENDWHILE

     if cosmo eq 0 then begin
      bin5 = bin5[where(bin5.pa/!dtor NE 0d)]
      bin5 = bin5[where(bin5.pa/!dtor LT 180d)]
     endif
     IF n_elements(bin5.quan) EQ 1 THEN bootstrap2.outerv = bin5.quan/sige ELSE bootstrap2.outerv = mean(bin5.quan/sige)
     IF n_elements(bin5.PA) EQ 1 THEN bootstrap2.opa_mean = bin5.PA/!dtor ELSE bootstrap2.opa_mean = mean(bin5[where(bin5.pa/!dtor GT -180d AND bin5.pa/!dtor NE 0d AND bin5.pa/!dtor LT 180d)].pa/!dtor) ;abs(median(bin4.PA/!dtor))
     IF n_elements(bin5.PA) EQ 1 THEN bootstrap2.outerpa = bin5.PA/!dtor ELSE bootstrap2.outerpa = mean(bin5[where(bin5.pa/!dtor GT -180d AND bin5.pa/!dtor NE 0d AND bin5.pa/!dtor LT 180d)].pa/!dtor) ;abs(median(bin4.PA/!dtor)) 

     IF i eq 0 THEN BEGIN 
     bootstrap_out = bootstrap2
     ENDIF ELSE BEGIN
     bootstrap_out = struct_append(bootstrap_out,bootstrap2)
     ENDELSE
  
     bin55 = bin5
     undefine,bin5
  ENDFOR

;...................................................
;write-out results

    innerouter = replicate({innerv:0d, outerv:0d, div:0d, dov:0d, delta:0d, derror:0d, gamma:0d, gamma2bin:0d, g2binerror:0d},1)
;    innerouter.innerv = abs(median(binin1.quan/sige))
    innerouter.innerv = abs(mean(bootstrap_in.innerv))
    innerouter.div = 2*stddev(bootstrap_in.innerv)
;    innerouter.outerv = abs(median(binout1.quan/sige))
    innerouter.outerv = abs(mean(bootstrap_out.outerv))
    innerouter.dov = 2*stddev(bootstrap_out.outerv)
    innerouter.delta = innerouter.outerv - innerouter.innerv
    innerouter.derror = sqrt(innerouter.div^2 + innerouter.dov^2)
;    innerouter.gamma = abs(median(binout.PA/!dtor))
;    innerouter.gamma2bin = abs(mean(bootstrap_out.outerpa) - mean(bininner2.PA/!dtor)) ;WAS using median for inner pa.
;    innerouter.gamma2bin = abs(mean(bootstrap_out.outerpa) - mean(bootstrap_in.innerpa)) ;WAS using median for inner pa.
    innerouter.gamma2bin = abs(mean(bootstrap_out.outerpa) - mean(bootstrap_in.innerpa)) ;WAS using median for inner pa.
    if innerouter.gamma2bin GT 180d then innerouter.gamma2bin = 360d - innerouter.gamma2bin
    innerouter.g2binerror = sqrt((2*stddev(bootstrap_out.opa_mean))^2 + (2*stddev(bootstrap_in.ipa_mean))^2)   ;stddev(bootstrap1.outerpa) ;stddev(bootstrap1.opa_mean)? We assume a NORMAL distribution so that 2*stddev is our 95% confidence interval (plotted)
    
    print,strcompress('inner rotation ='+string(innerouter.innerv)+'+/-'+string(innerouter.div))
    print,strcompress('outer rotation ='+string(innerouter.outerv)+'+/-'+string(innerouter.dov))
    print,strcompress('delta v ='+string(innerouter.delta))
    print,strcompress('delta error from bootstraps ='+string(innerouter.derror))
    print,strcompress('gamma ='+string(innerouter.gamma2bin))
    print,strcompress('gamma error from bootstraps ='+string(innerouter.g2binerror))
    F1 = innerouter
;...................................................

ENDIF ;end two-zone bins

;...................................................
;Plots
;  IF fullprofile EQ 1 THEN BEGIN
;  Kinemetry_profile, reff, sige, q, data, rrange, bin44, bin55, bin, binfullprof, twobin=twobin, fullprofile=fullprofile
;  ENDIF


END