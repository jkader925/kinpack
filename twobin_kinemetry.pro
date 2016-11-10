PRO TwoBin_Kinemetry, data, reff, vsys, diagnostics, name, inner, outer, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
    nopne=nopne, nogcs=nogcs, vs_fit, dvs_fit, sauron, smeagol, N0821, N1023, N1344, N2768, N3115, N3377, N4473, N4564, N4697, inplotdat, inmaxr, bin1, l, nbootstraps, twobin=twobin,$
    obser, peakplot=peakplot, qin, dqin, nbootstrapsq

diagnostics = 1


;check for vs_fitting
;check structures to hold kin results

if diagnostics GT 0 then window,0,xsize=500,ysize=500
if diagnostics GT 0 then plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3

;---------photometry results------------------------------------------------------------
if name eq 'NGC821' then begin 
  q = 0.61 
  sige = 183 
  vsys = vsys 
  ppa = 32.2d
  N0821.kqi = q
  N0821.vsysi = vsys
  N0821.sige = sige
  N0821.reff = reff
endif
if name eq 'NGC1023' then begin 
  q = 0.41 
  sige = 180 
  vsys = vsys 
  ppa = 87d
  N1023.kqi = q
  N1023.vsysi = vsys
  N1023.sige = sige
  N1023.reff = reff
endif
if name eq 'NGC1344' then begin 
  q = 0.62d
  sige = 170
  vsys = vsys 
  ppa = 165.4d
  N1344.kqi = q
  N1344.vsysi = vsys
  N1344.sige = sige
  N1344.reff = reff
endif
if name eq 'NGC2768' then begin 
  q = 0.53 
  sige = 202 
  vsys = vsys 
  ppa = 89.9d
  N2768.kqi = q
  N2768.vsysi = vsys
  N2768.sige = sige
  N2768.reff = reff
endif
if name eq 'NGC3115' then begin 
  q = 0.45 
  sige = 200 
  vsys = vsys 
  ppa = 43.5d
  N3115.kqi = q
  N3115.vsysi = vsys
  N3115.sige = sige
  N3115.reff = reff
endif
if name eq 'NGC3377' then begin 
  q = 0.50 
  sige = 132 
  vsys = vsys 
  ppa = 46.3d
  N3377.kqi = q
  N3377.vsysi = vsys
  N3377.sige = sige
  N3377.reff = reff
endif
if name eq 'NGC4473' then begin 
  q = 0.59d
  sige = 189d ;v_disp within **1KPC**, from Brodie et al. 2014
  vsys = vsys 
  ppa = 92.2d
  N4564.kqi = q
  N4564.vsysi = vsys
  N4564.sige = sige
  N4564.reff = reff
endif
if name eq 'ngc4564' then begin 
  q = 0.52 
  sige = 151 
  vsys = vsys 
  ppa = 48.5d
  N4564.kqi = q
  N4564.vsysi = vsys
  N4564.sige = sige
  N4564.reff = reff
endif
if name eq 'NGC4697' then begin 
  q = 0.64 
  sige = 180 
  vsys = vsys 
  ppa = 67.2d
  N4697.kqi = q
  N4697.vsysi = vsys
  N4697.sige = sige
  N4697.reff = reff
endif
;---------------------------------------------------------------------------------------
;print,strcompress('ppa ='+string(ppa)+'degrees')

;ppa=0d


;--------set inner annulus parameters---------------------------------------------------
if inner eq 1 then begin
  kaxisratio =  q ; Photometric Q
  kpa = 0d;ppa*!dtor ; Position angle of elliptical rolling bins
  dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort)
ind1 = [where(data.er GT 0.55*reff and data.er LE 1.2*reff)]
if name eq 'ngc4564' then ind1 = [where(data.er GT 0.85*reff and data.er LE 2.6*reff)]
endif
;---------------------------------------------------------------------------------------

;--------set outer annulus parameters---------------------------------------------------
if outer eq 1 then begin
  Nouter = 40 ;minimum number of data points to be included in outer bin
  kaxisratio =  1 ; Circular bin
  kpa = 0d ; Doesn't matter (it's a circle)
  dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort)

  k = 0
  while 1 do begin
;    ind4 = [where(data.er GT (3.5*reff - (k/10.)*reff) AND data.er LE (4.5*reff + (k/10.)*reff))]
    ind4 = [where(data.er GT (3d*reff) AND data.er LE (4.5*reff + (k/10.)*reff))]
    if n_elements(data[ind4]) GE Nouter then break
    k=k+1
    if 4.5d + (k/10.) GT 20d then break
  endwhile
;    ind4 = [where(data.er GT 3d*reff AND data.er LE 11d*reff)]
  undefine,k
  print,strcompress('N_data in outer annulus ='+string(n_elements(data[ind4])))
  print,strcompress('Outer annulus limits (R_eff) ='+string(minmax(data[ind4].er/reff)))
endif
;---------------------------------------------------------------------------------------


;-----read data into plot structures----------------------------------------------------
if name eq 'NGC821' then begin 
  if inner eq 1 then N0821.vradi = data[ind1].vel
  if inner eq 1 then N0821.dvradi = data[ind1].errvel
  if inner eq 1 then N0821.pai = data[ind1].pa
  if outer eq 1 then N0821.vrado = data[ind4].vel
  if outer eq 1 then N0821.dvrado = data[ind4].errvel
  if outer eq 1 then N0821.pao = data[ind4].pa
endif
if name eq 'NGC1023' then begin 
  if inner eq 1 then N1023.vradi = data[ind1].vel
  if inner eq 1 then N1023.dvradi = data[ind1].errvel
  if inner eq 1 then N1023.pai = data[ind1].pa
  if outer eq 1 then N1023.vrado = data[ind4].vel
  if outer eq 1 then N1023.dvrado = data[ind4].errvel
  if outer eq 1 then N1023.pao = data[ind4].pa
endif
if name eq 'NGC1344' then begin 
  if inner eq 1 then N1344.vradi = data[ind1].vel
  if inner eq 1 then N1344.dvradi = data[ind1].errvel
  if inner eq 1 then N1344.pai = data[ind1].pa
  if outer eq 1 then N1344.vrado = data[ind4].vel
  if outer eq 1 then N1344.dvrado = data[ind4].errvel
  if outer eq 1 then N1344.pao = data[ind4].pa
endif
if name eq 'NGC2768' then begin 
  if inner eq 1 then N2768.vradi = data[ind1].vel
  if inner eq 1 then N2768.dvradi = data[ind1].errvel
  if inner eq 1 then N2768.pai = data[ind1].pa
  if outer eq 1 then N2768.vrado = data[ind4].vel
  if outer eq 1 then N2768.dvrado = data[ind4].errvel
  if outer eq 1 then N2768.pao = data[ind4].pa
endif
if name eq 'NGC3115' then begin 
  if inner eq 1 then N3115.vradi = data[ind1].vel
  if inner eq 1 then N3115.dvradi = data[ind1].errvel
  if inner eq 1 then N3115.pai = data[ind1].pa
  if outer eq 1 then N3115.vrado = data[ind4].vel
  if outer eq 1 then N3115.dvrado = data[ind4].errvel
  if outer eq 1 then N3115.pao = data[ind4].pa
endif
if name eq 'NGC3377' then begin 
  if inner eq 1 then N3377.vradi = data[ind1].vel
  if inner eq 1 then N3377.dvradi = data[ind1].errvel
  if inner eq 1 then N3377.pai = data[ind1].pa
  if outer eq 1 then N3377.vrado = data[ind4].vel
  if outer eq 1 then N3377.dvrado = data[ind4].errvel
  if outer eq 1 then N3377.pao = data[ind4].pa
endif
if name eq 'NGC4473' then begin 
  if inner eq 1 then N4473.vradi = data[ind1].vel
  if inner eq 1 then N4473.dvradi = data[ind1].errvel
  if inner eq 1 then N4473.pai = data[ind1].pa
  if outer eq 1 then N4473.vrado = data[ind4].vel
  if outer eq 1 then N4473.dvrado = data[ind4].errvel
  if outer eq 1 then N4473.pao = data[ind4].pa
endif
if name eq 'ngc4564' then begin 
  if inner eq 1 then N4564.vradi = data[ind1].vel
  if inner eq 1 then N4564.dvradi = data[ind1].errvel
  if inner eq 1 then N4564.pai = data[ind1].pa
  if outer eq 1 then N4564.vrado = data[ind4].vel
  if outer eq 1 then N4564.dvrado = data[ind4].errvel
  if outer eq 1 then N4564.pao = data[ind4].pa
endif
if name eq 'NGC4697' then begin 
  if inner eq 1 then N4697.vradi = data[ind1].vel
  if inner eq 1 then N4697.dvradi = data[ind1].errvel
  if inner eq 1 then N4697.pai = data[ind1].pa
  if outer eq 1 then N4697.vrado = data[ind4].vel
  if outer eq 1 then N4697.dvrado = data[ind4].errvel
  if outer eq 1 then N4697.pao = data[ind4].pa
endif
;---------------------------------------------------------------------------------------


;-------bootstrap structures------------------------------------------------------------
nbootstraps = nbootstraps
bootstrapvs = replicate({vsys:0d},1)  ;v_systemic fitting
bootstrap1 = replicate({innerv:0d, innerpa:0d},1)
bootstrap2 = replicate({outerv:0d, outerpa:0d, vdisp:0d},1)
bootstrapq = replicate({qin:0d},1)
;---------------------------------------------------------------------------------------






;-------kinemetry to recover v_sys---------------------------------------------------------
IF inner EQ 1 AND vs_fitting EQ 1 THEN BEGIN
FOR i=0,nbootstraps-1 DO BEGIN
  
    ;v_sys fitting in bins appropriate for respective data sets
    sauron = 1
    smeagol = 0
    
    IF sauron EQ 1 THEN BEGIN
    sau_vsysfit = [where(data.er GT 0.5*reff and data.er LT 1.5*reff)]
    ind1 = sau_vsysfit
    randind = min(ind1) + ROUND((max(ind1)-min(ind1))*RANDOMU(seed,(max(ind1)-min(ind1))))
    data2=data[randind]
    data2=data2[sort(data2.er)]
    indresamp=[where(data2.er GT 0.5*reff and data.er LE 1.5*reff)]
    ENDIF
    
    IF smeagol EQ 1 THEN BEGIN
    sme_vsysfit = [where(data.er GT 0.5*reff and data.er LT 2.0*reff)]
    ind1 = sme_vsysfit
    randind = min(ind1) + ROUND((max(ind1)-min(ind1))*RANDOMU(seed,(max(ind1)-min(ind1))))
    data2=data[randind]
    data2=data2[sort(data2.er)]
    indresamp=[where(data2.er GT 0.5*reff and data.er LE 2.0*reff)]
    ENDIF
    
    ;for testing iterations on non-resampled data
;    data2=data
;    indresamp=ind1
    
    
    if diagnostics eq 1 then plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3    
    if diagnostics eq 1 then oplot,data2.x,data2.y,psym=3,color=fsc_color('powder blue')

    vsysfit=1
    vfit=0
    qfit=0
    pafit=0
    obs=1
    inner=1
    outer=0
    fullprofile=0
    
    kaxisratio = q
    
    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile, vsysfit
    result_vs = kinemetry(data2.pa, data2.vel, data2.errvel, pi, start, sige, quiet=1, diag=0)
    bootstrapvs.vsys = result_vs.sys

    IF i eq 0 THEN BEGIN 
      bootstrap_vs = bootstrapvs
    ENDIF ELSE BEGIN
      bootstrap_vs = struct_append(bootstrap_vs,bootstrapvs)
    ENDELSE

undefine,vel
ENDFOR

vs_fit = mean(bootstrap_vs.vsys)
dvs_fit = stddev(bootstrap_vs.vsys)
ENDIF
;---------------------------------------------------------------------------------------






;-------inner annulus kinemetry---------------------------------------------------------
IF inner EQ 1 THEN BEGIN
  IF vs_fitting eq 1 then vsys = vs_fit
  FOR i=0,nbootstraps-1 DO BEGIN
  
    ;regular kinemetry bin
    if i eq 0 then begin 
    indresamp = ind1
    data2 = data
    if name eq 'NGC4697' then data2 = data2[where(data2.errvel GT 0d)]
    endif
    if i gt 0 then begin
    randind = min(ind1) + ROUND((max(ind1)-min(ind1))*RANDOMU(seed,(max(ind1)-min(ind1))))
    data2=data[randind]
    data2=data2[sort(data2.er)]
    if name eq 'NGC4697' then data2 = data2[where(data2.errvel GT 0d)]
    indresamp=[where(data2.er GT 0.55*reff and data.er LE 1.2*reff)]
    if name eq 'ngc4564' then indresamp=[where(data2.er GT 0.15*reff and data.er LE 2.6*reff)]
    endif
    
    
    ;outlier rejection (2-sigma)
;    samp=data2[indresamp]
;    samp=samp[sort(samp.pa)]
;    
;    binn = replicate({VEL:0d,ERRVEL:0d,X:0d,Y:0d,R:0d,PA:0d,A:0d,B:0d,ER:0d},1)
;    holdbins = binn
;    p=plot(data2[indresamp].pa,data2[indresamp].vel,linestyle=6,symbol='o',sym_filled=1)
;    
;    for j=0,2 do begin
;      bin = samp[where(samp.pa GT j*1.5d and samp.pa LE (j+1)*1.5d)]
;      kk=0
;      k=0
;      for k = 0,n_elements(bin)-1 do begin
;        if abs(bin[k].vel) LT 2d*stddev(bin.vel) then begin
;          binn[kk] = bin[k]
;          kk=kk+1
;        endif
;      endfor
;      holdbins = struct_append(holdbins,binn)
;    endfor
;    data2[indresamp]=holdbins
;    stop
    
    
    ;for testing iterations on non-resampled data
;    data2=data
;    indresamp=ind1
    
    if diagnostics GT 0 then plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3    
    if diagnostics eq 1 then  oplot,data2[indresamp].x,data2[indresamp].y,psym=3,color=fsc_color('red')


  ;---------------------------------v_rot fit
    vsysfit=0
    vfit=0
    pafit=1
    obs=1
    qfit=0
    fullprofile=0

    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile, vsysfit
    result_v = kinemetry(data2[indresamp].pa, data2[indresamp].vel, data2[indresamp].errvel, pi, start, sige, quiet=1, diag=0)
    bootstrap1.innerv = result_v.quan/sige
    bootstrap1.innerpa = result_v.PA/!dtor


   if peakplot eq 1 and i eq 0 then begin ;this is inside a loop already, only do this loop on first iteration!
    ;---------------------------------Q_kin fit
    FOR kk = 0,nbootstrapsq-1 DO BEGIN
      kaxisratio =  q
      kpa = 0d;ppa*!dtor ; Position angle of elliptical rolling bins
      dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort)
      indhalf = [where(data.er GT 0.4*reff and data.er LE 0.6*reff)]
    
      ;regular kinemetry bin
      if kk eq 0 then begin 
        indresampq = indhalf
        data2 = data
        if name eq 'NGC4697' then data2 = data2[where(data2.errvel GT 0d)]
      endif
      if kk gt 0 then begin
        randindq = min(indhalf) + ROUND((max(indhalf)-min(indhalf))*RANDOMU(seed,(max(indhalf)-min(indhalf))))
        data2=data[randindq]
        data2=data2[sort(data2.er)]
        if name eq 'NGC4697' then data2 = data2[where(data2.errvel GT 0d)]
        indresampq=[where(data2.er GT 0.3*reff and data.er LE 0.7*reff)]
        if name eq 'ngc4564' then indresampq=[where(data2.er GT 0.3*reff and data.er LE 0.7*reff)]
      endif
    
    
      vsysfit=0
      vfit=0
      pafit=0
      qfit=1
      obs=1
      fullprofile=1
    
      IF kk EQ 0 THEN bootstrapq1 = bootstrapq
    
      window,1,xsize=500,ysize=500
      plot,data.x/reff,data.y/reff,psym=3,xrange=[-4,4],yrange=[-4,4]
      oplot,data2.x/reff,data2.y/reff,psym=3,color=fsc_color('orange')

      rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile, vsysfit
      result_q = kinemetry(data2.pa, data2.vel, data2.errvel, pi, start, sige, quiet=1, diag=0)
      bootstrapq1.qin = result_q.axisratio
    

      IF kk eq 0 THEN BEGIN
        bootstrapq = bootstrapq1
      ENDIF ELSE BEGIN
        bootstrapq = struct_append(bootstrapq,bootstrapq1)
      ENDELSE
    ENDFOR  
    qin = median(bootstrapq[where(bootstrapq.qin LT 1d)].qin)
    dqin = stddev(bootstrapq[where(bootstrapq.qin LT 1d)].qin)
    print,strcompress('inner Q_kin ='+string(qin)+'+/-'+string(dqin))
    
if name eq 'NGC821' then begin 
  N0821.kqi = qin
endif
if name eq 'NGC1023' then begin 
  N1023.kqi = qin
endif
if name eq 'NGC1344' then begin 
  N1344.kqi = qin
endif
if name eq 'NGC2768' then begin 
  N2768.kqi = qin
endif
if name eq 'NGC3115' then begin 
  N3115.kqi = qin
endif
if name eq 'NGC3377' then begin 
  N3377.kqi = qin
endif
if name eq 'NGC4473' then begin 
  N4564.kqi = qin
endif
if name eq 'ngc4564' then begin 
  N4564.kqi = qin
endif
if name eq 'NGC4697' then begin 
  N4697.kqi = qin
endif


ENDIF 

  
   ;------------------------------------------
   


    IF i eq 0 THEN BEGIN 
      bootstrap_in = bootstrap1
    ENDIF ELSE BEGIN
      bootstrap_in = struct_append(bootstrap_in,bootstrap1)
    ENDELSE

undefine,vel
ENDFOR

;careful -- need to disregard answers with very small vrot result!
uniquepa = bootstrap_in[uniq(bootstrap_in.innerpa)].innerpa
uniquev = bootstrap_in[uniq(bootstrap_in.innerv)].innerv
innerv = mean(uniquev[where(uniquev GT 0.001)])
;innerv = mean(abs(bootstrap_in[uniq(bootstrap_in.innerv)].innerv))
;dinnerv = stddev(abs(bootstrap_in[uniq(bootstrap_in.innerv)].innerv))
dinnerv = stddev(uniquev[where(uniquev GT 0.001)])
innerpa = mean(uniquepa[where(uniquepa GT -180d AND uniquepa NE 0d AND uniquepa LT 180)])
dinnerpa = stddev(uniquepa[where(uniquepa GT -180d AND uniquepa NE 0d AND uniquepa LT 180)])


print,strcompress('v_rot_inner ='+string(innerv)+'+/-'+string(dinnerv))
print,strcompress('kPA_inner ='+string(innerpa)+'+/-'+string(dinnerpa))


;---------------------------------------------------------------------------------------

IF twobin EQ 0 THEN BEGIN

  mininc=5
  inc=10
  incr=1

  bootstrap = replicate({er:0d,v:0d},1)


  j=0
  kpa = 0d;innerpa*!dtor
  kaxisratio = q
  
  dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort) ;Install new a and b values for each object. Also sort by elliptical radius.
  ind = lindgen(n_elements(data))             ;indices of the data
  rrange=15d
  if name eq 'NGC4697' then data = data[where(data.errvel GT 0d)]
  if diagnostics GT 0 then plot,data.x,data.y,xrange=[-4,4]*reff,yrange=[-4,4]*reff,psym=3    
  
  
  
  WHILE 1 DO BEGIN
    ;manually specified initial guesses
    kpa = innerpa*!dtor
    vsys = 0d
    kaxisratio = q
    fullprofile = 1
    inner = 1
    outer = 0
    qfit = 0
    vfit=1
    pafit=0
    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile, vsysfit
    
      if ((incr*j)) GE n_elements(ind) then break
      if ((incr*j+inc)) GE n_elements(ind) then break
      use = ind[incr*j:min([incr*j+inc,n_elements(data)])-1]            ;isolate the data for this elliptical bin
      if (n_elements(use) lt mininc) then break
      if max(data[use].er) GT rrange*reff then break
      ;.......................plot elliptical bins
         if diagnostics eq 1 then begin
      if j ne 0 then oplot,data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].x, $
        data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].y,psym=3
      oplot,data[use].x,data[use].y,psym=3,color=fsc_color('lime green')
      endif
      ;.......................
      for i = 0,nbootstraps-1 do begin
        if i eq 0 then randind = use
      if i gt 0 then begin
        randind = min(use) + ROUND((max(use)-min(use))*RANDOMU(seed,(max(use)-min(use))))
        data2=data[randind]
        data2=data2[sort(data2.er)]
      endif
      oplot,data[randind].x,data[randind].y,psym=3,color=fsc_color('purple')
      er = median(data[randind].er/reff)
      hold1 = kinemetry(data[randind].pa, data[randind].vel, data[randind].errvel, pi, start, sige, quiet=1, diag=0)
;      if size(hold1,/N_DIMENSIONS) EQ 0 then break
      hold1 = struct_addtags(hold1, {er:mean(data[use].er)})
      bootstrap.er = hold1.er
      bootstrap.v = hold1.quan
      IF i eq 0 THEN BEGIN 
        bootstrap_in = bootstrap
      ENDIF ELSE BEGIN
        bootstrap_in = struct_append(bootstrap_in,bootstrap)
      ENDELSE
      bin1 = struct_append(bin1, hold1) ; holds outer bin fitted vrot
      endfor

    j++
  ENDWHILE
  
        binfullprof = bin1
        fullprofile = 1
        twobin = 0
        stellar = 1 ;just for plot title
        discrete = 0
        cosmo = 0
        peakplot = 1
        
        
;        p=plot(bin1.er/reff,bin1.quan/sige,linestyle=6,symbol='dot',xrange=[0,6],yrange=[-0.3,2.5])
;        Kinemetry_profile, reff, sige, q, data, rrange, bin44, bin55, bin1, bin2, twobin=twobin, binfullprof, fullprofile=fullprofile, name, stellar, discrete, inplotdat, inmaxr, outx, outy
         Kinemetry_profile_adbin, reff, sige, q, data, rrange, bin44, bin55, bin1, bin2, twobin=twobin, binfullprof, fullprofile=fullprofile, name, $
          stellar, discrete, inplotdat, inmaxr, outx, outy, obs, cosmo, l, peakvQ, peakplot=peakplot, maxr
 ENDIF
;---------------------------------------------------------------------------------------
ENDIF
;---------------------------------------------------------------------------------------



;-------outer annulus kinemetry---------------------------------------------------------
IF outer EQ 1 AND sauron EQ 0 and smeagol EQ 0 THEN BEGIN
FOR i=0,nbootstraps-1 DO BEGIN
  
  
    if i eq 0 then begin 
    indresamp = ind4
    data2 = data
    endif
    if i gt 0 then begin
    randind = min(ind4) + ROUND((max(ind4)-min(ind4))*RANDOMU(seed,(max(ind4)-min(ind4))))
    data2=data[randind]
    data2=data2[sort(data2.er)]
    indresamp=[where(data2.er GT min(data[ind4].er) and data2.er LE max(data[ind4].er))]
    endif
    
;    indresamp=[where(data2.er GT 3.0*reff and data2.er LE 4.9*reff)]
;     indresamp=[where(data2.er GT 2.5*reff and data2.er LE 5.5*reff)]
;     indresamp=[where(data2.er GT 3.5*reff and data2.er LE 6.5*reff)]
     
;     for testing iterations on non-resampled data
;    data2=data
;    indresamp=ind4
    
    if diagnostics GT 0 then plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3    
    if diagnostics eq 1 then oplot,data2[indresamp].x,data2[indresamp].y,psym=3,color=fsc_color('red')    
    
;    vfit=1
;    pafit=0
;    obs=1
;;    kaxisratio = 1
;    kpa = ppa*!dtor
;; initial vel guess from previous iteration
;    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, obs, inner, outer
;    result_v = kinemetry_max(data2[indresamp].pa, data2[indresamp].vel, data2[indresamp].errvel, pi, start, lam=lam, verbose=verbose)

    vfit=0
    pafit=1
    obs=1
    fullprofile=0
    qfit = 0
    outer = 1
    inner =0

    kpa = 30*!dtor
    kaxisratio = 1
    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile, vsysfit
    window,0,xsize=500,ysize=500
    k=0
    while 1 do begin
    result_p = kinemetry_max(data2[indresamp].pa, data2[indresamp].vel, data2[indresamp].errvel, pi, start, lam=lam, verbose=verbose, sige, reff)
    if size(result_p,/N_DIMENSIONS) EQ 1 then break
    if k GT 1000 then break
    k=k+1
    endwhile
    bootstrap2.outerv = result_p.quan/sige
    bootstrap2.outerpa = result_p.PA/!dtor
    bootstrap2.vdisp = result_p.disp/sige


    IF i eq 0 THEN BEGIN 
      bootstrap_out = bootstrap2
    ENDIF ELSE BEGIN
      bootstrap_out = struct_append(bootstrap_out,bootstrap2)
    ENDELSE

undefine,result_v
;undefine,vel
;undefine,kpa
ENDFOR

;unique = bootstrap_out[uniq(bootstrap_out.outerpa)].outerpa
unique = bootstrap_out[where(bootstrap_out.outerpa GT -180 AND bootstrap_out.outerpa NE 0d AND bootstrap_out.outerpa LT 180)].outerpa
outerv = mean(abs(bootstrap_out[uniq(bootstrap_out.outerv)].outerv))
douterv = stddev(abs(bootstrap_out[uniq(bootstrap_out.outerv)].outerv))
;outerpa = mean(unique[where(unique NE 0d AND unique LT 180)])
outerpa = mean(unique)
;douterpa = stddev(bootstrap_out[uniq(bootstrap_out.outerpa)].outerpa)
douterpa = stddev(unique)
vdisp = mean(bootstrap_out.vdisp)
dvdisp = stddev(bootstrap_out.vdisp)
;N0821.vdisp = vdisp
;N0821.dvdisp = dvdisp

print,strcompress('v_rot_outer ='+string(outerv)+'+/-'+string(douterv))
print,strcompress('kPA_outer ='+string(outerpa)+'+/-'+string(douterpa))


;---------------------------------------------------------------------------------------
IF twobin EQ 0 THEN BEGIN
  bootstrap = replicate({er:0d,v:0d},1)

  mininc=5
  inc=10
  incr=1


  kpa = 0d
  j=0
  dummy = ellippar(data.x, data.y, q, kpa, struc=data, /sort) ;Install new a and b values for each object. Also sort by elliptical radius.
  ind = lindgen(n_elements(data))             ;indices of the data
  rrange=15d
  
  WHILE 1 DO BEGIN
    ;manually specified initial guesses
    kpa = ppa*!dtor;outerpa*!dtor;0d;30d*!dtor
    vsys = 0d
    kaxisratio = q
    fullprofile=1
    obs=1
    outer = 1

    vfit=1
    pafit=0
    rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile, vsysfit
    
      if ((incr*j)) GE n_elements(ind) then break
      if ((incr*j+inc)) GE n_elements(ind) then break
      use = ind[incr*j:min([incr*j+inc,n_elements(data)])-1]            ;isolate the data for this elliptical bin
      if (n_elements(use) lt mininc) then break
      if max(data[use].er) GT rrange*reff then break
      ;.......................plot elliptical bins
         if diagnostics eq 1 then begin
      if j ne 0 then oplot,data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].x, $
        data[ind[incr*(j-1):min([incr*(j-1)+inc,n_elements(data)])-1]].y,psym=3
      oplot,data[use].x,data[use].y,psym=3,color=fsc_color('lime green')
      endif
      ;.......................
        for i = 0,nbootstraps-1 do begin
        if i eq 0 then randind = use
      if i gt 0 then begin
        randind = min(use) + ROUND((max(use)-min(use))*RANDOMU(seed,(max(use)-min(use))))
        data2=data[randind]
        data2=data2[sort(data2.er)]
      endif
      oplot,data[randind].x,data[randind].y,psym=3,color=fsc_color('purple')
    hold1 = kinemetry_max(data[randind].pa, data[randind].vel, data[randind].errvel, pi, start, lam=lam, verbose=verbose, sige, reff)
      if size(hold1,/N_DIMENSIONS) EQ 0 then break
      hold1 = struct_addtags(hold1, {er:mean(data[use].er)})
;      bin2 = struct_append(bin2, hold1) ; holds outer bin fitted vrot
          bootstrap.er = hold1.er
      bootstrap.v = hold1.quan
      IF i eq 0 THEN BEGIN 
        bootstrap_in = bootstrap
      ENDIF ELSE BEGIN
        bootstrap_in = struct_append(bootstrap_in,bootstrap)
      ENDELSE
      bin2 = struct_append(bin2, hold1) ; holds outer bin fitted vrot
      endfor
    
    j++
  ENDWHILE
        binfullprof = bin2
        fullprofile = 1
        twobin = 0

        stellar = 0 ;just for plot title
        discrete = 1
        cosmo = 0
        obs = 1
        peakplot = 1
        
        
        
;        p=plot(bin1.er/reff,bin1.quan/sige,linestyle=6,symbol='dot',xrange=[0,8],yrange=[0,2])
;        p=plot(bin1.er/reff,bin1.pa/!dtor,linestyle=6,symbol='dot',xrange=[0,8],yrange=[-180,180])
;        Kinemetry_profile, reff, sige, q, data, rrange, bin44, bin55, bin1, bin2, twobin=twobin, binfullprof, fullprofile=fullprofile, name, stellar, discrete, inplotdat, inmaxr, outx, outy

        Kinemetry_profile_adbin, reff, sige, q, data, rrange, bin44, bin55, bin1, bin2, twobin=twobin, binfullprof, fullprofile=fullprofile, name, $
          stellar, discrete, inplotdat, inmaxr, outx, outy, obs, cosmo, l, peakvQ, peakplot=peakplot, maxr, bin

        


        IF name EQ 'NGC821' THEN BEGIN
          N0821 = struct_addtags(N0821, {rawer:bin.er/reff, rawv:bin.quan/sige})
          N0821.plotx = outx
          N0821.ploty = outy
          obser = N0821
        ENDIF
        IF name EQ 'NGC1023' THEN BEGIN
          N1023 = struct_addtags(N1023, {rawer:bin.er/reff, rawv:bin.quan/sige})        
          N1023.plotx = outx
          N1023.ploty = outy
          obser = N1023
        ENDIF
        IF name EQ 'NGC2768' THEN BEGIN
          N2768 = struct_addtags(N2768, {rawer:bin.er/reff, rawv:bin.quan/sige})        
          N2768.plotx = outx
          N2768.ploty = outy
          obser = N2768
        ENDIF
        IF name EQ 'NGC3115' THEN BEGIN
          N3115 = struct_addtags(N3115, {rawer:bin.er/reff, rawv:bin.quan/sige})        
          N3115.plotx = outx
          N3115.ploty = outy
          obser = N3115
        ENDIF
        IF name EQ 'NGC3377' THEN BEGIN
          N3377 = struct_addtags(N3377, {rawer:bin.er/reff, rawv:bin.quan/sige})        
          N3377.plotx = outx
          N3377.ploty = outy
          obser = N3377
        ENDIF
        IF name EQ 'ngc4564' THEN BEGIN
          N4564 = struct_addtags(N4564, {rawer:bin.er/reff, rawv:bin.quan/sige})        
          N4564.plotx = outx
          N4564.ploty = outy
          obser = N4564
        ENDIF
        IF name EQ 'NGC4697' THEN BEGIN
          N4697 = struct_addtags(N4697, {rawer:bin.er/reff, rawv:bin.quan/sige})        
          N4697.plotx = outx
          N4697.ploty = outy
          obser = N4697
        ENDIF
        
        undefine,bin
        undefine,bin1
        undefine,bin2
        
  ENDIF      
ENDIF
;---------------------------------------------------------------------------------------







END


