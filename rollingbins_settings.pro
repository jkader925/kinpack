PRO rollingbins_settings, vel, kpa, kaxisratio, vsys, veldisp, veldisp_sys, start, pi, vfit=vfit, pafit=pafit, qfit=qfit, obs, inner=inner, outer=outer, fullprofile=fullprofile, vsysfit

IF inner EQ 1 THEN BEGIN
  ;generate automated initial guesses if none provided above
  if (n_elements(vel) eq 0) then start = randomu(seed)*1000d else start = vel
  if (n_elements(kpa) eq 0) then start = [start, abs(randomu(seed)*!pi - !pi/2d)] else start = [start, kpa]
  if (n_elements(kaxisratio) eq 0) then start = [start, 1d] else start = [start, kaxisratio]
  if (n_elements(vsys) eq 0) then start = [start, 0d] else start = [start, vsys]
  if (n_elements(veldisp) eq 0) then start = [start, 150d] else start = [start, veldisp]
  if (n_elements(veldisp) eq 0) then start = [start, 150d] else start = [start, veldisp]

  if vsysfit EQ 1 then begin
    undefine,vel
    undefine,kpa
    undefine,vsys
  endif  

  if qfit EQ 1 then begin
    undefine,vel
    undefine,kaxisratio
    undefine,kpa
  endif

  if vfit GT 0 then begin
;     undefine,kpa ;=30
;     undefine,kaxisratio
      undefine,vel
;     if obs eq 1 then undefine,veldisp
;      if obs eq 0 then undefine,vsys ;=0
  endif
  
  if pafit GT 0 then begin
      undefine,kpa
;     undefine,kaxisratio
     undefine,vel
     if obs eq 0 then undefine,veldisp
;     if obs eq 0 then undefine,vsys ;=0
  endif
  
  ;model parameter fixed/free
  IF fullprofile EQ 1 THEN BEGIN
  ;model parameter fixed/free
  pi = [  {fixed:n_elements(vel)      , limited:[1,1], limits:[-10000d,10000d], step:0d}, $ ;rotational velocity
        {fixed:n_elements(kpa)      , limited:[1,1], limits:[-2d*!pi/2d,2d*!pi/2d], step:0d},  $ ;kinematic position angle
        {fixed:n_elements(kaxisratio) , limited:[1,1], limits:[0.0001d,1.00001d], step:0d},   $ ;kinematic axisratio
        {fixed:n_elements(vsys)     , limited:[1,1], limits:[-1d,10000d], step:0d},  $    ;systemic velocity
        {fixed:n_elements(veldisp)  , limited:[1,1], limits:[-3000d,3000d], step:0d},  $   ;velocity dispersion
        {fixed:n_elements(veldisp_sys)  , limited:[1,1], limits:[-1d,10000d], step:0d}   ]  ;systemic vdisp
  ENDIF ELSE BEGIN
    pi = [  {fixed:n_elements(vel)      , limited:[1,1], limits:[0d,10000d], step:0d}, $ ;rotational velocity
        {fixed:n_elements(kpa)      , limited:[1,1], limits:[-2d*!pi/2d,2d*!pi/2d], step:0d},  $ ;kinematic position angle
        {fixed:n_elements(kaxisratio) , limited:[1,1], limits:[0.0001d,1.00001d], step:0d},   $ ;kinematic axisratio
        {fixed:n_elements(vsys)     , limited:[1,1], limits:[-1d,10000d], step:0d},  $    ;systemic velocity
        {fixed:n_elements(veldisp)  , limited:[1,1], limits:[-3000d,3000d], step:0d},  $   ;velocity dispersion
        {fixed:n_elements(veldisp_sys)  , limited:[1,1], limits:[-1d,10000d], step:0d}   ]  ;systemic vdisp
  ENDELSE
ENDIF

IF outer EQ 1 THEN BEGIN

  ;initial guess for kPA > 0 degrees. Provide guess here, below would possibly give kPA < 0 as initial guess!
 
  ;generate automated initial guesses if none provided above
  if (n_elements(vel) eq 0) then start = randomu(seed)*100d else start = vel
  if (n_elements(kpa) eq 0) then start = [start, abs(randomu(seed)*!pi - !pi/2d)] else start = [start, kpa] ;had to change so that random guess is >0
  if (n_elements(kaxisratio) eq 0) then start = [start, 1d] else start = [start, kaxisratio]
  if (n_elements(vsys) eq 0) then start = [start, 0d] else start = [start, vsys]
  if (n_elements(veldisp) eq 0) then start = [start, 150d] else start = [start, veldisp]
  if (n_elements(veldisp) eq 0) then start = [start, 150d] else start = [start, veldisp]


  if qfit EQ 1 then begin
    undefine,vel
    undefine,veldisp
    undefine,kaxisratio
    undefine,kpa
  endif


  if vfit EQ 1 then begin
;     undefine,kpa ;=30
;     undefine,kaxisratio
      undefine,vel
     undefine,veldisp
;     undefine,vsys
      if obs eq 0 then undefine,vsys ;=0
  endif
  
    if vfit EQ 2 then begin
;     undefine,kpa ;=30
;     undefine,kaxisratio
      undefine,vel
     undefine,veldisp
;     undefine,vsys
;      if obs eq 0 then undefine,vsys ;=0
  endif
  
  if pafit GT 0 then begin
      undefine,kpa
     ;undefine,kaxisratio
     undefine,vel
     undefine,veldisp
;     if obs eq 1 then undefine,veldisp
;     if obs eq 0 then undefine,vsys ;=0
  endif
  
  IF fullprofile EQ 1 THEN BEGIN
  ;model parameter fixed/free
  pi = [  {fixed:n_elements(vel)      , limited:[1,1], limits:[-1000d,1000d], step:0d}, $ ;rotational velocity
        {fixed:n_elements(kpa)      , limited:[1,1], limits:[-2d*!pi/2d,2d*!pi/2d], step:0d},  $ ;kinematic position angle
        {fixed:n_elements(kaxisratio) , limited:[1,1], limits:[0.0001d,1.00001d], step:0d},   $ ;kinematic axisratio
        {fixed:n_elements(vsys)     , limited:[1,1], limits:[-1d,10000d], step:0d},  $    ;systemic velocity
        {fixed:n_elements(veldisp)  , limited:[1,1], limits:[-3000d,3000d], step:0d},  $   ;velocity dispersion
        {fixed:n_elements(veldisp_sys)  , limited:[1,1], limits:[-1d,10000d], step:0d}   ]  ;systemic vdisp
  ENDIF ELSE BEGIN
    pi = [  {fixed:n_elements(vel)      , limited:[1,1], limits:[0d,10000d], step:0d}, $ ;rotational velocity
        {fixed:n_elements(kpa)      , limited:[1,1], limits:[-2d*!pi/2d,2d*!pi/2d], step:0d},  $ ;kinematic position angle
        {fixed:n_elements(kaxisratio) , limited:[1,1], limits:[0.0001d,1.00001d], step:0d},   $ ;kinematic axisratio
        {fixed:n_elements(vsys)     , limited:[1,1], limits:[-1d,10000d], step:0d},  $    ;systemic velocity
        {fixed:n_elements(veldisp)  , limited:[1,1], limits:[-3000d,3000d], step:0d},  $   ;velocity dispersion
        {fixed:n_elements(veldisp_sys)  , limited:[1,1], limits:[-1d,10000d], step:0d}   ]  ;systemic vdisp

  ENDELSE
  
  
  
ENDIF

END