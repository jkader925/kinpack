;04/15/15 -- multiply V_rad by a factor of 1.5 to reduce underestimation of v_rot




;+
;PURPOSE
; To fit rotational amplitude and calculate velocity dispersion about that fit
; for a dataset consisting of (position angle, velocity) pairs
;SYNTAX
; kinemetry(structure{pa, vel}, mpfit parameters, initial guess, weigthing)
;Written by Jacob A. Arnold, 8-22-09, UCSC
;-



function maxkin, X, df, _EXTRA=extra
common share, pangle, v, errv, lambda, sigcorr

vrot = x[0]
pa = x[1]
q = x[2]
vsys = x[3]
sig = x[4]
sig_sys = x[5]
;h3 = x[6]
;h3_sys = x[7]
;h4 = x[8]
;h4_sys = x[9]
vmod = kinem_func(pangle, [vrot, pa, q, vsys])
;h3mod =  kinem_func(pangle, [h3, pa, q, h3_sys]) 
;h4mod = kinem_func(pangle, [h4, pa, q, h4_sys])
;sig = kinem_func(pangle, [sig, pa, q, sig_sys])

sig = sig / sigcorr
;.......................plot fit on top of (pa,vel) measurements
;linex=-indgen(7)
;liney=indgen(7)*0.
;plot,pangle,v,psym=5,yrange=[-400d,400d],xrange=[0d,6d],xtitle='PA [rad]',ytitle='V_rot [km/s]'
;oplot,linex,liney,linestyle=2,color=fsc_color('dodger blue')
;oplot,pangle[sort(pangle)],vmod[sort(pangle)],color=fsc_color('lime green'),linestyle=0,thick=2
;oplot,pangle[sort(pangle)],vmod[sort(pangle)]+2.*sig,color=fsc_color('lime green'),linestyle=1
;oplot,pangle[sort(pangle)],vmod[sort(pangle)]-2.*sig,color=fsc_color('lime green'),linestyle=1
;xyouts,-7.5,250,strcompress('PA ='+string(PA/!dtor)+'degrees'),color=fsc_color('lime green'),/data
;xyouts,-7.5,300,strcompress('vrot ='+string(vrot)),color=fsc_color('lime green'),/data
;xyouts,-7.5,350,strcompress('Sigma ='+string(sig)),color=fsc_color('lime green'),/data
;erase
;p=plot(pangle,v,linestyle=6,symbol='o',sym_filled=1,color='black',yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='V_LOS [km/s]',sym_size=0.1)
;p=plot(linex,liney,linestyle=2,color='black',/overplot)
;p=plot(pangle[sort(pangle)],vmod[sort(pangle)],color='red',linestyle=0,thick=2,/overplot)
;p=plot(pangle[sort(pangle)],vmod[sort(pangle)]+2.*sig,color='red',linestyle=1,thick=2,/overplot)
;p=plot(pangle[sort(pangle)],vmod[sort(pangle)]-2.*sig,color='red',linestyle=1,thick=2,/overplot)
;t=text(4.5,300.,strcompress('Vrot='+string(round(vrot))+'km/s'),color='orange red',target=p,/data)
;t=text(4.5,250.,strcompress('Sigma='+string(round(sig))+'km/s'),color='orange red',target=p,/data)
;stop
;.......................

;F = total(  0.5d*( (((v - vmod)^2d) / (sig^2d + errv^2d)) + alog(2d*!dpi) + alog(sig^2d + errv^2d) ) - alog(1d + h3mod*gausshermite((v - vmod)/sqrt(sig^2d + errv^2d),3) + h4mod*gausshermite((v - vmod)/sqrt(sig^2d + errv^2d),4)) )
;F = total(  0.5d*( (((v - vmod)^2d) / (sig^2d + errv^2d)) + alog(2d*!dpi) + alog(sig^2d + errv^2d) ))
;F = total(  0.5d*( (((v - vmod)^2d) / (sig^2d + errv^2d)) + alog(sig^2d + errv^2d) ))
;F = total(  0.5d*( (((v - vmod)^2d) / (sig^2d + errv^2d)) + alog(sig^2d + errv^2d) ))
F = total(  ( (((v - vmod)^2d) / (sig^2d + errv^2d)) + alog(sig^2d + errv^2d) ))

;f = f*(1d + (lambda^2d)*(h3^2d + h4^2d))
if (finite(f) eq 0) then f = 1d10

if (q eq 0) then begin
  print, f
  f = 1d10
  ;stop
endif

return, F
END









function kinemetry_max2, pa, quan, err, pi, start, lam=lam, verbose=verbose, sige
common share

errv = err
v = quan
pangle = pa

;for nn = 0,1 do begin
npnts = double(n_elements(v))
sigcorr = kn(npnts)*(npnts/(npnts - 1d))^0.5d
sigcorr = 1.0037696

if (n_elements(lam) eq 0) then lam = 0
lambda = lam

if (n_elements(verbose) eq 0) then quiet = 1


;ALERT  ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT 
;ALERT  ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT 
;ALERT  ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT 

;this program can fail without any warning if you don't make sure that your starting values are 
; within the parinfo limits, need to check this here.
;Also, make sure there are no infinities or NaN

;ALERT  ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT 
;ALERT  ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT 
;ALERT  ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT ALERT 

diag = 1
FOR nn = 0,1 DO BEGIN

  p = tnmin('maxkin', start, bestmin=bestmin, autoderivative=1, nprint=1, parinfo=pi, status=status, quiet=1)

  if (n_elements(p) eq 1) then return, p

  result = {quan:p[0], pa:p[1], axisratio:p[2], sys:p[3], disp:p[4], $
    disp_sys:p[5], err_vrot:errv,   $
    bestnorm:bestmin, dof:n_elements(pa)-n_elements(where(pi.fixed eq 0))}


  vmod = kinem_func(pa, [P[0],P[1],P[2],P[3]])


;error rejection
;---------------------------------------------------

;  holder = fltarr(n_elements(quan))
;  for i = 0,n_elements(quan)-1 do begin
;    holder[i] = abs(quan[i] - vmod[i])
;  endfor
;  quan1 = quan
;  pa1 = pa
;  err1 = err
;
;  quan = quan[where(holder LT 1.5*P[4],complement=rquan)]
;  pa = pa[where(holder LT 1.5*P[4],complement=rpa)]
;  err = err[where(holder LT 1.5*P[4],complement=rerr)]
;
;  vmod = kinem_func(pa, [P[0],P[1],P[2],P[3]])


;  if diag eq 1 and nn eq 0 then begin 
;  linex=indgen(7)
;  liney=indgen(7)*0.
;  window,2,xsize=500,ysize=500
;  plot,pa,quan,psym=5,yrange=[-400d,400d],xrange=[0d,6d],xtitle='PA [rad]',ytitle='V_rot [km/s]'
;  oplot,linex,liney,linestyle=2,color=fsc_color('dodger blue')
;  oplot,pa[sort(pa)],vmod[sort(pa)],color=fsc_color('dodger blue'),linestyle=0,thick=2
;  oplot,pa1[rpa],quan1[rpa],psym=5,color=fsc_color('red')
;  stop
;  endif
;---------------------------------------------------
  


  pangle = pa
  v = quan
  errv = err

ENDFOR





linex=indgen(7)
liney=indgen(7)*0.
pp=errorplot(pangle,v,errv,linestyle=6,symbol='o',sym_filled=1,color='black',yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='Vrot [km/s]',sym_size=0.3,title='Wet 1:1 Orbit i')
pp=plot(linex,liney,linestyle=2,color='black',/overplot)
pp=plot(pangle[sort(pangle)],vmod[sort(pangle)],color='red',linestyle=0,thick=2,/overplot)
pp=plot(pangle[sort(pangle)],vmod[sort(pangle)]+2.*p[4],color='red',linestyle=1,thick=2,/overplot)
pp=plot(pangle[sort(pangle)],vmod[sort(pangle)]-2.*p[4],color='red',linestyle=1,thick=2,/overplot)
t=text(4.5,350.,strcompress('Vrot='+string(round(p[0]))+'km/s'),color='red',font_size=16,target=pp,/data)
t=text(4.5,300.,strcompress('kPA='+string(round(p[1]/!dtor))+'degrees'),color='red',font_size=16,target=pp,/data)
t=text(4.5,250.,strcompress('Sigma='+string(round(p[4]))+'km/s'),color='red',font_size=16,target=pp,/data)
t=text(4.5,200.,strcompress('Q_kin='+string(p[2])),color='red',font_style='bold',font_size=16,target=pp,/data)
;t=text(4.5,200.,strcompress('V/sig_e='+string(p[0]/sige)),color='orange red',target=pp,/data)
;t=text(0.25,300,'~4 R_e',target=pp,/data)
;t=text(0.25,300,'~1 R_e',target=pp,/data)





return, result



end








