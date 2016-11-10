;+
;PURPOSE
; To fit rotational amplitude and calculate velocity dispersion about that fit
; for a dataset consisting of (position angle, velocity) pairs
;SYNTAX
; kinemetry(structure{pa, vel}, mpfit parameters, initial guess, weigthing)
;Written by Jacob A. Arnold, 8-22-09, UCSC
;-

function kinemetry2, pa, quan, err, pi, start, key=key, res=res, quiet=quiet, diag=diag, sige, er

for nn = 0,1 do begin

  if (n_elements(pi) lt n_elements(start)) then pi1 = struct_append(pi,pi[1]) else pi1 = pi

  result = mpfitfun('kinem_func', pa, quan, err, start, parinfo=pi1, perror=perror, dof=dof, bestnorm=bestnorm, quiet=quiet, maxiter = 1000)
  
  xpar = findgen(360)*!pi/180d
  ypar = kinem_func(xpar, [result[0],result[1],result[2],result[3]])
  linex=indgen(7)
  liney=indgen(7)*0
  plot, xpar, ypar, psym=1,xrange=[0,6.2],yrange=[-500,500]
  oplot,pa,quan,psym=5,color=fsc_color('lime green')
  oplot,linex,liney,linestyle=1,color=fsc_color('dodger blue')
  xyouts,0.5,-400,strcompress('Vrot/Sigma ='+string(result[0]/sige)),color=fsc_color('red'),/data
  xyouts,0.5,-300,strcompress('R/R_e ='+string(er)),color=fsc_color('red'),/data
;  erase
; stop
  if (n_elements(result) eq 1) then stop
  
  
  ;calculate a dispersion about the solution
  ;disp = sqrt(total((quan - kinem_func(pa, [result[0],result[1],result[2],result[3]]))^2d)/(n_elements(pa) - 1l))
  disp = robust_sigma(quan - kinem_func(pa, [result[0],result[1],result[2],result[3]]), /zero)      ;ultimately not used for anything?
  
  if (n_elements(key) gt 0) then begin
    plot, data.pa, data.vel, psym=8
    xx = findgen(100)/100.*2*!dpi
    oplot, xx, kinem_func(xx, [result[0],result[1],result[2],result[3]])
    pause
  endif
  
  
  
vmod = kinem_func(pa, [result[0],result[1],result[2],result[3]])

;error rejection
;---------------------------------------------------

;holder = fltarr(n_elements(quan))
;for i = 0,n_elements(quan)-1 do begin
;    holder[i] = abs(quan[i] - vmod[i])
;endfor
;quan1 = quan
;pa1 = pa
;err1 = err
;
;quan = quan[where(holder LT 1.5*disp,complement=rquan)]
;pa = pa[where(holder LT 1.5*disp,complement=rpa)]
;err = err[where(holder LT 1.5*disp,complement=rerr)]
;
;vmod = kinem_func(pa, [result[0],result[1],result[2],result[3]])

;---------------------------------------------------



;if nn eq 0 then pp=errorplot(pa,quan,err,linestyle=6,symbol='o',sym_filled=1,color='lime green',sym_size=1,/overplot)
;diag = 1
;nn = 0
if diag eq 1 and nn eq 0 then begin 
linex=indgen(7)
liney=indgen(7)*0.
window,2,xsize=500,ysize=500
plot,pa,quan,psym=5,yrange=[-400d,400d],xrange=[0d,6d],xtitle='PA [rad]',ytitle='V_rot [km/s]'
oplot,linex,liney,linestyle=2,color=fsc_color('dodger blue')
oplot,pa[sort(pa)],vmod[sort(pa)],color=fsc_color('dodger blue'),linestyle=0,thick=2
;oplot,pa1[rpa],quan1[rpa],psym=5,color=fsc_color('red')
stop
endif









endfor
  
  
  
;linex=indgen(7)
;liney=indgen(7)*0.
;pp=errorplot(pa,quan,err,linestyle=6,symbol='o',sym_filled=1,color='black',yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='Vrot [km/s]',sym_size=0.3)
;;;pp=errorplot(pa,quan,err,linestyle=6,symbol='o',sym_filled=1,color='black',xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='Vrot [km/s]',sym_size=0.3,title='NGC 821')
;pp=plot(linex,liney,linestyle=2,color='black',/overplot)
;pp=plot(pa[sort(pa)],vmod[sort(pa)],color='dodger blue',linestyle=0,thick=2,/overplot)
;pp=plot(pa[sort(pa)],vmod[sort(pa)]+2d*disp,color='dodger blue',linestyle=3,thick=2,/overplot)
;pp=plot(pa[sort(pa)],vmod[sort(pa)]-2d*disp,color='dodger blue',linestyle=3,thick=2,/overplot)
;t=text(4.5,350.,strcompress('Vrot='+string(round(result[0]))+'km/s'),color='dodger blue',target=pp,/data)
;t=text(4.5,300.,strcompress('V/sig='+string(result[0]/sige)),color='dodger blue',target=pp,/data)
;t=text(4.5,250.,strcompress('kPA='+string(round(result[1]/!dtor))+'degrees'),color='dodger blue',target=pp,/data)
;;t=text(4.5,200.,strcompress('Sigma='+string(round(result[4]))+'km/s'),color='dodger blue',target=pp,/data)
;t=text(0.25,300,'~1 R_e',target=pp,/data)
;stop
  
  
  
  res = result
  return, {quan:result[0], pa:result[1], axisratio:result[2], sys:result[3], disp:disp, bestnorm:bestnorm, dof:dof}
  
end

















