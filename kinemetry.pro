;+
;PURPOSE
;	To fit rotational amplitude and calculate velocity dispersion about that fit
;	for a dataset consisting of (position angle, velocity) pairs
;SYNTAX
;	kinemetry(structure{pa, vel}, mpfit parameters, initial guess, weigthing)
;Written by Jacob A. Arnold, 8-22-09, UCSC
;-

function kinemetry, pa, quan, err, pi, start, key=key, res=res, quiet=quiet

	if (n_elements(pi) lt n_elements(start)) then pi1 = struct_append(pi,pi[1]) else pi1 = pi

	result = mpfitfun('kinem_func', pa, quan, err, start, parinfo=pi1, perror=perror, dof=dof, bestnorm=bestnorm, quiet=quiet)
	
	;xpar = findgen(360)*!pi/180d
	;ypar = kinem_func(xpar, start[0:4])
	;plot, xpar, ypar, psym=1
	
	if (n_elements(result) eq 1) then stop
	
	
	;calculate a dispersion about the solution
	;disp = sqrt(total((quan - kinem_func(pa, [result[0],result[1],result[2],result[3]]))^2d)/(n_elements(pa) - 1l))
	disp = robust_sigma(quan - kinem_func(pa, [result[0],result[1],result[2],result[3]]), /zero)			;ultimately not used for anything?
	if (n_elements(key) gt 0) then begin
		plot, data.pa, data.vel, psym=8
		xx = findgen(100)/100.*2*!dpi
		oplot, xx, kinem_func(xx, [result[0],result[1],result[2],result[3]])
		stop
		pause
	endif
	
vmod = kinem_func(pa, [result[0],result[1],result[2],result[3]])

linex=indgen(7)
liney=indgen(7)*0.
pp=errorplot(pa,v,errv,linestyle=6,symbol='o',sym_filled=1,color='black',yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='Vrot [km/s]',sym_size=0.3,title='NGC 821')
pp=plot(linex,liney,linestyle=2,color='black',/overplot)
pp=plot(pa[sort(pa)],vmod[sort(pa)],color='red',linestyle=0,thick=2,/overplot)
pp=plot(pa[sort(pa)],vmod[sort(pa)]+2.*p[4],color='red',linestyle=1,thick=2,/overplot)
pp=plot(pa[sort(pa)],vmod[sort(pa)]-2.*p[4],color='red',linestyle=1,thick=2,/overplot)
t=text(4.5,350.,strcompress('Vrot='+string(round(p[0]))+'km/s'),color='orange red',target=pp,/data)
t=text(4.5,300.,strcompress('kPA='+string(round(p[1]/!dtor))+'degrees'),color='orange red',target=pp,/data)
t=text(4.5,250.,strcompress('Sigma='+string(round(p[4]))+'km/s'),color='orange red',target=pp,/data)
t=text(0.25,300,'~1 R_e',target=pp,/data)
	
	stop
	
	
	res = result
	return, {quan:result[0], pa:result[1], axisratio:result[2], sys:result[3], disp:disp, bestnorm:bestnorm, dof:dof}
	
end

















