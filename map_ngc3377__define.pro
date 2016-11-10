

;-----------------------------------------------------------------------------------------------

pro map_ngc3377::get_variogram_info, vs_fitting
	compile_opt idl2, hidden
	d = self->get_smeag(vs_fitting)
	d = struct_trimtags(d, select=['X','Y','VEL','VELDISP','H3','H4','RMS'])
	dummy = execute('d = struct_addtags(d, arr_struct({z:d.'+self.kinType+'}))')	
	case self.kinType of 
		'VEL'		: 	BEGIN
							svar_info = {range:100., ftype:'POLY', axisratio:0.54, inc:3}
						END
		'VELDISP'	: 	BEGIN
							svar_info = {range:200., ftype:'POLY', axisratio:1.0d, inc:5}
						END
		'RMS'		: 	BEGIN
							svar_info = {range:200., ftype:'POLY', axisratio:1.0d, inc:5}
						END
		'H3'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:90., ftype:'POLY', axisratio:1.0d, inc:10}
						END
		'H4'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:90., ftype:'POLY', axisratio:1.0d, inc:3}
						END
		else		: stop
	endcase	
	if (n_elements(self.maxdist) eq 0) then self.maxDist = svar_info.range*0.5
	svar_info = struct_addtags(svar_info, {maxDist:self.maxDist})
	svar_info = struct_addtags(svar_info, {data:d})
	*(self.svar_info) = svar_info
end

;-----------------------------------------------------------------------------------------------

pro map_ngc3377::get_gal_info,vs_fitting
	compile_opt idl2, hidden
	self.name = 'NGC3377'
	self.reff = 33.0 ;32.8
	self.distance = 10.9d6
	self.arc2kpc = self.distance/206265d/1000d
	self.q = 0.54
	self.parot = 46.3;41.3d
	;self.parot = 46d
	;self.parot = 52d
	self.galra = 161.9263808d  
	self.galdec = 13.9859158d
	self->get_variogram_info,vs_fitting
end

;-----------------------------------------------------------------------------------------------

function map_ngc3377::get_PNe
	compile_opt idl2, hidden
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC3377/N3377_PNdata1a.txt', id, x, y, r, pa, a, a, vel, dv, a, format='d,d,d,d,d,d,d,d,d,d', /silent
	d = arr_struct({ra:x/3600d + self.galra, dec:y/3600d + self.galdec, vel:vel - 686.5, errvel:dv})
	d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))
;  d.x = -d.x
	d = d[where(abs(d.vel) le 350)]
	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc3377::get_GCs, red=red, blue=blue
	compile_opt idl2, hidden	
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC3377/NGC3377_FINAL.dat', id, ra, dec, vel, errvel, g, gerr, r, rerr, i, ierr, format='a,d,d,d,d,d,d,d,d,d,d', /silent
;	d = arr_struct({ra:ra, dec:dec, vel:vel - median(vel,/even), errvel:errvel})
  d = arr_struct({ra:ra, dec:dec, vel:vel - 686.5, errvel:errvel})
	d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))	
;  d.x = -d.x
  if (n_elements(red) gt 0) then d = d[where((g - i) gt 1.05)] else if (n_elements(blue) gt 0) then d = d[where((g - i) le 1.05)]
	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc3377::get_smeag, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
;	d = self->get_skims('E:/Discrete Tracers/NGC3377/obj1_n3377.fits')
;	d = self->get_skims('E:/Discrete Tracers/NGC3377/comb_deimos_ppxf_v0_N3377.fits');, medianing=medianing)
  d = self->get_skims('I:/Brutus_InternalDrive/Discrete Tracers/NGC3377/deimos_ppxf_NGC3377.fits',vs_fitting);, medianing=medianing)
	;d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))
	d = d[where( (d.sn gt 4.5) )]
	d = d[where( (d.r lt 250) )]
	d.errvel = 35d
	d.errveldisp = 35d
;	d.errh3 = 0.01
;	d.errh4 = 0.01
;d.vel = sqrt(d.vel^2d + d.veldisp^2d)
;d.errvel *= sqrt(2d)
;plot, d.x, d.vel, psym=1
;pause
  d.x = -d.x

  d=struct_trimtags(d,except=['PA'])
  d = struct_addtags(d, arr_struct({pa:pacalc(d.x, d.y, self.parot,/radians)}))


	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc3377::get_coccato_fors
	compile_opt idl2, hidden	
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC3377/N3377_FORSdata.txt', id, pa, r, vel, dv, sig, dsig, h3, dh3, h4, dh4, /silent, format='d,d,d,d,d,d,d,d,d,d,d'			
	x = fltarr(n_elements(id))
	y = x
	x[(ind = where(pa eq 35, complement=ind1))] = r[ind]
	y[ind1] = r[ind1]
	x[0:130] = -x[0:130]
	vel[0:130] = vel[0:130] + 15.
	vel[131:258] = vel[131:258] - 15.
	vel[131:258] = -vel[131:258]
	y[259:389] = -y[259:389]
	ls = arr_struct({xi:-x, yi:y, pa:pa, vel:vel, errvel:dv, veldisp:sig, errveldisp:dsig, h3:h3, errh3:dh3, h4:h4, errh4:dh4})
;	xyrot, ls.xi, ls.yi, (35d - 90d)*!dtor, xp, yp			;35d - self.parot
  xyrot, ls.xi, ls.yi, (self.parot-90d)*!dtor, xp, yp      ;35d - self.parot

	ls = struct_addtags(ls, arr_struct({x:xp, y:yp}))
	ls = struct_addtags(ls, arr_struct({ra:self.galra - ls.x/3600d, dec:self.galdec + ls.y/3600d}))
	hold = arr_struct(projdist2(ls.ra, ls.dec, self.galra, self.galdec, self.parot))
	ls.x = hold.x
	ls.y = hold.y
	ls.x = -ls.x

  ls=struct_trimtags(ls,except=['PA'])
  ls = struct_addtags(ls, arr_struct({pa:pacalc(ls.x, ls.y, self.parot,/radians)}))
  

;ls.vel = sqrt(ls.vel^2d + ls.veldisp^2d)
;ls.errvel *= sqrt(2d)
;plot, ls.x, ls.vel, psym=1
;pause

	self->calc_rms, ls
	return, ls
end

;-----------------------------------------------------------------------------------------------

function map_ngc3377::get_sauron, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
	d = self->map_galaxy::get_sauron(pxf = 'I:/Brutus_InternalDrive/Discrete Tracers/NGC3377/PXF_bin_MS_NGC3377_r1_idl.fits', $
									ms='I:/Brutus_InternalDrive/Discrete Tracers/NGC3377/MS_NGC3377_r1.fits',vs_fitting);, medianing=medianing)	
;d.vel = sqrt(d.vel^2d + d.veldisp^2d)
;d.errvel *= sqrt(2d)
;plot, d.x, d.vel, psym=1
;pause
  d=struct_trimtags(d,except=['PA'])
  d.x = -d.x
  d = struct_addtags(d, arr_struct({pa:pacalc(d.x, d.y, self.parot,/radians)}))
	return, d
end

;-----------------------------------------------------------------------------------------------
PRO map_ngc3377::rolling_bins_kinemetry, d, nodiscrete=nodiscrete, nostellar=nostellar, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
    nopne=nopne, nogcs=nogcs, nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, N0821, N1023, N1344, N2768, N3115, N3377, N4564, N4697, N4473, inplotdat, inmaxr, bin1, l, nbootstraps, twobin=twobin,$
    obser, peakplot=peakplot, qin, dqin, nbootstrapsq
  compile_opt idl2, hidden
   
    diagnostics=1
    data = d
    reff = self.reff
    vsys = self.vsys
    name = self.name
    smeagol = 0 
    sauron = 0
    vs_fitting = vs_fitting

if n_elements(nostellar) eq 1 then outer = 1
if n_elements(nostellar) eq 1 then inner = 0
if n_elements(nodiscrete) eq 1 then outer = 0
if n_elements(nodiscrete) eq 1 then inner = 1

 TwoBin_Kinemetry, data, reff, vsys, diagnostics, name, inner, outer, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
    nopne=nopne, nogcs=nogcs, vs_fit, dvs_fit, sauron, smeagol, N0821, N1023, N1344, N2768, N3115, N3377, N4473, N4564, N4697, inplotdat, inmaxr, bin1, l, nbootstraps, twobin=twobin,$
    obser, peakplot=peakplot, qin, dqin, nbootstrapsq


END
;-----------------------------------------------------------------------------------------------
pro map_ngc3377::get_data, d, nopne=nopne, nogcs=nogcs, nodiscrete=nodiscrete, nostellar=nostellar, $
	nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, red=red, blue=blue, medianing=medianing,vs_fitting=vs_fitting
	compile_opt idl2, hidden

	if ((n_elements(nopne) eq 1) and (n_elements(nogcs) eq 1)) then nodiscrete = 1
	if (n_elements(nodiscrete) eq 0) then begin
		if (n_elements(nopne) eq 0) then discrete = self->get_PNe()								;retrieve the PN data
		if (n_elements(nogcs) eq 0) then discrete = struct_append(discrete, self->get_GCs(red=red, blue=blue))	;retrieve the GC data
		
		if 0 then begin
			print, '/Users/jaaarnol/NGC3377, using point symmetry on the bottom left corner!!!!'
			hold = discrete
			hold = hold[where((hold.x/self.reff gt 4.5) and (hold.y/self.reff gt 0.45))]
			hold.x = -1d*hold.x
			hold.y = -1d*hold.y
			hold.vel = -1d*hold.vel
			discrete = struct_append(discrete, hold)
		endif
		*(self.discrete) = discrete									;put the discrete data in a pointer
	endif
	
	if ((n_elements(nosauron) eq 1) and (n_elements(nolong) eq 1) and (n_elements(nosmeag) eq 1)) then nostellar = 1
	if (n_elements(nostellar) eq 0) then begin
		if (n_elements(nosauron) eq 0) then stellar = self->get_sauron(vs_fitting);medianing=medianing)		
		if (n_elements(nolong) eq 0) then begin
			ls = struct_append(ls, self->get_coccato_fors())
			if ((self.kinType ne 'H3') and (self.kinType ne 'H4')) then begin
			endif
			stellar = struct_append(stellar, ls)
		endif
		if (n_elements(nosmeag) eq 0) then stellar = struct_append(stellar, self->get_smeag(vs_fitting));medianing=medianing))
		*(self.stellar) = stellar				
		s=stellar
	endif
    if n_elements(nosauron) eq 0 or n_elements(nosmeag) eq 0 or n_elements(nolong) eq 0 or n_elements(nostellar) eq 0 then d = struct_selecttags(s,select=['x','y','r','pa','vel','errvel']) ;subarray to hold stellar+discrete data (start by adding stellar data)
    if n_elements(nodiscrete) eq 0 then begin
      dummy = struct_selecttags(discrete,select=['x','y','r','pa','vel','errvel']) 
      d = struct_append(d,dummy)
    endif
end

;-----------------------------------------------------------------------------------------------
pro map_ngc3377::isophote, p2
  compile_opt idl2, hidden

  s = *(self.plot)
  pos = self.pos
  
nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
q = 0.67 ; 2Re
;maxa = 1d/sqrt(q)*s.xsiz/float(difference(minmax(s.xrange)))
;maxa = 2*self.reff/sqrt(1-q)
maxa = (4/7d)*(max(s.xsiz)/2)
xpa = 0d*!dtor

i = 0
  pa = dindgen(nbound)*!dpi*2d/nbound - !dpi
  a = dblarr(nbound) + maxa[i]
  xhigh = a * sign(cos(pa))/sqrt(1d + tan(pa)^2d / q^2d)
  yhigh = -xhigh * tan(pa)
  xh = xhigh*cos(xpa) - yhigh*sin(xpa)
  yhigh = xhigh*sin(xpa) + yhigh*cos(xpa)
  xhigh = xh + xmid
  yhigh = yhigh + ymid
  
  
  xhigh1 = [xhigh, xhigh[0]]
  yhigh1 = [yhigh, yhigh[0]]
  
nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
q = 0.54 ; 1Re
;maxa = 1d/sqrt(q)*s.xsiz/float(difference(minmax(s.xrange)))
;maxa = 156
xpa = 0d*!dtor

i = 0
  pa = dindgen(nbound)*!dpi*2d/nbound - !dpi
  a = dblarr(nbound) + maxa[i]
  xhigh = a * sign(cos(pa))/sqrt(1d + tan(pa)^2d / q^2d)
  yhigh = -xhigh * tan(pa)
  xh = xhigh*cos(xpa) - yhigh*sin(xpa)
  yhigh = xhigh*sin(xpa) + yhigh*cos(xpa)
  xhigh = xh + xmid
  yhigh = yhigh + ymid

  isophote2 = 4.
  
  outerxhigh = xhigh/isophote2 + (xmid - xmid/isophote2)
  outeryhigh = yhigh/isophote2 + (ymid - ymid/isophote2)
  
  xhigh2 = [outerxhigh, outerxhigh[0]]
  yhigh2 = [outeryhigh, outeryhigh[0]]

  p2 = plot(xhigh1, yhigh1, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)
  p2 = plot(xhigh2, yhigh2, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)

  

  
end
;-----------------------------------------------------------------------------------------------
;pro map_ngc3377::isophote, p2
;  compile_opt idl2, hidden
;
;  s = *(self.plot)
;  pos = self.pos
;
;	files = 'E:/Discrete Tracers/NGC3377/n3377_'+roundx(s.isophote_level,0)+'Re.dat'
;
;	k = 0
;	foreach confile, files do begin
;		readcol, confile, con, format='d', /silent
;		ncon = n_elements(con)
;		racon = con[(indgen(ncon)*2)[0:ncon/2-1]]
;		deccon = con[(indgen(ncon)*2 + 1)[0:ncon/2-1]]
;		;con = arr_struct(projdist2(racon, deccon, galra, galdec, xpa))
;		con = arr_struct({xi:(racon - self.galra)*3600d, yi:(deccon - self.galdec)*3600d})
;		con = struct_append(con, con[0])
;		xyrot, con.xi, con.yi, (self.parot - 90d)*!dtor, xp, yp
;		con = struct_addtags(con, arr_struct({x:xp, y:yp}))
;		if (k eq 0) then smcon = 0 else smcon = 15
;		if self.doreff then begin
;			con.x /= self.reff
;			con.y /= self.reff
;		endif else begin
;			con.x *= self.arc2kpc
;			con.y *= self.arc2kpc
;		endelse
;		con.x = smooth(con.x*s.xsiz/float(difference(minmax(s.xrange))) + self.xmid, smcon);, /edge)
;		con.y = smooth(con.y*s.xsiz/float(difference(minmax(s.xrange))) + self.ymid, smcon);, /edge)
;		con = struct_append(con, con[0])
;;		p2 = plot(con.x, con.y, position=pos, /overplot, thick=s.c_thick)
;		++k
;	endforeach
;end

;-----------------------------------------------------------------------------------------------

pro map_ngc3377__define
	compile_opt idl2, hidden

	void={map_ngc3377, inherits map_galaxy}

end






 
 
