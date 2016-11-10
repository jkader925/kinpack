;-----------------------------------------------------------------------------------------------

pro map_ngc2768::get_variogram_info,vs_fitting
	compile_opt idl2, hidden
	d = self->get_smeag(vs_fitting)
	d = struct_trimtags(d, select=['X','Y','VEL','VELDISP','H3','H4','RMS'])
	dummy = execute('d = struct_addtags(d, arr_struct({z:d.'+self.kinType+'}))')	
	case self.kinType of 
		'VEL'		: 	BEGIN
							svar_info = {range:280., ftype:'POLY', axisratio:0.42d, inc:5}
						END
		'VELDISP'	: 	BEGIN
							svar_info = {range:150., ftype:'LINEAR', axisratio:1.0d, inc:10}
						END
		'RMS'		: 	BEGIN
							svar_info = {range:150., ftype:'LINEAR', axisratio:1.0d, inc:10}
						END
		'H3'		: 	BEGIN
							svar_info = {range:220., ftype:'POLY', axisratio:1.0d, inc:4}
						END
		'H4'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:120., ftype:'POLY', axisratio:1.0d, inc:4}
						END
		else		: stop
	endcase	
	if (n_elements(self.maxdist) eq 0) then self.maxDist = svar_info.range*0.5
	svar_info = struct_addtags(svar_info, {maxDist:self.maxDist})
	svar_info = struct_addtags(svar_info, {data:d})
	*(self.svar_info) = svar_info
end

;-----------------------------------------------------------------------------------------------

pro map_ngc2768::get_gal_info,vs_fitting
	compile_opt idl2, hidden
	self.name = 'NGC2768'
	self.reff = 64.
	self.distance = 20.1d6
	self.q = 0.42
	self.parot = 89.9;91.6d (value is by default restricted by kinemetry_max parameter info limits)
	self.galra = 137.906250d  
	self.galdec = 60.037222d
	self.arc2kpc = self.distance/206265d/1000d
	self->get_variogram_info,vs_fitting
end

;-----------------------------------------------------------------------------------------------

function map_ngc2768::get_smeag, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
;	d = self->get_skims('I:/Discrete Tracers/NGC2768/NGC2768_skims.fits')
;	d = self->get_skims('I:/Discrete Tracers/NGC2768/comb_deimos_ppxf_v0.fits', medianing=medianing)
 d = self->get_skims('I:/Discrete Tracers/NGC2768/deimos_ppxf_NGC2768.fits',vs_fitting);, medianing=medianing)

	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc2768::get_illingworth_long
  compile_opt idl2, hidden
  readcol, 'I:/Discrete Tracers/NGC2768/illingworth_long.txt', r, pa, vel, errvel, veldisp, errveldisp
  relPA = 90d - self.parot
  d = arr_struct({r:r, pa:pa*0d, dx:r*0d, dy:r*0d, vel:vel, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp});what v_sys did illingworth use?
  d.pa = relPA

  d.dx = d.r*cos((d.pa)*!dtor)
  d.dy = d.r*sin((d.pa)*!dtor)

;  ind = where(d.pa ne 0)
;  d1 = d[where(d.pa ne 0)]
;  d1.dx = -d1.dx
;  d1.dy = -d1.dy
;  d = struct_append(d, d1)  
      
;  mdm = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
;  mdm = struct_addtags(mdm, struct_trimtags(d, except=['PA','R','SYM']))
;  mdm = struct_addtags(mdm, arr_struct({x:mdm.dx, y:mdm.dy}))
;;  mdm = struct_addtags(mdm, arr_struct(projdist2(mdm.ra, mdm.dec, self.galra, self.galdec, self.parot)))
;  self->calc_rms, mdm 
;  ;self->get_Flux, mdm
;  return, mdm
  het = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
  het = struct_addtags(het, struct_trimtags(d, except=['PA','R','DX','DY']))
  het = struct_addtags(het, arr_struct(projdist2(het.ra, het.dec, self.galra, self.galdec, self.parot)))
  het=struct_trimtags(het,except=['PA'])
  het.x = -het.x
  het = struct_addtags(het, arr_struct({pa:pacalc(het.x, het.y, self.parot,/radians)}))
  self->calc_rms, het
  ;;self->get_Flux, het
  
  return, het

end

;-----------------------------------------------------------------------------------------------

function map_ngc2768::get_PNe
  compile_opt idl2, hidden
  
  readcol, 'I:/Discrete Tracers/NGC2768/n2768_pndata.txt', ra, dec, vel, velerr, format='d,d,d,d', /silent
  d = arr_struct({ra:ra, dec:dec, vel:vel - 1376.5, errvel:velerr})
  d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))  

  return, d
end

;-----------------------------------------------------------------------------------------------
function map_ngc2768::get_sauron,medianing=medianing,vs_fitting
	compile_opt idl2, hidden
	d = self->map_galaxy::get_sauron(pxf = 'I:/Discrete Tracers/NGC2768/PXF_bin_MS_NGC2768_r5_idl.fits', $
									ms='I:/Discrete Tracers/NGC2768/MS_NGC2768_r5.fits',vs_fitting);, medianing=medianing)

	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc2768::get_GCs, red=red, blue=blue
  compile_opt idl2, hidden  
;readcol, 'I:/Discrete Tracers/NGC2768/N2768gcvels.txt', ra, dec, vel, errvel, Rc, Rcerr, i, ierr, z, zerr, format='d,d,d,d,d,d,d,d,d,d', /silent ;Pota et al. 2012 catalog of GCs
readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC2768/n2768_forbesetal2012.txt', ra, dec, r, pa, vel, verr, sig, sigerr,h3, h3err, h4, h4err, format='d,d,d,d,d,d,d,d,d,d,d,d', /silent ;Pota et al. 2012 catalog of GCs
  d = arr_struct({ra:ra, dec:dec, vel:vel - 1376.5, errvel:verr})
  d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))  
  if (n_elements(red) gt 0) then d = d[where((g - i) gt 0.88)] else if (n_elements(blue) gt 0) then d = d[where((g - i) le 0.88)]

  return, d
end

;-----------------------------------------------------------------------------------------------

pro map_ngc2768::get_data, d, nopne=nopne, nogcs=nogcs, nodiscrete=nodiscrete, nostellar=nostellar, $
	nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, red=red, blue=blue, medianing=medianing,vs_fitting=vs_fitting
	compile_opt idl2, hidden

	if ((n_elements(nopne) eq 1) and (n_elements(nogcs) eq 1)) then nodiscrete = 1
	if (n_elements(nodiscrete) eq 0) then begin
		if (n_elements(nopne) eq 0) then discrete = self->get_PNe()								;retrieve the PN data
		if (n_elements(nogcs) eq 0) then discrete = struct_append(discrete, self->get_GCs(red=red, blue=blue))	;retrieve the GC data
		*(self.discrete) = discrete									;put the discrete data in a pointer
	endif
	
	if ((n_elements(nosauron) eq 1) and (n_elements(nolong) eq 1) and (n_elements(nosmeag) eq 1)) then nostellar = 1
	if (n_elements(nostellar) eq 0) then begin
		if (n_elements(nosauron) eq 0) then stellar = self->get_sauron(vs_fitting);medianing=medianing)		
		if (n_elements(nolong) eq 0) then begin
		  ls = struct_append(ls, self->get_illingworth_long())
			stellar = struct_append(stellar, ls)
		endif
    if (n_elements(nosmeag) eq 0) then stellar = struct_append(stellar, self->get_smeag(vs_fitting));medianing=medianing))
    *(self.stellar) = stellar
    s = stellar
      svar_info = *(self.svar_info)			
	endif	
    if n_elements(nosauron) eq 0 or n_elements(nosmeag) eq 0 or n_elements(nolong) eq 0 or n_elements(nostellar) eq 0 then d = struct_selecttags(s,select=['x','y','r','pa','vel','errvel']) ;subarray to hold stellar+discrete data (start by adding stellar data)
    if n_elements(nodiscrete) eq 0 then begin
      dummy = struct_selecttags(discrete,select=['x','y','r','pa','vel','errvel']) 
      d = struct_append(d,dummy)
    endif
end


;-----------------------------------------------------------------------------------------------

PRO map_ngc2768::rolling_bins_kinemetry, d, nodiscrete=nodiscrete, nostellar=nostellar, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
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
pro map_ngc2768::isophote, p2
  compile_opt idl2, hidden

  s = *(self.plot)
  pos = self.pos

nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
q = 0.48 ; outer isophote flattening
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
q = 0.42  ;inner isophote flattening
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
  
  outerxhigh = xhigh/isophote2 + (xmid - xmid/isophote2)  ;inner isophote!
  outeryhigh = yhigh/isophote2 + (ymid - ymid/isophote2)  ;inner isophote!
  
  xhigh2 = [outerxhigh, outerxhigh[0]]
  yhigh2 = [outeryhigh, outeryhigh[0]]

  p2 = plot(xhigh1, yhigh1, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)
  p2 = plot(xhigh2, yhigh2, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)

  

  
end

;-----------------------------------------------------------------------------------------------
pro map_ngc2768__define
	compile_opt idl2, hidden

	void={map_ngc2768, inherits map_galaxy}

end






 
 
