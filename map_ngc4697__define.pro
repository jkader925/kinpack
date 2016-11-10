;-----------------------------------------------------------------------------------------------

;PRO map_ngc4697::get_variogram_info
;	compile_opt idl2, hidden
;	svar_info = [$
;        {type:'VEL',        range:220., ftype:'POLY1',      axisratio:1.0d, inc:3, varplot:0}, 	$
;        {type:'VELDISP',    range:200., ftype:'GAUSSIAN',   axisratio:1.0d, inc:5, varplot:0}, 	$
;        {type:'RMS',        range:300., ftype:'POLY',       axisratio:1.0d, inc:3, varplot:1}, 	$
;        {type:'H3',         range:300., ftype:'GAUSSIAN',   axisratio:1.0d, inc:5, varplot:0}, 	$
;        {type:'H4',         range:300., ftype:'POLY1',      axisratio:1.0d, inc:3, varplot:0}, 	$
;        {type:'FLUX',       range:180., ftype:'POLY',       axisratio:1.0d, inc:3, varplot:0}	]
;	self->map_galaxy::get_variogram_info, svar_info
;	stop
;END
;
;-----------------------------------------------------------------------------------------------

pro map_ngc4697::get_variogram_info,vs_fitting
  compile_opt idl2, hidden
  d = self->get_smeag(vs_fitting)
  d = struct_trimtags(d, select=['X','Y','VEL','VELDISP','H3','H4','RMS','FLUX'])
  dummy = execute('d = struct_addtags(d, arr_struct({z:d.'+self.kinType+'}))')  
  case self.kinType of 
    'VEL'   :   BEGIN
              svar_info = {range:220., ftype:'POLY', axisratio:0.64d, inc:3}
            END
    'VELDISP' :   BEGIN
              svar_info = {range:200., ftype:'GAUSSIAN', axisratio:1.0d, inc:5}
            END
    'RMS'   :   BEGIN
              svar_info = {range:300., ftype:'POLY', axisratio:1.0d, inc:3}
            END
    'H3'    :   BEGIN
              svar_info = {range:300., ftype:'GAUSSIAN', axisratio:1.0d, inc:5}
            END
    'H4'    :   BEGIN
              svar_info = {range:300., ftype:'POLY1', axisratio:1.0d, inc:3}
            END
    'FLUX'    :   BEGIN
              svar_info = {range:180., ftype:'POLY', axisratio:1.0d, inc:3}
            END
    else    : stop
  endcase 
  stop
  if (n_elements(self.maxdist) eq 0) then self.maxDist = svar_info.range*0.5
  svar_info = struct_addtags(svar_info, {maxDist:self.maxDist})
  svar_info = struct_addtags(svar_info, {data:d})
  *(self.svar_info) = svar_info

end
;-----------------------------------------------------------------------------------------------

PRO map_ngc4697::get_gal_info,vs_fitting
	compile_opt idl2, hidden
	self.name = 'NGC4697'
	self.reff = 62d  ;Arnold 2014/Brodie 2014
	self.distance = 12.3d6
	self.q = 0.64 ;figureS6_Progress.ods
	self.parot = 67.2d   ;NED K_s (LGA/2MASS)
	self.galra = 192.1494908d      ;NED
	self.galdec = -5.8007419d      ;NED
  self.arc2kpc = self.distance/206265d/1000d
  self->get_variogram_info,vs_fitting


  
;	self.morph = -4.5              ;HyperLeda
; self->calc_galaxy_properties	
END

;-----------------------------------------------------------------------------------------------

FUNCTION map_ngc4697::get_GCs, red=red, blue=blue
	compile_opt idl2, hidden
	return, d
END

;-----------------------------------------------------------------------------------------------

FUNCTION map_ngc4697::get_smeag, nosyscorr=nosyscorr, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
;	file = 'I:/Discrete Tracers/NGC4697/comb_deimos_ppxf_NGC4697_v0.fits'
	file = 'I:/Discrete Tracers/NGC4697/deimos_ppxf_NGC4697.fits'
;	d = self->get_skims(file, nosyscorr=nosyscorr, maxVelDisp=maxVelDisp)	
  d = self->get_skims(file,vs_fitting);, medianing=medianing)
  d.x = -d.x
  d=struct_trimtags(d,except=['PA'])
  d = struct_addtags(d, arr_struct({pa:pacalc(d.x, d.y, self.parot,/radians)}))  
  
	return, d
END


;-----------------------------------------------------------------------------------------------

;FUNCTION map_ngc4697::get_deimos_longslit, srn=srn
;;	d = mrdfits(self.basedir+'NGC4697/longslit/deimos_ppxf_d0411_0092_0.fits',1,/silent)
;;	d1 = mrdfits(self.basedir+'NGC4697/longslit/deimos_ppxf_d0411_0095_0.fits',1,/silent)
;  d = mrdfits('I:/Discrete Tracers/NGC4697/N4697.kindata1',1,/silent)
;;  d1 = mrdfits('I:/Discrete Tracers/NGC4697/deimos_ppxf_d0411_0095_0.fits',1,/silent)
;;	if ~self.noSysCorr then begin
;;		d1.vel -= median(d.vel,/even)
;;		d.vel -= median(d.vel,/even)
;;	endif
;	d = struct_append(d,d1)
;	d = struct_trimtags(d, except=['PA','X','Y','R'])
;	d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))
;	if (self.kinType eq 'RMS') then begin
;		self->calc_rms, d
;	endif
;	;self->get_Flux, d
;	self->slitScatter, d, width=width
;	return, d
;END
;-----------------------------------------------------------------------------------------------

FUNCTION map_ngc4697::get_deimos_longslit, srn=srn
  readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC4697/N4697.kindata1',a, dx, dy, r, pa, b, c, dvel, vel, dd, e, f, format='d,d,d,d,d,d,d,d,d,d,d,d'
  d = arr_struct({dx:dx, dy:dy, r:r, pa:pa, vel:vel, dvel:dvel})
  het = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
  het = struct_addtags(het, struct_trimtags(d, except=['PA','R','DX','DY']))
  het = struct_addtags(het, arr_struct(projdist2(het.ra, het.dec, self.galra, self.galdec, self.parot)))
;  a=het[where(het.vel GE median(het.vel))]
;  b=het[where(het.vel LT median(het.vel))]
;  p=plot(a.x,a.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='red',title='N4697 LS')
;  p=plot(b.x,b.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='dodger blue',/overplot)
  ;self->get_Flux, het
  return, het
END

;-----------------------------------------------------------------------------------------------


function map_ngc4697::get_PNe
  compile_opt idl2, hidden

    readcol, 'I:/Discrete Tracers/NGC4697/Mendez.data1', id, x, y, r, pa, a, a, a, dv, vel, format='d,d,d,d,d,d,d,d,d,d', /silent
    d = arr_struct({ra:x/3600d + self.galra, dec:y/3600d + self.galdec, vel:vel-1261.2, errvel:dv});might need vel:vel-1241. (NED Heliocentric Radial Velocity)
    d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot))) 
    d.x = -d.x
    d=struct_trimtags(d,except=['PA'])
    d = struct_addtags(d, arr_struct({pa:pacalc(d.x, d.y, self.parot,/radians)}))  
    return, d

;    readcol, 'I:/Discrete Tracers/NGC4697/Mendez.data1_newformat', id, ra, dec, vel, format='d,d,d,d', /silent
;    d = arr_struct({ra:ra, dec:dec, vel:vel-1241.});, errvel:dv});might need vel:vel-1241. (NED Heliocentric Radial Velocity)
;    d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot))) 
;
;  return, d
  
end
;-----------------------------------------------------------------------------------------------

FUNCTION map_ngc4697::get_sauron, nosyscorr=nosyscorr, sum=sum, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
	d = self->map_galaxy::get_sauron($
;		pxf = '/Users/jaaarnol/sauron/cappellari/PXF_bin_MS_NGC4697_r1_idl.fits', $
    pxf = 'I:/Discrete Tracers/NGC4697/PXF_bin_MS_NGC4697_r1_idl.fits',vs_fitting);, medianing=medianing);,$
;		nosyscorr=nosyscorr, minVel=minVel, mxVel=mxVel, maxVeldisp=maxVeldisp, sum=sum)		
  d.x = -d.x
  d=struct_trimtags(d,except=['PA'])
  d = struct_addtags(d, arr_struct({pa:pacalc(d.x, d.y, self.parot,/radians)}))
	return, d
END

;-----------------------------------------------------------------------------------------------

PRO map_ngc4697::get_data, d, nopne=nopne, nogcs=nogcs, nodiscrete=nodiscrete, nostellar=nostellar, $
	nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, red=red, blue=blue, h3=h3, h4=h4, $
	veldisp=veldisp, mdist=mdist, medianing=medianing,vs_fitting=vs_fitting
	compile_opt idl2, hidden
	if ((n_elements(nopne) eq 1) and (n_elements(nogcs) eq 1)) then nodiscrete = 1
	if (n_elements(nodiscrete) eq 0) then begin
		if (n_elements(nopne) eq 0) then discrete = struct_append(discrete, self->get_pne())
		if (n_elements(nogcs) eq 0) then discrete = struct_append(discrete, self->get_GCs(red=red, blue=blue))
		*(self.discrete) = discrete									;put the discrete data in a pointer
		d = discrete
	endif

    if 0 then begin
        d = self->get_smeag(medianing=medianing)
        d1 = self->get_deimos_longslit()
        d2 = self->get_pinkney_longslit()
        stop
    endif
	
	if ((n_elements(nosauron) eq 1) and (n_elements(nolong) eq 1) and (n_elements(nosmeag) eq 1)) then nostellar = 1
	if (n_elements(nostellar) eq 0) then begin
;		self->srn_corr, mdist=mdist, sum=sum
		if (n_elements(nosauron) eq 0) then stellar = self->get_sauron(sum=sum,vs_fitting);, medianing=medianing)
		if (n_elements(nolong) eq 0) then begin
			ls = struct_append(ls, self->get_deimos_longslit())
;			if (n_elements(ls) gt 0) then stellar = struct_append(stellar, ls) ;not correct LS file
		endif
		if (n_elements(nosmeag) eq 0) then stellar = struct_append(stellar, self->get_smeag(vs_fitting));medianing=medianing))
		if (n_elements(stellar) gt 0) then *(self.stellar) = stellar
		s = stellar
	endif
    if n_elements(nosauron) eq 0 or n_elements(nosmeag) eq 0 or n_elements(nolong) eq 0 or n_elements(nostellar) eq 0 then d = struct_selecttags(s,select=['x','y','r','pa','vel','errvel']) ;subarray to hold stellar+discrete data (start by adding stellar data)
    if n_elements(nodiscrete) eq 0 then begin
      dummy = struct_selecttags(discrete,select=['x','y','r','pa','vel','errvel']) 
      d = struct_append(d,dummy)
    endif
END

;-----------------------------------------------------------------------------------------------
PRO map_ngc4697::rolling_bins_kinemetry, d, nodiscrete=nodiscrete, nostellar=nostellar, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
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

pro map_ngc4697::isophote, p2
  compile_opt idl2, hidden

  s = *(self.plot)
  pos = self.pos

nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
q = 0.70 ; outer isophote flattening
;maxa = 1d/sqrt(q)*s.xsiz/float(difference(minmax(s.xrange)))
;maxa = 2*self.reff/sqrt(1-q)-25d
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
q = 0.64  ;inner isophote flattening
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
PRO map_ngc4697__define
	compile_opt idl2, hidden

	void={map_ngc4697, inherits map_galaxy}

END






 
 
