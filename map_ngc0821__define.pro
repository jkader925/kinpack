
;-----------------------------------------------------------------------------------------------

pro map_ngc0821::get_variogram_info, vs_fitting
	compile_opt idl2, hidden
	d = self->get_smeag2(vs_fitting)
	d = struct_trimtags(d, select=['X','Y','VEL','VELDISP','H3','H4','RMS','FLUX'])
		dummy = execute('d = struct_addtags(d, arr_struct({z:d.'+self.kinType+'}))')	

	case self.kinType of 
		'VEL'		: 	BEGIN
							svar_info = {range:150., ftype:'POLY', axisratio:0.677d, inc:3}
						END
		'VELDISP'	: 	BEGIN
							svar_info = {range:150., ftype:'LINEAR', axisratio:1.0d, inc:10}
						END
		'RMS'		: 	BEGIN
							svar_info = {range:150., ftype:'LINEAR', axisratio:1.0d, inc:10}
						END
		'H3'		: 	BEGIN
							svar_info = {range:120., ftype:'POLY', axisratio:1.0d, inc:4}
						END
		'H4'		: 	BEGIN
							svar_info = {range:120., ftype:'POLY', axisratio:1.0d, inc:4}
						END
		'FLUX'		: 	BEGIN
							svar_info = {range:60., ftype:'POLY', axisratio:1.0d, inc:4}
						END
		else		: stop
	endcase	
	if (n_elements(self.maxdist) eq 0) then self.maxDist = svar_info.range*0.5
	svar_info = struct_addtags(svar_info, {maxDist:self.maxDist})
	svar_info = struct_addtags(svar_info, {data:d})
	*(self.svar_info) = svar_info
end

;-----------------------------------------------------------------------------------------------

;pro map_ngc0821::get_gal_info
;	compile_opt idl2, hidden
;	self.name = 'NGC821'
;	self.reff = 40. ;sluggssurvey paper
;	self.distance = 23.4d6 ;sluggssurvey paper 
;	self.q = 0.65  ;sluggssurvey paper
;	self.parot = 31.2d ;sluggssurvey paper
;	self.galra = 32.088083d  
;	self.galdec = 10.994917d
;	self.arc2kpc = self.distance/206265d/1000d
;	self->get_variogram_info
;
;end


pro map_ngc0821::get_gal_info,vs_fitting
  compile_opt idl2, hidden
  self.name = 'NGC821'
  self.reff = 40.
  self.distance = 24.4d6
  self.q = 0.677
  self.parot = 30.5d;90d;32.2d
  self.galra = 32.088083d  
  self.galdec = 10.994917d
  self.arc2kpc = self.distance/206265d/1000d
  self->get_variogram_info,vs_fitting
end

;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_weijmans_sauron
	compile_opt idl2, hidden
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/weijmans_sauron_NGC821.dat', name, dx, dy, vel, errvel, veldisp, errveldisp, h3, errh3, h4, errh4, vmean, $
		errvmean, veldispmean, errveldispmean, format='a,d,d,d,d,d,d,d,d,d,d,d,d,d,d', /silent
	print, 'Need to know what PA Weijmans used for the minor-axis pointing.'
	relPA = [25d,32d,32d,32d]

;	stop
	xsiz = 43d
	ysiz = 35d
	tmp = arr_struct({$
		x:(cmreplicate((findgen(xsiz)/(xsiz - 1d) - 0.5d)*41d, ysiz))[indgen(xsiz*ysiz)], $
		y:(cmreplicate((findgen(ysiz)/(ysiz - 1d) - 0.5d)*33d, xsiz))[indgen(xsiz*ysiz)]})
	
	tmp = replicate({x:0d, y:0d},4)
	tmp = {x:0d, y:0d}
	
	
	for i=0,3 do begin
		d = arr_struct({dx:tmp.x, dy:tmp.y})
		if (i eq 0) then ang = 0d else ang = 90d
		xyrot, d.dx, d.dy, ang*!dtor, xp, yp
		;xp += sign(dx[i])*sqrt(dx[i]^2d + dy[i]^2d)
		xyrot, xp, yp, (90d - relPA[i])*!dtor, xp, yp
		xp += dx[i]
		yp += dy[i]
		d.dx = xp
		d.dy = yp
		foreach tag, ['VEL','ERRVEL','VELDISP','ERRVELDISP','H3','ERRH3','H4','ERRH4'] do begin
			if (n_elements(tmp) gt 1) then dummy = execute('d = struct_addtags(d, arr_struct({'+tag+':tmp.x*0d + '+tag+'[i]}))') $
				else dummy = execute('d = struct_addtags(d, {'+tag+':tmp.x*0d + '+tag+'[i]})')
		endforeach
		d.errvel *= sqrt(n_elements(d))
		hold = struct_append(hold, d)
	endfor
	d = hold
	srn = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
	srn = struct_addtags(srn, struct_trimtags(d, except=['DX','DY']))
	srn = struct_addtags(srn, arr_struct(projdist2(srn.ra, srn.dec, self.galra, self.galdec, self.parot)))
	self->calc_rms, srn
;	self->get_Flux, srn
	return, srn
end

;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_pinkney_mdm
	compile_opt idl2, hidden
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/pinkney_mdm_NGC821.dat', name, pa, sym, r, vel, errvel, veldisp, errveldisp, h3, errh3, h4, errh4, format='a,d,a,d,d,d,d,d,d,d,d,d', /silent
	print, 'Need to know what major-axis PA Pinkney used.'
	relPA = 25d
	d = arr_struct({name:name, pa:pa, r:r, dx:r*0d, dy:r*0d, vel:vel, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp, h3:h3, errh3:errh3, h4:h4, errh4:errh4})
	d.dx = d.r*cos((d.pa + 90d - relPA)*!dtor)
	d.dy = d.r*sin((d.pa + 90d - relPA)*!dtor)
	ind = where(d.pa ne 0)
	d1 = d[where(d.pa ne 0)]
	d1.dx = -d1.dx
	d1.dy = -d1.dy
	d = struct_append(d, d1)			
	mdm = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
	mdm = struct_addtags(mdm, struct_trimtags(d, except=['PA','R','DX','DY','SYM']))
	mdm = struct_addtags(mdm, arr_struct(projdist2(mdm.ra, mdm.dec, self.galra, self.galdec, 29.0)))
	self->calc_rms, mdm	
;	self->get_Flux, mdm
	return, mdm
end

;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_proctor_gmos
	compile_opt idl2, hidden
;	stop
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/proctor_gmos_NGC821.dat', pa, r, vel, errvel, veldisp, errveldisp, format='d,d,d,d,d,d', /silent
	d = arr_struct({pa:pa, r:r, dx:r*0d, dy:r*0d, vel:vel, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp})
	d.dx = d.r*cos((90d - d.pa)*!dtor)
	d.dy = d.r*sin((90d - d.pa)*!dtor)
	gmos = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
	gmos = struct_addtags(gmos, struct_trimtags(d, except=['PA','R','DX','DY','SYM']))
	gmos = struct_addtags(gmos, arr_struct(projdist2(gmos.ra, gmos.dec, self.galra, self.galdec, self.parot)))
	gmos.vel -= 1714.
	self->calc_rms, gmos
;	self->get_Flux, gmos	
	return, gmos
end

;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_forestell_het
	compile_opt idl2, hidden
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/forestell_HET_major_NGC821.dat', r, vel, errvel, veldisp, errveldisp, h3, errh3, h4, errh4, format='d,d,d,d,d,d,d,d,d', /silent
	print, 'Need to know what major-axis PA Forestell used.'
	relPA = 32d
	d = arr_struct({pa:r*0d + relPA, r:r, dx:r*0d, dy:r*0d, vel:vel, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp, h3:h3, errh3:errh3, h4:h4, errh4:errh4})
	d.dx = d.r*cos((90d - d.pa)*!dtor)
	d.dy = d.r*sin((90d - d.pa)*!dtor)
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/forestell_HET_minor_NGC821.dat', r, vel, errvel, veldisp, errveldisp, h3, errh3, h4, errh4, format='d,d,d,d,d,d,d,d,d', /silent
	d1 = arr_struct({pa:r*0d + relPA + 90d, r:r, dx:r*0d, dy:r*0d, vel:vel, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp, h3:h3, errh3:errh3, h4:h4, errh4:errh4})
	d1.dx = d1.r*cos((90d - d1.pa)*!dtor)
	d1.dy = d1.r*sin((90d - d1.pa)*!dtor)
	d2 = d1
	d2.dx *= -1
	d2.dy *= -1
	d1 = struct_append(d1, d2)
	d = struct_append(d, d1)
	het = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
	het = struct_addtags(het, struct_trimtags(d, except=['PA','R','DX','DY']))
	het = struct_addtags(het, arr_struct(projdist2(het.ra, het.dec, self.galra, self.galdec, self.parot)))
	self->calc_rms, het
;	self->get_Flux, het	
	return, het
end

;-----------------------------------------------------------------------------------------------
function map_ngc0821::get_discrete_combine
   compile_opt idl2, hidden
   
   readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/N821.PNGCdata1', id, x, y, r, pa, a, a, vel, dv, a, format='d,d,d,d,d,d,d,d,d,d', /silent
;   readcol, 'E:/Discrete Tracers/NGC821/N821_v2.GCdata3c', id, x, y, r, pa, a, a, vel, dv, a, format='a,d,d,d,d,d,d,d,d,d', /silent
   d = arr_struct({id:id, ra:x/3600d + self.galra, dec:y/3600 + self.galdec, vel:vel-1696.0, errvel:dv})
   d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))  
  return, d
end
;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_PNe
	compile_opt idl2, hidden
	
	if 0 then begin
		readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/N821.PNdata1a', id, x, y, r, pa, a, a, vel, dv, a, format='d,d,d,d,d,d,d,d,d,d', /silent
		d = arr_struct({ra:x/3600d + self.galra, dec:y/3600d + self.galdec, vel:vel - 1696d, errvel:dv})
	endif else begin
		readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/N821.PNdata2b', id, x, y, r, pa, a, a, vel, dv, a, format='d,d,d,d,d,d,d,d,d,d', /silent
;		d = arr_struct({ra:x/3600d + self.galra, dec:y/3600d + self.galdec, vel:vel - 1696., errvel:dv})
    d = arr_struct({ra:x/3600d + self.galra, dec:y/3600d + self.galdec, vel:vel - 1696d, errvel:dv})
	endelse
;stop
	;d = arr_struct({xi:x, yi:y, vel:vel - 1700.7})					;x and y seem to just be delta(RA/DEC) from the center
	;xyrot, d.xi, d.yi, (0*32d - 90d)*!dtor, xp, yp
	;d = struct_addtags(d, arr_struct({x:xp, y:yp}))
	;plot, x, y, psym=1, xrange=max(abs(minmax(x)))*1.2*[1,-1], /xstyle, /iso
	d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))	
	;plot, d.x, d.y, psym=1, /iso
	;plot, d.ra, d.dec, xrange=reverse(minmax(d.ra)), /xstyle, /iso, yrange=minmax(d.dec), psym=1

	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_GCs, red=red, blue=blue, green=green, medianing=medianing
	compile_opt idl2, hidden	
	readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/NGC821_GCs_pota.dat', id, ra, dec, vel, errvel, g, errg, r, errr, i, erri, v, errv, II, errII, format='a,d,d,d,d,d,d,d,d,d,d,d,d,d,d', /silent
;  d = arr_struct({ra:ra, dec:dec, vel:vel, errvel:errvel})
;	readcol, 'E:/Discrete Tracers/NGC821/Kin_NGC821_J_2.dat', id, ra, dec, vel, errvel, g, errg, r, errr, i, erri, type, format='a,d,d,d,d,d,d,d,d,d,d,i', /silent

;	if n_elements(medianing) eq 1 then begin
;	d = arr_struct({ra:ra, dec:dec, vel:vel - median(vel,/even), errvel:errvel})
;	endif else begin
	d = arr_struct({ra:ra, dec:dec, vel:vel - 1735d, errvel:errvel}) ;NEDvalue=1735km/s, sluggs=1718km/s, twozonepaper=1696km/s
;  endelse
;	d = struct_append(d,d1)
	d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))	
	;plot, d.x, d.vel, psym=1
 if (n_elements(red) gt 0) then d = d[where((g - i) gt 1.01)] else if (n_elements(blue) gt 0) then d = d[where((g - i) le 0.78)] $
  else if (n_elements(green) gt 0) then d = d[where((g-i) gt 0.78 and (g-i) lt 1.01)] ;for pota GC data
;	if (n_elements(red) gt 0) then d = d[where(type eq 1)] else if (n_elements(blue) gt 0) then d = d[where(type eq 0)] ; for other GC data
	d = d[where(d.vel gt -400)]
;	self->get_Flux, d	

	return, d
end

;-----------------------------------------------------------------------------------------------
function map_ngc0821::get_smeag2, vs_fitting
  compile_opt idl2, hidden
  
 medvel = 1712d
 file = 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/NGC821_SKiMS.dat5c'
   readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/NGC821_SKiMS.dat5c', ID, x, y, R, PA, R_mkin, R_a, vel, dvel, sig, dsig, h_3, dh_3, h_4, dh_4, R_mphot, format='a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', /silent
data = arr_struct({ra:x/3600d + self.galra, dec:y/3600 + self.galdec, vel:vel, errvel:dvel})
data = struct_addtags(data, arr_struct(projdist2(data.ra, data.dec, self.galra, self.galdec, self.parot)))
;
;if vs_fitting eq 1 then begin
;  ;----------------find v_sys using kinemetry
; inner = 1
; outer = 0
; sauron = 0
; smeagol = 1
; reff = self.reff
; name = self.name
; TwoBin_Kinemetry, data, reff, vsys, diagnostics, name, inner, outer, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, dvdisp, dinnerpa,  $
;    nopne=nopne, nogcs=nogcs, vs_fitting, vs_fit, dvs_fit, sauron, smeagol
; undefine,data
;endif
;;-------------------------------- 

; if vs_fitting eq 1 then offsetVel = vs_fit
; if vs_fitting eq 0 then offsetVel = medvel
offsetVel = medvel

;  if vs_fitting eq 1 then print,strcompress('SMEAGOL fitted v_sys ='+string(vs_fit)+'+/-'+string(dvs_fit))
;  if vs_fitting eq 0 then print,strcompress('SMEAGOL photometric v_sys='+string(medvel))
  
   readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/NGC821_SKiMS.dat5c', ID, x, y, R, PA, R_mkin, R_a, vel, dvel, sig, dsig, h_3, dh_3, h_4, dh_4, R_mphot, format='a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', /silent
   d = arr_struct({ra:x/3600d + self.galra, dec:y/3600 + self.galdec, vel:vel-offsetVel, errvel:dvel})
   d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))
  return, d
  
end

;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_smeag, medianing=medianing, vs_fitting
	compile_opt idl2, hidden
;	if self.vor then file = '/Volumes/data4/skimspaper/NGC0821/comb_deimos_ppxf_v0.fits' $
; if self.vor then file = 'E:/Discrete Tracers/NGC821/comb_deimos_ppxf_v0.fits' $
;		else file = 'E:/Discrete Tracers/NGC821/comb_deimos_ppxf_v0.fits'
;	d = self->get_skims(file)
file = 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/comb_deimos_ppxf_NGC821.fits'
; file = 'E:/Discrete Tracers/NGC821/deimos_ppxf_NGC821.fits'

d = self->get_skims(file, vs_fitting);,medianing=medianing)

	if 0 then begin
		plot, d.pa, d.h3, psym=1
		
		c = self->get_skims('I:/Brutus_InternalDrive/Discrete Tracers/NGC821/comb_deimos_ppxf_v0.fits')
		d = self->get_skims('I:/Brutus_InternalDrive/Discrete Tracers/NGC821/comb_deimos_ppxf_v1.fits')
		f = self->get_skims('I:/Brutus_InternalDrive/Discrete Tracers/NGC821/deimos_ppxf_n821_1_1.fits')
		d = d[sort(d.pa)]
	
		ploterror, c.pa, c.vel, c.errvel, psym=1, yrange=[-1,1]*150
		oplot, d.pa, d.vel, psym=-6, color=fsc_color('red')
		
		ploterror, c.y, c.veldisp, c.errveldisp, psym=1, yrange=[0,1]*350, /nohat
		oploterror, d.y, d.veldisp, d.errveldisp, psym=1, errcolor=fsc_color('red'), /nohat, color=fsc_color('red')
		
		ploterror, c.pa, c.h3, c.errh3, psym=1, yrange=[-1,1]*0.5
		oplot, d.pa, d.h3, psym=-6, color=fsc_color('red')

		
		ploterror, c.pa, c.h4, c.errh4, psym=1, yrange=[-1,1]*0.5
		oplot, d.pa, d.h4, psym=-6, color=fsc_color('red')

	
	endif
	;d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))
	;smg = smg[where( (smg.sn gt 4.5) )]
	;smg = smg[where( (smg.r lt 250) )]
	return, d
end

;-----------------------------------------------------------------------------------------------
PRO map_ngc0821::rolling_bins_kinemetry, d, nodiscrete=nodiscrete, nostellar=nostellar, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
    nopne=nopne, nogcs=nogcs, nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, N0821, N1023, N1344, N2768, N3115, N3377, N4564, N4697, N4473, inplotdat, inmaxr, bin1, l, nbootstraps, twobin=twobin,$
    obser, peakplot=peakplot, qin, dqin, nbootstrapsq
compile_opt idl2, hidden
   
    diagnostics=1
    data = d
    reff = self.reff
    vsys = self.vsys
    name = self.name
    vs_fitting = vs_fitting
    smeagol = 0 
    sauron = 0

if n_elements(nostellar) eq 1 then outer = 1
if n_elements(nostellar) eq 1 then inner = 0
if n_elements(nodiscrete) eq 1 then outer = 0
if n_elements(nodiscrete) eq 1 then inner = 1

 TwoBin_Kinemetry, data, reff, vsys, diagnostics, name, inner, outer, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
    nopne=nopne, nogcs=nogcs, vs_fit, dvs_fit, sauron, smeagol, N0821, N1023, N1344, N2768, N3115, N3377, N4473, N4564, N4697, inplotdat, inmaxr, bin1, l, nbootstraps, twobin=twobin,$
    obser, peakplot=peakplot, qin, dqin, nbootstrapsq

END

;-----------------------------------------------------------------------------------------------

function map_ngc0821::get_sauron, medianing=medianing, vs_fitting
	compile_opt idl2, hidden
;	stop
;	d = self->map_galaxy::get_sauron(pxf = '/Users/jaaarnol/sauron/PXF_NGC821_r2/PXF_bin_MS_NGC821_r2.fits', $
  d = self->map_galaxy::get_sauron(pxf = 'I:/Brutus_InternalDrive/Discrete Tracers/NGC821/PXF_bin_MS_NGC0821_r2_idl.fits',$;PXF_bin_MS_NGC821_r2.fits', $
									ms='I:/Brutus_InternalDrive/Discrete Tracers/NGC821/MS_NGC821_r2.fits', vs_fitting);, medianing=medianing)	
;  d = self->map_galaxy::get_sauron(pxf = 'E:/Discrete Tracers/NGC821/PXF_bin_MS_NGC0821_r2_idl.fits', $
;                  ms='E:/Discrete Tracers/NGC821/MS_NGC821_r2.fits', medianing=medianing) 



	return, d
end

;-----------------------------------------------------------------------------------------------

pro map_ngc0821::get_data, d, nopne=nopne, nogcs=nogcs, nodiscrete=nodiscrete, nostellar=nostellar, $
	nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, red=red, blue=blue, green=green, medianing=medianing, $
	vs_fitting=vs_fitting
	
	compile_opt idl2, hidden
	if ((n_elements(nopne) eq 1) and (n_elements(nogcs) eq 1)) then nodiscrete = 1
	if (n_elements(nodiscrete) eq 0) then begin
	  discrete = self->get_discrete_combine()

;		if (n_elements(nopne) eq 0) then discrete = self->get_PNe()								;retrieve the PN data
;		if (n_elements(nogcs) eq 0) then discrete = struct_append(discrete, self->get_GCs(red=red, blue=blue))	;retrieve the GC data

		;print, 'Need errvels for the NGC 821 GC and PN data!'
		;discrete = struct_addtags(discrete, arr_struct({errvel:discrete.x*0d + 25d}))
		*(self.discrete) = discrete																;put the discrete data in a pointer
		dis = discrete
	endif

	if ((n_elements(nosauron) eq 1) and (n_elements(nolong) eq 1) and (n_elements(nosmeag) eq 1)) then nostellar = 1
	if (n_elements(nostellar) eq 0) then begin
		if (n_elements(nosauron) eq 0) then stellar = self->get_sauron(vs_fitting);medianing=medianing)
		if (n_elements(nolong) eq 0) then begin
			ls = struct_append(ls, 	self->get_pinkney_mdm())
			ls = struct_append(ls, self->get_forestell_het())
			ls = struct_append(ls, self->get_weijmans_sauron())					
			if ((self.kinType ne 'H3') and (self.kinType ne 'H4')) then begin
				ls = struct_append(ls, self->get_proctor_gmos())			
			endif
			stellar = struct_append(stellar, ls)
		endif

		if (n_elements(nosmeag) eq 0) then stellar = struct_append(stellar, self->get_smeag2(vs_fitting));medianing=medianing))
		*(self.stellar) = stellar
    s = stellar
			svar_info = *(self.svar_info)
    endif    

    if n_elements(nosauron) eq 0 or n_elements(nosmeag) eq 0 or n_elements(nolong) eq 0 or n_elements(nostellar) eq 0 then d = struct_selecttags(s,select=['x','y','r','pa','vel','errvel']) ;subarray to hold stellar+discrete data (start by adding stellar data)
    if n_elements(nodiscrete) eq 0 then begin
      dummy = struct_selecttags(dis,select=['x','y','r','pa','vel','errvel']) 
      d = struct_append(d,dummy)
    endif
end

;--------------------------------------------------------------------------------------------------------------------    
; Make Sigma vs. R/R_e plots to investigate "sigma bump" e.g. Schauer, A.T.P. et al. 2013
;    stellar = struct_append(stellar, discrete)  ;test to see about splicing together discrete/stellar data
;    
;    stellar.R = stellar.R/39.
;    pl = plot(stellar.R, stellar.VELDISP,linestyle=6,symbol="+",xrange=[0,8])
;--------------------------------------------------------------------------------------------------------------------    


	
	      ;Scale to Vmax
;      discrete.vel = discrete.vel/max(stellar.vel)
;      stellar.vel = stellar.vel/max(stellar.vel)
;      *(self.discrete) = discrete 
;      *(self.stellar) = stellar
;      print,self.name
;      print,'maximum stellar velocity =',max(stellar.vel)
;      print,'Effective Radius [arcsec] =',self.reff


;-----------------------------------------------------------------------------------------------

pro map_ngc0821::isophote, p2
	compile_opt idl2, hidden

	s = *(self.plot)
	pos = self.pos
;s.xsiz = 500;357.;500;
;s.ysiz = 400;285.7
;s.xrange = [-7,7]
;s.yrange = [-4,4]
;s.xrange = [-5,5]
;-----------------------------------------------------------------------------------------------
nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
;q = 0.68;0.8 ; outer isophote flattening
q = 0.733 ;joel outer photometry
maxa = (4/7d)*(max(s.xsiz)/2)
;maxa = 1d/sqrt(q)*s.xsiz/float(difference(minmax(s.xrange)))
;maxa = 4*(s.xsiz/float(difference(minmax(s.xrange))))
xpa = -1.5d*!dtor ;outer isophote position angle
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
q = 0.677  ;inner isophote flattening sluggssurvey paper
;maxa = 1d/sqrt(q)*s.xsiz/float(difference(minmax(s.xrange)))

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
  outeryhigh = yhigh/isophote2 + (ymid - ymid/isophote2)  
  
  xhigh2 = [outerxhigh, outerxhigh[0]]
  yhigh2 = [outeryhigh, outeryhigh[0]]

  p2 = plot(xhigh1, yhigh1, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)
  p2 = plot(xhigh2, yhigh2, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)

  
;  p2 = plot(xhigh2, yhigh2, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)


;  files = 'E:/Discrete Tracers/NGC3377/n3377_'+roundx(s.isophote_level,0)+'Re.dat'

;	k = 0
;	foreach confile, files do begin
;	 stop
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
;;		con.x = smooth(con.x*s.xsiz/float(difference(minmax(s.xrange))) + self.xmid, smcon, /edge)
;;		con.y = smooth(con.y*s.xsiz/float(difference(minmax(s.xrange))) + self.ymid, smcon, /edge)
;    con.x = smooth(con.x*s.xsiz/float(difference(minmax(s.xrange))) + self.xmid, smcon)
;    con.y = smooth(con.y*s.xsiz/float(difference(minmax(s.xrange))) + self.ymid, smcon)
;		con = struct_append(con, con[0])
;		stop
;		p2 = plot(con.x, con.y, /overplot, /current);, thick=s.c_thick)
;		++k
;	endforeach
;xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
;ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
;
;maxa = self.reff/sqrt(self.Q)*self.arc2kpc*s.xsiz/float(difference(minmax(s.xrange)))
;nbound = 100
;xpa = 0d*!dtor

	
end

;-----------------------------------------------------------------------------------------------

pro map_ngc0821__define
	compile_opt idl2, hidden

	void={map_ngc0821, inherits map_galaxy}

end






 
 
