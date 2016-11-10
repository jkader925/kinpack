

function map_galaxy::init, s, doreff=doreff, norm=norm, gsm=gsm, veldisp=veldisp, h3=h3, $
	h4=h4, maxDist=maxDist, vor=vor, rms=rms, dss=dss
	compile_opt idl2, hidden
	if (n_elements(doreff) gt 0) then self.doreff = 1
	if (n_elements(norm) eq 0) then self.norm = 1d else self.norm = norm
	self.gsm = ptr_new(/allocate_heap)
	if (n_elements(gsm) gt 0) then *(self.gsm) = gsm
	self.discrete = ptr_new(/allocate_heap)
	self.stellar = ptr_new(/allocate_heap)
	self.obj = ptr_new(/allocate_heap)
	self.plot = ptr_new(/allocate_heap)
	if (n_elements(s) gt 0) then *(self.plot) = s
	self.image = ptr_new(/allocate_heap)
	self.current = ptr_new(/allocate_heap)
	self.dss = ptr_new(/allocate_heap)
	self.svar_info = ptr_new(/allocate_heap)
	self.kinType = 'VEL'
	if (n_elements(dss) gt 0) then self.kinType = 'FLUX'
	if (n_elements(rms) gt 0) then self.kinType = 'RMS'
	if (n_elements(veldisp) gt 0) then self.kinType = 'VELDISP'
	if (n_elements(h3) gt 0) then self.kinType = 'H3'
	if (n_elements(h4) gt 0) then self.kinType = 'H4'
	if (n_elements(maxdist) gt 0) then self.maxdist = maxDist
	if (n_elements(vor) gt 0) then self.vor = 1 else self.vor = 0
	return, 1
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::cleanup
	compile_opt idl2, hidden
	ptr_free, self.discrete
	ptr_free, self.stellar
	ptr_free, self.obj
	ptr_free, self.plot
	ptr_free, self.image
	ptr_free, self.current
	ptr_free, self.dss
	ptr_free, self.svar_info
	ptr_free, self.gsm
end

;-----------------------------------------------------------------------------------------------

function map_galaxy::get, coplist
	compile_opt idl2, hidden
	all = create_struct(name=obj_class(self))
	struct_assign, self, all
	if (n_elements(coplist) gt 0) then all = struct_trimtags(all, select=coplist)
	return, all
end

;-----------------------------------------------------------------------------------------------

function map_galaxy::get_sauron, pxf=pxf, ms=ms, pz=pz, medianing=medianing, vs_fitting
  compile_opt idl2, hidden
  if (n_elements(pxf) eq 0) then stop
  if self.name NE 'NGC4697' then if (n_elements(ms) eq 0) then stop  ;comment out for 4697
  epxf = file_search(pxf,count=nfile)
  if (nfile eq 0) then stop
  if self.name NE 'NGC4697' then ems = file_search(ms,count=nfile) ;comment out for 4697
  if (nfile eq 0) then stop
  dpxf = mrdfits(pxf,1,/silent)
  if self.name NE 'NGC4697' then dms = mrdfits(ms,1,/silent) ;comment out for 4697
  

  if self.name NE 'NGC4697' then match = match_2d(dms.a, dms.d, dpxf.xs, dpxf.ys, 20d, match_distance=mindist)    ;find closest bin to each pixel  ;comment out for 4697
  if self.name NE 'NGC4697' then new_dpxf = dpxf[match]                               ;match up those bins to the set of pixels ;comment out for 4697

if (n_elements(pz) gt 0) then stop

; new_dpxf = struct_trimtags(new_dpxf, except=['VPXF','SPXF'])  ;this was already commented
  if self.name NE 'NGC4697' then new_dms = struct_trimtags(dms, except=['NO','FLUX','VPXF','SPXF'])  ;comment out for 4697
  new_dpxf = struct_trimtags(new_dpxf, except=['A','D'])
  
  if self.name NE 'NGC4697' then new = struct_addtags(new_dms, new_dpxf) ;comment out for 4697

  if self.name EQ 'NGC4697' then new = dpxf  ;comment out if NOT 4697
  srn = new
  srn = struct_addtags(srn, arr_struct({ra:self.galra - srn.a/3600d, dec:self.galdec + srn.d/3600d}))
  
  srn = dpxf
  srn = struct_addtags(srn, arr_struct({ra:self.galra - srn.xs/3600d, dec:self.galdec + srn.ys/3600d}))
  
  
  srn = struct_addtags(srn, arr_struct({vel:srn.vpxf, errvel:srn.evpxf, veldisp:srn.SPXF, $
    errveldisp:srn.ESPXF, h3:srn.H3PXF, errh3:srn.eH3PXF, h4:srn.H4PXF, errh4:srn.eH4PXF}))

 if n_elements(medianing) eq 1 then begin
  
 medVel = median(smg.vel,/even) ;calculate median velocity
 srn.vel -= medVel  ; offset velocity
  
 endif else begin

 if self.name eq 'NGC821' then offsetVel =   1707.9
 if self.name eq 'NGC1023' then offsetVel =  604.9
 if self.name eq 'NGC2768' then offsetVel =  1379.2
 if self.name eq 'NGC3115' then offsetVel =  644.8
 if self.name eq 'NGC3377' then offsetVel =  680.0
 if self.name eq 'NGC4473' then offsetVel =  2244.0
 if self.name eq 'ngc4564' then offsetVel = 1155.0 ;SLUGGS Survey Paper Brodie+2014
 if self.name eq 'NGC4697' then offsetVel =  1265.0
 
medvel = offsetVel
 
 endelse
 
  print,self.name
;  print,'SAURON kinemetry v_sys =',offsetVel
  
  srn = struct_addtags(srn, arr_struct(projdist2(srn.ra, srn.dec, self.galra, self.galdec, self.parot)))
  srn = struct_trimtags(srn, except='flux')
;if vs_fitting eq 1 then begin
;;----------------find v_sys using kinemetry
; inner = 1
; outer = 0
; sauron = 1
; smeagol = 0
;; vs_fitting = 1
; data = srn
; reff = self.reff
; name = self.name
; TwoBin_Kinemetry, data, reff, vsys, diagnostics, name, inner, outer, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, dvdisp, dinnerpa,  $
;    nopne=nopne, nogcs=nogcs, vs_fitting, vs_fit, dvs_fit, sauron, smeagol
; undefine,data
;endif
;-------------------------------- 
; if vs_fitting eq 1 then offsetVel = vs_fit
; if vs_fitting eq 0 then offsetVel = medvel
offsetVel = medvel

 srn.vel -= offsetVel
; if vs_fitting eq 1 then print,strcompress('SAURON fitted v_sys ='+string(vs_fit)+'+/-'+string(dvs_fit))
; if vs_fitting eq 0 then print,strcompress('SAURON photometric v_sys='+string(medvel))


;srn.x = -srn.x ;everyone else is using +x -> E, but SAURON used +x -> W

; if self.name eq 'NGC3377' or self.name eq 'NGC4697' then srn.x = -srn.x
;--------------------------------------------------------------------------------------kinemetry
  


;--------------------------------------------------------------------------------------

  self->calc_rms, srn
;  self->get_Flux, srn
  
; srn.vel += medVel
  return, srn

  efile = file_search(file,count=nfile)
  if (nfile eq 0) then stop
  srn = mrdfits(file,1,/silent)
  srn = struct_addtags(srn, arr_struct({ra:self.galra - srn.xs/3600d, dec:self.galdec + srn.ys/3600d, $
  vel:srn.vpxf, errvel:srn.evpxf, veldisp:srn.SPXF, errveldisp:srn.ESPXF, h3:srn.H3PXF, errh3:srn.eH3PXF, $
  h4:srn.H4PXF, errh4:srn.eH4PXF}))
; srn.vel -= median(srn.vel) ;should this be here?
  srn = struct_addtags(srn, arr_struct(projdist2(srn.ra, srn.dec, self.galra, self.galdec, self.parot)))
  return, srn
end
;-----------------------------------------------------------------------------------------------

pro map_galaxy::get_flux, d
	compile_opt idl2, hidden
	if ~tag_exist(d, 'FLUX') then begin
		if (n_elements(*(self.dss)) eq 0) then self->get_dss
		dss = *(self.dss)
		match = match_2d(d.ra, d.dec, dss.ra, dss.dec, 5d/3600d, match_distance=mindist)
		ibad = where(match eq -1,nbad)
		if (nbad gt 0) then stop
		d = struct_addtags(d, arr_struct({flux:double(dss[match].val), errflux:double(sqrt(dss[match].val))}))
	endif
end

;-----------------------------------------------------------------------------------------------

function map_galaxy::get_skims, file, vs_fitting
	compile_opt idl2, hidden
	
	
	
	 
;  if self.name NE 'ngc4564' then begin
  efile = file_search(file,count=nfile)
  if (nfile eq 0) then stop
  smg = mrdfits(file,1,/silent) 
;  endif
;
;  if self.name eq 'ngc4564' then begin
;  readcol,'E:/Discrete Tracers/NGC4564/skims_data_4564.txt',ra,dec,vel,errvel,veldisp,errveldisp,flag,$
;  format='d,d,d,d,d,d,a',/silent
;  smg = arr_struct({ra:ra, dec:dec, vel:vel, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp, flag:flag})
;  smg = smg[where(smg.flag NE 'D')]
;  endif

	medVel = median(smg.vel,/even)
	if self.name eq 'NGC0821' then medVel = 1712d
	if self.name eq 'NGC1023' then medVel = 604.9
	if self.name eq 'NGC2768' then medVel = 1379.2
	if self.name eq 'NGC3115' then medVel = 644.8
	if self.name eq 'NGC3377' then medVel = 680.8
	if self.name eq 'ngc4564' then medVel = 1155.0;1309.5144;1155.0
	if self.name eq 'ngc4697' then medVel = 1265.0
;	smg.vel -= medvel 
	;smg.vel = smg.vel - median(smg.vel,/even)
	self->calc_rms, smg
;	smg.vel += medVel

	;observatory, 'keck', obs																		;get the observatory information
	;jd = ?
	;helio = heliocentric(smg.ra, smg.dec, jd=jd, longitude=obs.longitude, $					;compute the heliocentric correction
	;	latitude=obs.latitude, altitude=obs.altitude)
	;smg.vel -= helio
	;if tag_exist(smg, 'helio_corr') then stop 	;i'm adding this tag to deimos_spec__define.pro, when it gets propagated, this should trigger and I will have to deal with it	
	smg = struct_trimtags(smg, except=['X','Y','R','PA'])
	smg = struct_addtags(smg, arr_struct(projdist2(smg.ra, smg.dec, self.galra, self.galdec, self.parot)))
;  if self.name eq 'NGC3377' or self.name eq 'NGC4697' then smg.x = -smg.x
 
	if ((self.kinType eq 'H3') or (self.kinType eq 'H4')) then begin					;if h3 or h4 is selected
		dummy = execute('igood = where(abs(smg.'+self.kinType+') ne 0.3,ngood)')		;make sure that a limit is not hit
		if (ngood eq 0) then stop
		smg = smg[igood]
	endif
	
if vs_fitting eq 1 then begin	
	;----------------find v_sys using kinemetry
 inner = 1
 outer = 0
; vs_fitting = 1
 data = smg
 sauron = 0
 smeagol = 1
 reff = self.reff
 name = self.name
 TwoBin_Kinemetry, data, reff, vsys, diagnostics, name, inner, outer, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
    nopne=nopne, nogcs=nogcs, vs_fit, dvs_fit, sauron, smeagol, N0821, N1023, N1344, N2768, N3115, N3377, N4473, N4564, N4697, inplotdat, inmaxr, bin1, l, nbootstraps, twobin=twobin,$
    obser
    
    
 undefine,data
endif

 if vs_fitting eq 1 then offsetVel = vs_fit
 if vs_fitting eq 0 then offsetVel = medvel
;-------------------------------- 
offsetVel = medvel
smg.vel -= offsetVel
if vs_fitting eq 1 then print,strcompress('SMEAGOL fitted v_sys ='+string(vs_fit)+'+/-'+string(dvs_fit))
if vs_fitting eq 0 then print,strcompress('SMEAGOL photometric v_sys ='+string(medvel))

	
	
;	self->get_Flux, smg
	return, smg
end
;-----------------------------------------------------------------------------------------------

pro map_galaxy::point_symmetrize
	compile_opt idl2, hidden
	types = ['discrete','stellar']
	foreach type, types do begin
		dummy = execute('data = *(self.'+type+')')
		if (n_elements(data) gt 0) then begin
			d0 = data
			ind = indgen(n_elements(d0))
			d1 = d0[ind]
			d1.x = -d1.x
			d1.y = -d1.y
		;	d1.vel = -d1.vel
			if tag_exist(d1,'H3') then d1.h3 = -d1.h3
			if tag_exist(d1,'H4') then d1.h4 = -d1.h4	
			d = struct_append(d0, d1)
			dummy = execute('*(self.'+type+') = d')
		endif
	endforeach
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::get_dss, invert=invert, log=log
	compile_opt idl2, hidden
	
	s = *(self.plot)
	
	d = getdss(self.name)
	getrot, d.dsshdr, imrot, cdelt											;	extract the plate scale
	;adxy, d.dsshdr, d.ra, d.dec, xcen, ycen
	
	mmm, d.dssimage, mode
	img = d.dssimage
	;img = d.dssimage - mode
	dim = size(img,/dimensions)
	xcoord = rebin(indgen(dim[0]), dim)
	ycoord =  transpose(rebin(indgen(dim[1]), reverse(dim)))
	xyad, d.dsshdr, xcoord, ycoord, ra, dec
	ind = lindgen(n_elements(ra))
	obj = arr_struct({ra:ra[ind], dec:dec[ind], val:img[ind]})
	obj = struct_addtags(obj, arr_struct(projdist2(obj.ra, obj.dec, self.galra, self.galdec, self.parot)))
	if (n_elements(invert) gt 0) then obj.val = -obj.val
	if (n_elements(log) gt 0) then obj.val = alog10(obj.val) else obj.val -= 1.*mode
	
	if self.doreff then norm = self.reff else norm = 1d/self.arc2kpc

if 0 then begin
	mult = 1.2
	ind = where( (obj.x/norm gt mult*s.xrange[0]) and (obj.x/norm lt mult*s.xrange[1]) and (obj.y/norm gt mult*s.yrange[0]) and (obj.y/norm lt mult*s.yrange[1]), count)
	if (count eq 0) then stop
	obj = obj[ind]	
endif

	obj.x /= norm
	obj.y /= norm
	*(self.dss) = obj
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::change_units
	compile_opt idl2, hidden

	obj = *(self.obj)
	if (n_elements(obj) eq 0) then begin
		discrete = *(self.discrete)
		if (n_elements(discrete) gt 0) then obj = struct_append(obj, discrete)
		stellar = *(self.stellar)
		if (n_elements(stellar) gt 0) then obj = struct_append(obj, stellar)
	endif

	if self.doreff then begin
		obj.x /= self.reff
		obj.y /= self.reff
		self.maxDist /=self.reff
	endif else begin
		obj.x *= self.arc2kpc
		obj.y *= self.arc2kpc
		self.maxDist *= self.arc2kpc
	endelse
	if tag_exist(obj, 'VAL') then obj.val /= self.norm	
	
	*(self.obj) = obj

end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::map_obj, sneb=sneb, dneb=dneb
	compile_opt idl2, hidden

	if (n_elements(sneb) eq 0) then strneb = 0d else strneb = sneb
	if (n_elements(dneb) eq 0) then disneb = 10d else disneb = dneb	
	if (self.kintype ne 'VEL') then dummy = execute(self.kinType+' = 1')
	

	dss = *(self.dss)
	if (n_elements(dss) gt 0) then begin
		obj = dss																;take the point-symmetrizing out of map_obj!
		obj = struct_addtags(obj, replicate({type:'stellar'},n_elements(obj)))
	endif else begin
		stellar = *(self.stellar)
		if (n_elements(stellar) gt 0) then $
			obj = map_obj(stellar, strneb, veldisp=veldisp, h3=h3, h4=h4)		;take the point-symmetrizing out of map_obj!

		discrete = *(self.discrete)
		if (n_elements(discrete) gt 0) then $
			obj = struct_append(obj, map_obj(discrete, disneb))
	endelse

	obj.x += randomn(seed,n_elements(obj))*1d-3									;this ensures that triangulate will not fail
	obj.y += randomn(seed,n_elements(obj))*1d-3	

	*(self.obj) = obj
	
	stop

end

;-----------------------------------------------------------------------------------------------


;pro map_galaxy::get_variogram, d, range, ftype, axisratio, inc, varplot=varplot
pro map_galaxy::get_variogram, varplot=varplot

	compile_opt idl2, hidden

	;generate semi-variogram
	;svar_info = semi_variogram(d.x, d.y, d.z, inc=inc, range=range, ftype=ftype, axisratio=axisratio, plot=varplot)	

	svar_info = *(self.svar_info)
	d = svar_info.data
	svar_info = struct_trimtags(svar_info, except=['data'])

;if self.name eq 'NGC821' then svar_info.ftype = 'linear' & svar_info.range = 250.0  ;set sigma_e, sigma_kpc from sluggs survey paper


	
	;self->map_galaxy::get_variogram, d, range, ftype, axisratio, inc, varplot=varplot

;varplot = 1	
	
	if (self.varplot eq 1) then varplot = 1
	svar_info = struct_addtags(svar_info, $
		semi_variogram(d.x, d.y, d.z, svar_info=svar_info, plot=varplot));inc=inc, range=range, ftype=ftype, axisratio=axisratio, plot=varplot)	
	*(self.svar_info) = svar_info


	;generate error bars for the discrete data
	discrete = *(self.discrete)			
	if (n_elements(discrete) gt 0) then begin
		;discrete = struct_addtags(discrete, arr_struct({err:(discrete.errvel^2d), type:strarr(n_elements(discrete))+'discrete'}))	;add the velocity dispersion in quadrature, to act as effective error, turn up or down to see effect
	
		;discrete = struct_addtags(struct_trimtags(discrete, select=['X','Y','ERRVEL']), arr_struct({z:discrete.vel, err:discrete.errvel^2d, type:strarr(n_elements(discrete))+'discrete'}))
		dummy = execute('hold = arr_struct({z:discrete.'+self.kinType+', err:discrete.err'+self.kinType+'^2d, type:strarr(n_elements(discrete))+"discrete"})')
		discrete = struct_addtags(struct_trimtags(discrete, select=['X','Y','ERRVEL','ERRVELDISP','ERRH3','ERRH4','ERRRMS','FLUX','ERRFLUX']), hold)
		k = 0
		while 0 do begin
			;dummy = var_locmean(discrete, svar_info, mapsize=mapsize, nthreads=0, $
			;	maxdist=self.maxDist*50., maxPoints=self.maxPoints, inter_data=inter_data, /justdata)
			dummy = var_locmean(discrete, svar_info, mapsize=mapsize, nthreads=0, $
				maxdist=1d10, maxPoints=1d10, inter_data=inter_data, /justdata)
			if ((min(inter_data.varz) eq 0) or (min(finite(inter_data.varz)) eq 0)) then stop
			discrete.err = inter_data.varz
			;plot, discrete.x, sqrt(discrete.err), psym=1
			;pause
			++k
			if (k gt 10) then break
		endwhile	
	endif
	;generate error bars for the stellar data
	stellar = *(self.stellar)						
	if (n_elements(stellar) gt 0) then begin	
		;stellar = struct_addtags(struct_trimtags(stellar, select=['X','Y','ERRVEL']), arr_struct({z:stellar.vel, err:0^2d + stellar.errvel^2d, type:strarr(n_elements(stellar))+'stellar'}))
		dummy = execute('hold = arr_struct({z:stellar.'+self.kinType+', err:stellar.err'+self.kinType+'^2d, type:strarr(n_elements(stellar))+"stellar"})')
		stellar = struct_addtags(struct_trimtags(stellar, select=['X','Y','ERRVEL','ERRVELDISP','ERRH3','ERRH4','ERRRMS','FLUX','ERRFLUX']), hold)
		k = 0
		while 0 do begin
			;dummy = var_locmean(stellar, svar_info, mapsize=mapsize, nthreads=0, $
			;	maxdist=self.maxDist, maxPoints=self.maxPoints, inter_data=inter_data, /justdata)
			dummy = var_locmean(stellar, svar_info, mapsize=mapsize, nthreads=0, $
				maxdist=1d10, maxPoints=1d5, inter_data=inter_data, /justdata)
			if ((min(inter_data.varz) eq 0) or (min(finite(inter_data.varz)) eq 0)) then stop			
			stellar.err = inter_data.varz
			plot, stellar.x, sqrt(stellar.err), psym=1
			pause
			++k
			if (k gt 10) then break
		endwhile
	endif
	
	
	if (n_elements(discrete) gt 0) then obj = discrete
	if (n_elements(stellar) gt 0) then obj = struct_append(obj, stellar)
	nobj = n_elements(obj)
	obj = struct_addtags(obj, arr_struct({pa:pacalc(obj.x, obj.y)*!dtor}))

	if 0 then begin
		if self.doreff then begin
			obj.x /= self.reff
			obj.y /= self.reff
			self.maxDist /=self.reff
		endif else begin
			obj.x *= self.arc2kpc
			obj.y *= self.arc2kpc
			;self.maxDist *= self.arc2kpc
		endelse
	endif
	obj.z /= self.norm

	*(self.obj) = obj

end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::plotting_init, s, pos=pos, current=current
	compile_opt idl2, hidden

	self.pos = pos
	if (n_elements(current) gt 0) then *(self.current) = 1
	
	if (n_elements(s) eq 0) then stop
	
	obj = *(self.obj)

	vals = create_struct( $
			'xrange'		, [-1,1], $;max(abs(obj.x))*[-1,1], $
			'yrange'		, [-1,1], $	'max(abs(obj.y))*[-1,1], $
			'xsiz'			, 100		, $
			'rnorm'			, 1.		, $
			'xticknames'	, ''		, $
			'yticknames'	, ''		, $
			'xtitle'		, ''		, $
			'ytitle'		, ''		, $
			'title'			, ''		, $
			'xminor'		, 1			, $
			'yminor'		, 1			, $
			'min_points'	, 150		, $
			'method'		, 'kriging'	, $
			'smooth'		, 15		, $
			'yticklen'		, 0.015		, $
			'xticklen'		, 0.015		, $
			'fill'			, 255.		, $
			'tickfont_size'	, 15		  $
		)

	if (n_elements(reff) eq 0) then reff = 1.			;????
	
	if tag_exist(s, 'dimensions') then self.dimensions = s.dimensions
	valtags = tag_names(vals)
	for i=0,(n_tags(vals)-1) do if ~tag_exist(s, valtags[i]) then s = create_struct(valtags[i], vals.(i), s)
	
	if tag_exist(obj,'VAl') then begin
		if ~tag_exist(s,'range') then begin
			range = max(abs(obj.val))*[-1,1]
			if (min(obj.val) gt 0) then range[0] = 0
			s = struct_addtags(s, {range:range})
			undefine, range
		endif
	endif
	
  if self.name eq 'NGC821' then sige= 183d  ;set sigma_e, sigma_kpc from sluggs survey paper
  if self.name eq 'NGC1023' then sige = 180d
  if self.name eq 'NGC2768' then sige = 202d
  if self.name eq 'NGC3115' then sige = 200d
  if self.name eq 'NGC3377' then sige = 132d
  if self.name eq 'ngc4564' then sige = 151d
  if self.name eq 'NGC4697' then sige = 180d
  
  s.range = [-1,1]*sige

	ysiz = fix(difference(s.yrange)/float(difference(s.xrange))*s.xsiz)
	if tag_exist(s, 'ysiz') then s.ysiz = ysiz else s = struct_addtags(s, {ysiz:ysiz})
	undefine, ysiz
	
	s.min_points = min([n_elements(obj)-1,s.min_points])

	self.xmid = round((s.xsiz-1)/2.)												;x-midpoint of the map
	self.ymid = round((s.ysiz-1)/2.)												;y-midpoint of the map

	
	*(self.plot) = s

	
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::make_grid, new=new, z=z, varz=varz, minerr=minerr, noextrap=noextrap, holdim=holdim, fill=fill
  compile_opt idl2, hidden

  obj = *(self.obj)
  s = *(self.plot)
  s.fill = s.range[1]
  *(self.plot) = s
  gsm = *(self.gsm)
;  s.xsiz=500d
;  s.ysiz=285.7d
;  s.xsiz = 50d  
;  s.ysiz = 28.6d  ;maintain aspect ratio: ysiz/xsiz = 286/500 = 4/7 = yrange/xrange
  ;s.ysiz = fix(difference(s.yrange)/float(difference(s.xrange))*s.xsiz)
  s.ysiz = s.xsiz * (4/7d)
  s.smooth = 0
  
  if (n_elements(holdim) eq 0) then begin
    
    holdim = fltarr(s.xsiz, s.ysiz) + s.fill                ;create an array to hold the map, and make the background gray

    if (n_elements(fill) gt 0) then begin
      *(self.image) = holdim
      return
    endif
    
;   if ((n_elements(new) gt 0) and (n_elements(*(self.dss)) eq 0)) then begin
    if ((n_elements(new) gt 0) ) then begin

      ;min_nrange = self.min_nrange
      
      svar_info = *(self.svar_info)

      s.xrange = [-7,7];
;      s.xrange = [-5,5]
      s.yrange = [-4,4];-4,4
;s.xrange=[-6,6]
;s.yrange=[-4.8,4.8]
      
      xd = findgen(s.xsiz)/double(s.xsiz - 1d)*difference(s.xrange) + s.xrange[0]
      yd = findgen(s.ysiz)/double(s.ysiz - 1d)*difference(s.yrange) + s.yrange[0]
      xd = rebin(xd, s.xsiz, s.ysiz)                        ;convert to a 2d image
      yd = transpose(rebin(yd, s.ysiz, s.xsiz))
      
;      xd = xd*2.
;      yd = yd*2.


      if self.doreff then begin
        xd *= self.reff
        yd *= self.reff
      endif else begin
        xd /= self.arc2kpc
        yd /= self.arc2kpc
      endelse
      indx = lindgen(n_elements(xd))
      ;loc = [[[xd]],[[yd]]]
      loc = replicate({x:0d, y:0d}, n_elements(xd))
      loc.x = xd[indx]
      loc.y = yd[indx]
  
        
      obj = *(self.obj)

      if 0*(self.kintype eq 'FLUX') then begin
        if self.doreff then begin
          obj.x *= self.reff
          obj.y *= self.reff
        endif else begin
          obj.x /= self.arc2kpc
          obj.x /= self.arc2kpc
          stop
        endelse
        ibad = where(match eq -1,nbad)
        if (nbad gt 0) then stop
        ind = where((match ne -1),ngood)
        obj = obj[match]
      
        grid0 = fltarr(s.xsiz,s.ysiz)                         ;create an array of the output size
        for i=0l,(s.ysiz-1) do grid0[0,i] = double(obj[i*s.xsiz:(i+1l)*s.xsiz-1l].val)    ;fill it with the interpolated values
      endif else begin
        grid0 = var_locmean_4(obj, svar_info, mapsize=[s.xsiz,s.ysiz],xrange=2.*minmax(obj.x),yrange=2.*minmax(obj.y), nthreads=0, varim=varim, $
          maxdist=1d10, maxPoints=1d4, nrange=nrange, loc=loc, minerr=minerr, noextrap=noextrap, gsm=gsm)
      endelse
      
      
      
  
      grid = grid0
      
      if (n_elements(varim) gt 0) then varz = varim
      ;bPix = where( (varz le -1) or ~finite(varz) or ~finite(grid) or ((nrange[*,*,0] lt min_nrange[0]) and (nrange[*,*,1] lt min_nrange[1])), nbad)
      bPix = where( (grid eq 0.0), nbad)
      if (n_elements(varz) gt 0) then if (nbad gt 0) then varz[bPix] = !values.f_nan
      
      s.smooth = s.smooth < (min(size(grid,/dim)) - 2)

      dgrid = smooth(grid,s.smooth,/nan);,/edge)
;     if (nbad gt 0) then dgrid[bPix] = !values.f_nan

      holdim = dgrid
      
      z = dgrid
      
      
  
;     display, z, /aspect, min=-75, max=75, top=254, /silent
;     pause
      
      if 0 then begin
        window, 2, xpos=2900
        cleanplot, /silent
        erase
        device, decomposed=0
        ;sauron_colormap, rgb_table=rgb_table
        loadct, 33, /silent
        tvlct, r, g, b, /get
        rgb_table = [[r],[g],[b]] 
        ;loadct, 0, /silent
        
        display, dgrid, /aspect, min=s.range[0], max=s.range[1], top=254
        
        window, 0, xpos=1950
        cleanplot, /silent
        erase
        display, sqrt(varz), /aspect, min=0, max=160, top=254
      
        ;stop
      endif   
    
      ;return
  
    endif else begin
      ;s.xrange = [-1,1]*100.

      grid = map_grid(obj, holdim, s.xsiz, s.xrange, dim=dim, smooth=s.smooth, fill=s.fill, out=out, method=s.method, min_points=s.min_points)
    
      imloc = [n_elements(holdim[*,0]) - n_elements(grid[*,0]),n_elements(holdim[0,*]) - n_elements(grid[0,*])]/2 
      holdim[imloc[0],imloc[1]] = grid*1.
    endelse
  endif

 

  *(self.image) = holdim
    
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::plot, p2=p2, nox=nox, noy=noy, s=s, fill=fill, galname
  compile_opt idl2, hidden

  holdim = *(self.image)
  s = *(self.plot)
  if (n_elements(fill) gt 0) then holdim = holdim*0. + s.fill
  pos = self.pos
  current = *(self.current)
  print,strcompress('range ='+string(s.range))
  s.xrange=[-7,7]
;  s.xrange=[-5,5]
  s.yrange=[-4,4]
;s.xrange=[-6,6]
;s.yrange=[-4.8,4.8]

  ihigh = where(holdim gt (s.range[1] - (s.range[1] - s.range[0])/255.),nhigh)

  ibad = where(~finite(holdim),nbad)
  if (nbad gt 0) then holdim[ibad] = s.range[1]
  
  if (nHigh gt 0) then holdim[ihigh] = s.range[1] - (s.range[1] - s.range[0])/255.

; p2 = image(holdim, min_value=s.range[0], max_value=s.range[1], yticks=1, xticks=1, xminor=1, xtickname=[' ',' '], yminor=1, $
;   ytickname=[' ',' '], xsubticklen=0, xticklen=s.xticklen, rgb_table=s.rgb_table, axis_style=0, position=pos, current=current, $
;   dimension=self.dimensions, tickfont_size=s.tickfont_size, yticklen=s.yticklen, title=s.title)
;value of self.dimensions?
 p2 = image(holdim, min_value=s.range[0], max_value=s.range[1], xminor=1, xtickname=[' ',' '], yminor=1, $
   ytickname=[' ',' '], xsubticklen=0, xticklen=s.xticklen, rgb_table=s.rgb_table, axis_style=0, position=pos, current=current, $
   dimension=self.dimensions,yticklen=s.yticklen, title=s.title)

;s.xsiz = 500;357.;
;s.ysiz = 285.7
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
  epsilon1=sqrt(1-0.677^2);eccentricity=sqrt(1-q)
  epsilon4=sqrt(1-0.733^2)
;  major1=(self.reff)/sqrt(0.677)
;  major4=(4*self.reff)/sqrt(0.733)
if galname eq 'ngc0821' then begin
  major1=43.292609069
  minor1=29.4389741669
  major4=166.7923855084
  minor4=122.2588185776
endif
if galname eq 'ngc1023' then begin
  major1=55.7540329942
  minor1=22.8591535276
  major4=238
  minor4=85.68
endif
if galname eq 'ngc1344' then begin
  major1=92d
  minor1=57d
  major4=184d
  minor4=184d
endif
if galname eq 'ngc2768' then begin
  major1=55.0862959365
  minor1=23.1362442933
  major4=206.1140461007
  minor4=98.9347421283
endif
if galname eq 'ngc3115' then begin
  major1=66.2932357292
  minor1=19.2250383615
  major4=238
  minor4=85.68
endif
if galname eq 'ngc3377' then begin
  major1=48.5815465652
  minor1=26.2340351452
  major4=174.4579665408
  minor4=116.8868375823
endif
if galname eq 'ngc4564' then begin
  major1=52.0738019647
  minor1=24.4746869234
  major4=174.4579665408
  minor4=116.8868375823
endif
if galname eq 'ngc4697' then begin
  major1=44.625
  minor1=28.56
  major4=170.678645413
  minor4=119.4750517891
endif

;  p=ellipse(xmid,ymid,minor=minor1,major=major1,theta=0,fill_background=0,target=p2,/data)  ;inner isophote
;  if galname eq 'ngc0821' then p=ellipse(xmid,ymid,minor=minor4,major=major4,theta=-1.5d,fill_background=0,target=p2,/data) $ ;outer isophote
;  else p=ellipse(xmid,ymid,minor=minor4,major=major4,fill_background=0,target=p2,/data)  ;outer isophote
;  p2.Refresh, /DISABLE




 ticknames = s.xticknames_latex
;ticknames = s.xticknames
 pr = strpos(ticknames[0],'.')
 if (strpos(ticknames[0],'.') eq -1) then pr = 0 else pr = strlen(ticknames[0]) - strsplit(ticknames[0],'.') - 1
; s.xsiz = 500.;357.1;
; s.ysiz = 285.71
ticknames = strtrim(['-6','-4','-2','0','2','4','6'],1)
;ticknames = strtrim(['-4','-2','0','2','4'],1)
tickvalues = indgen(n_elements(ticknames))
for i=0,n_elements(tickvalues)-1 do begin
tickvalues[i] = (s.xsiz/2)*(1+ticknames[i]/s.xrange[1])
endfor

;  s.xticknames_latex = roundx(float(ticknames) > s.xrange[0] < s.xrange[1],pr)
;  tickvalues = s.xticknames_latex*float(s.xsiz)/float(difference(s.xrange)) + self.xmid

;  s.xticknames = roundx(float(ticknames) > s.xrange[0] < s.xrange[1],pr)
;  tickvalues = s.xticknames*float(s.xsiz)/float(difference(s.xrange)) + self.xmid
;  tickvalues = ticknames*float(s.xsiz)/float(difference(s.xrange)) + self.xmid

;  tickvalues = [0,83.333328,       166.66666 ,      250.00000,       333.33334,       416.66667, 500]
;  tickvalues = [83.333328,       166.66666 ,      250.00000,       333.33334,       416.66667]

;tickvalues = [1.16,4.14,7.12,10.10,13.08,16.0600, 19.0400,22.0200,25, 27.9800, 30.9600, 33.9400,36.9200,39.9000,42.8800,45.8600,48.8400]
;tickvalues = [3.57,10.71,17.86,25,32.19,39.29,46.43]*10.
;fr = (2./5.)*178.5
;tickvalues = [178.5-2*fr,178.5-fr,178.5,178.5+fr,178.5+2.*fr]

;v = indgen(7)
;for i=0,6 do begin
;    v[i] = (0.143)*250 + i*(0.286)*250
;endfor

;   tickvalues = 30.*(0.5*s.xsiz/(minmax(d.x)-self.galra)*0.5)*ticknames
;s.xsiz = 500.;357.1;
;s.ysiz = 285.71
topaxis = axis('X', location=[0,s.ysiz], showtext=0, tickdir=1, target=p2, textpos=0, ticklen=s.xticklen, subticklen=0.5, $
         tickfont_size=s.tickfont_size, tickvalues=tickvalues, tickname=strarr(n_elements(tickvalues))+' ', minor=s.xminor)
;topaxis = axis('X', location=[0,0], tickdir=0, target=p2, textpos=0, ticklen=s.xticklen, subticklen=0.5, $
;          tickfont_size=s.tickfont_size, tickvalues=tickvalues, title = '$\itX/$$R_e$', tickname=ticknames, minor=s.xminor)
;
  if (nox) then begin
    topaxis = axis('X', location=[0,0], tickdir=0, showtext=0, target=p2, $
          ticklen=s.xticklen, subticklen=0.5, tickfont_size=s.tickfont_size, tickvalues=tickvalues, $
          minor=s.xminor, title='') 
  endif else begin   
    topaxis = axis('X', location=[0,0], tickdir=0, target=p2, textpos=0, ticklen=s.xticklen, $
          subticklen=0.5, tickfont_size=s.tickfont_size, tickvalues=tickvalues, tickfont_style=1, tickname=ticknames, minor=s.xminor, title='$\itX/R_{\rm e}$')
  endelse
  
;topaxis = axis('X', location=[0,0], tickdir=0, target=p2, textpos=0, ticklen=s.xticklen, $
;                  subticklen=0.5, tickfont_size=s.tickfont_size, tickvalues=tickvalues, tickfont_style=1, tickname=ticknames, minor=s.xminor, title='$\itX/$$R_e$') 
             



 pr = strpos(ticknames[0],'.')
 if (strpos(ticknames[0],'.') eq -1) then pr = 0 else pr = strlen(ticknames[0]) - strsplit(ticknames[0],'.') - 1
; s.yticknames_latex = roundx(float(ticknames) > s.yrange[0] < s.yrange[1],pr)
; tickvalues = s.yticknames_latex*float(s.ysiz)/float(difference(s.yrange)) + self.ymid
;  s.yticknames = roundx(float(ticknames) > s.yrange[0] < s.yrange[1],pr)
;  tickvalues = s.yticknames*float(s.ysiz)/float(difference(s.yrange)) + self.ymid


;   ticknames = ['-4','-2','0','2','4']
;   tickvalues = [0.5*s.ysiz+((ticknames[0]/s.yrange[1])*(s.ysiz/2.)),0.5*s.ysiz+((ticknames[1]/s.yrange[1])*(s.ysiz/2.)),$
;    0.5*s.ysiz+((ticknames[2]/s.yrange[1])*(s.ysiz/2.)),0.5*s.ysiz+((ticknames[3]/s.yrange[1])*(s.ysiz/2.)), $
;    0.5*s.ysiz+((ticknames[4]/s.yrange[1])*(s.ysiz/2.))]
;   tickvalues = [0+7.13/2.,7.13/2.,14.25,21.38/2.,28.5-7.13/2.]*10.
;  tickvalues = [3, 6, 9, 10., 12.96, 15.92, 18.88];*10
;  tickvalues = [0.25*(s.ysiz/2.), (s.ysiz/2.), (s.ysiz/2.)+0.75*(s.ysiz/2.)]

ticknames = strtrim(['-4','-2','0','2','4'],1)
;ticknames = strtrim(['-3','0','3'],1)
tickvalues = indgen(n_elements(ticknames))
for i=0,n_elements(tickvalues)-1 do begin
tickvalues[i] = (s.ysiz/2)*(1+ticknames[i]/s.yrange[1])
endfor

;  yaxis = axis('Y', location=[s.xsiz,0], tickdir=1, position=pos, target=p2, ticklen=s.yticklen, subticklen=0.5, tickfont_size=s.tickfont_size, $
  yaxis = axis('Y', location=[s.xsiz,0], showtext=0, tickdir=1, target=p2, ticklen=s.yticklen, subticklen=0.5, tickfont_size=s.tickfont_size, $
         tickvalues=tickvalues, tickname=strarr(n_elements(tickvalues))+' ', major=1, minor=s.yminor)
;  if (n_elements(noy) eq 0) then yaxis = axis('Y', location=[0,0], tickdir=0, position=pos, target=p2, ticklen=s.yticklen, subticklen=0.5, $
  if (noy) then begin
    yaxis = axis('Y', location=[0,0], tickdir=0, target=p2, ticklen=s.yticklen, subticklen=0.5, $
         tickvalues=tickvalues, showtext=0, tickfont_size=s.tickfont_size, major=1, title=s.ytitle, minor=s.yminor) 
  endif else begin
    yaxis = axis('Y', location=[0,0], tickdir=0, target=p2, textpos=0, ticklen=s.xticklen, subticklen=0.5, $
         tickfont_size=s.tickfont_size, tickvalues=tickvalues, tickfont_style=1, title = '$\itY/R_{\rm e}$', tickname=strtrim(ticknames,1), minor=s.xminor)
  endelse
  

end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::outline_data, p2=p2, bb=bb
	compile_opt idl2, hidden

	obj = *(self.obj)
	s = *(self.plot)
	pos = self.pos




	triangulate, obj.x, obj.y, tr, b

	if (n_elements(bb) gt 0) then begin
		obj = bb
		b = indgen(n_elements(obj))
	endif

	xbnd = obj[b].x*s.xsiz/float(difference(minmax(s.xrange))) + self.xmid
	;xbnd = [xbnd, xbnd[0]]
	ybnd = obj[b].y*s.xsiz/float(difference(minmax(s.xrange))) + self.ymid	
	;ybnd = [ybnd, ybnd[0]]
	
	;xbnd = 0 > xbnd < s.xsiz
	;ybnd = 0 > ybnd < s.ysiz
		
	;cycle through the points, if a point is outside the boundary, then spawn an extra point in front and behind on the boundary
	;then remove the offending point
	;there are probably unaccounted for situations where this will fail
	;just need to add in extra if statements to take care of those as they arise
	
	d = arr_struct({x:xbnd, y:ybnd})
	d = struct_addtags(d, arr_struct({pa:pacalc(d.x - self.xmid, d.y - self.ymid, 0., /radians)}))
	d = d[reverse(sort(d.pa))]
	d = struct_trimtags(d, except='PA')
	d = [d,d[0]]
	d0 = d
	
	loadct, 0
	flag = 1
	while flag do begin
		for i=0,(n_elements(d) - 1) do begin
			if ((d[i].x lt 0) or (d[i].y lt 0) or (d[i].x gt s.xsiz) or (d[i].y gt s.ysiz)) then begin
				flag = 1
	
				if (d[i].y lt 0) then begin														;if the point is below the bottom line
					if (i eq (n_elements(d)-1)) then ind = 0 else ind = i + 1
					loc = intlin([0,s.xsiz,d[ind].x,d[i].x], [0,0,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = d[i] else newd = d[0:i]							;isolate the points before and including the out of bounds point
						newd = struct_append(newd, {x:loc[0], y:loc[1]})						;add on the new point on the boundary
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif
					if (i eq 0) then ind = -1 else ind = i - 1
					loc = intlin([0,s.xsiz,d[ind].x,d[i].x], [0,0,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = {x:loc[0], y:loc[1]} else newd = struct_append(d[0:i-1], {x:loc[0], y:loc[1]})		;isolate the points before the out of bounds point
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif else d = d[where(indgen(n_Elements(d)) ne i)]
					break
				endif
				
				if (d[i].y gt s.ysiz) then begin														;if the point is above the top line
					if (i eq (n_elements(d)-1)) then ind = 0 else ind = i + 1
					loc = intlin([0,s.xsiz,d[ind].x,d[i].x], [s.ysiz,s.ysiz,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = d[i] else newd = d[0:i]							;isolate the points before and including the out of bounds point
						newd = struct_append(newd, {x:loc[0], y:loc[1]})						;add on the new point on the boundary
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif
					if (i eq 0) then ind = -1 else ind = i - 1
					loc = intlin([0,s.xsiz,d[ind].x,d[i].x], [s.ysiz,s.ysiz,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = {x:loc[0], y:loc[1]} else newd = struct_append(d[0:i-1], {x:loc[0], y:loc[1]})		;isolate the points before the out of bounds point
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif else d = d[where(indgen(n_Elements(d)) ne i)]
					break
				endif			
	
				if (d[i].x lt 0) then begin														;if the point is left of the left line
					if (i eq (n_elements(d)-1)) then ind = 0 else ind = i + 1
					loc = intlin([0,0,d[ind].x,d[i].x], [0,s.ysiz,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = d[i] else newd = d[0:i]							;isolate the points before and including the out of bounds point
						newd = struct_append(newd, {x:loc[0], y:loc[1]})						;add on the new point on the boundary
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif
					if (i eq 0) then ind = -1 else ind = i - 1
					loc = intlin([0,0,d[ind].x,d[i].x], [0,s.ysiz,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = {x:loc[0], y:loc[1]} else newd = struct_append(d[0:i-1], {x:loc[0], y:loc[1]})		;isolate the points before the out of bounds point
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif else d = d[where(indgen(n_Elements(d)) ne i)]
					break
				endif	
	
				if (d[i].x gt s.xsiz) then begin														;if the point is left of the left line
					if (i eq (n_elements(d)-1)) then ind = 0 else ind = i + 1
					loc = intlin([s.xsiz,s.xsiz,d[ind].x,d[i].x], [0,s.ysiz,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = d[i] else newd = d[0:i]							;isolate the points before and including the out of bounds point
						newd = struct_append(newd, {x:loc[0], y:loc[1]})						;add on the new point on the boundary
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif
					if (i eq 0) then ind = -1 else ind = i - 1
					loc = intlin([s.xsiz,s.xsiz,d[ind].x,d[i].x], [0,s.ysiz,d[ind].y,d[i].y])
					if (((loc[0] ge d[ind].x) and (loc[0] le d[i].x)) or  ((loc[0] le d[ind].x) and (loc[0] ge d[i].x))) then begin
						if (i eq 0) then newd = {x:loc[0], y:loc[1]} else newd = struct_append(d[0:i-1], {x:loc[0], y:loc[1]})		;isolate the points before the out of bounds point
						if (i lt (n_Elements(d) - 1)) then newd = struct_append(newd, d[i+1:*]) ;if this is not the last point, then add on the remaining points
						d = newd
					endif else d = d[where(indgen(n_Elements(d)) ne i)]
					break
				endif	
			endif else flag = 0
		endfor
	endwhile
	
	p2 = plot(d.x, d.y, position=pos, /overplot, thick=s.c_thick)
	
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::plot_data_loc, p2=p2
	compile_opt idl2, hidden
	
	obj = *(self.obj)
	s = *(self.plot)
	pos = self.pos

	xobj = obj.x*s.xsiz/float(difference(minmax(s.xrange))) + self.xmid
	yobj = obj.y*s.xsiz/float(difference(minmax(s.xrange))) + self.ymid
		
	inside = where((xobj gt 0) and (yobj gt 0) and (xobj lt s.xsiz) and (yobj lt s.ysiz),count)
	if (count gt 0) then begin
		xobj = xobj[inside]
		yobj = yobj[inside]		
		p2 = plot(xobj, yobj, position=pos, /overplot, symbol='.', linestyle='')
	endif
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::sharpcorner, p2=p2
	compile_opt idl2, hidden
	sharpcorner, p2, thick=1, position=self.pos, dimension=self.dimensions
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::ellipsync, kpa=kpa, kq=kq, nodiscrete=nodiscrete, nostellar=nostellar, thisdevice=thisdevice
	compile_opt idl2, hidden

  if n_elements(nodiscrete) eq 1 then data = *(self.stellar)
  if n_elements(nostellar) eq 1 then data = *(self.discrete)

	;data.vel += 700.

	kpa = self.parot*!dtor
;	if (kpa gt !pi/2d) then kpa -= !pi
	vsys = 0.
;	kaxisratio = 0.55d
	kaxisratio = 1d

	if (n_elements(vel) eq 0) then start = randomu(seed)*300d else start = vel
	if (n_elements(kpa) eq 0) then start = [start, randomu(seed)*!pi - !pi/2d] else start = [start, kpa]
	if (n_elements(kaxisratio) eq 0) then start = [start, randomu(seed)] else start = [start, kaxisratio]
	if (n_elements(vsys) eq 0) then start = [start, 0d] else start = [start, vsys]
	
	start[1] = kpa
	
	undefine, kpa
;	undefine, kaxisratio ; Turning this off might eliminate anomalous/unstable solutions!
	
	pi = [	{fixed:n_elements(vel)			, limited:[1,1], limits:[-1000d,10000d]},	$	;rotational velocity
			{fixed:n_elements(kpa)			, limited:[1,1], limits:[-!pi/2d,!pi/2d]},	$	;kinematic position angle
			{fixed:n_elements(kaxisratio)	, limited:[1,1], limits:[0.1d,1.00001d]},		$	;kinematic axisratio
			{fixed:n_elements(vsys)			, limited:[1,1], limits:[-1d,10000d]}		]	;systemic velocity
		
	minNumSync = 20
	synciters = 10
	loc = where(tag_names(data) eq 'VEL')
	d = kinemetry_ellipsync(data, loc, synciters, start, pi, minNumSync=minNumSync, smkq=smkq)
	
	kpa = median(d.kpa,/even)
	kq = median(d.kq,/even)

;	undefine,pi
	
end






;-----------------------------------------------------------------------------------------------

pro map_galaxy::kinemetry, pp1=pp1, thisdevice=thisdevice, nosync=nosync, nodiscrete=nodiscrete, nostellar=nostellar, $
  dstellar=dstellar, ddiscrete=ddiscrete, fullrange=fullrange, innervel=innervel, innerpa=innerpa
	compile_opt idl2, hidden

	if n_elements(nodiscrete) eq 1 then begin
	 data = *(self.stellar)
	 nsamples = 1
  endif
	if n_elements(nostellar) eq 1 then begin 
	 data = *(self.discrete)
   nsamples=1
  endif
;	data = data[where(finite(data.pa))]
	vsys = 0d

  for k=0,nsamples-1 do begin ;start boostrap loop
  print,k
  if k ne 0 then begin
    bootind = boot_indices(n_elements(data))
    data = data[bootind]
  endif

	if (n_elements(nosync) gt 0) then begin
		kpa = self.parot*!dtor  
		if (kpa gt !pi/2d) then kpa -= !pi
		kaxisratio = self.q ;regular
	endif else begin
		self->ellipsync, kpa=kpa, kq=kaxisratio, nodiscrete=nodiscrete, nostellar=nostellar, thisdevice=thisdevice
		kaxisratio = kaxisratio < 0.95
	endelse

if self.name eq 'NGC821' then sige= 183d  ;set sigma_e, sigma_kpc from sluggs survey paper
if self.name eq 'NGC1023' then sige = 183d
if self.name eq 'NGC2768' then sige = 202d
if self.name eq 'NGC3115' then sige = 248d
if self.name eq 'NGC3377' then sige = 132d
if self.name eq 'ngc4564' then sige = 153d
if self.name eq 'NGC4697' then sige = 180d

if n_elements(nodiscrete) eq 1 then begin

  kaxisratio = 0.64d
  kpa = self.parot*!dtor
  vsys = 0 

	if (n_elements(vel) eq 0) then start = randomu(seed)*300d else start = vel
	if (n_elements(kpa) eq 0) then start = [start, randomu(seed)*!pi - !pi/2d] else start = [start, kpa]
	if (n_elements(kaxisratio) eq 0) then start = [start, 1d] else start = [start, kaxisratio]
	if (n_elements(vsys) eq 0) then start = [start, 0d] else start = [start, vsys]
	
	start[1] = kpa
	
	undefine,vel
	undefine,kPA
;	undefine,vsys
; undefine,kaxisratio
  
  
	pi = [	{fixed:n_elements(vel)			, limited:[1,1], limits:[-1000d,10000d]},	$	;rotational velocity
;			{fixed:n_elements(kpa)			, limited:[1,1], limits:[!pi/2d,!pi/2d]},	$	;kinematic position angle
      {fixed:n_elements(kpa)      , limited:[1,1], limits:[-4*!pi/2d,4*!pi/2d]},  $ ;kinematic position angle
			{fixed:n_elements(kaxisratio)	, limited:[1,1], limits:[0.1d,1.00001d]},		$	;kinematic axisratio
			{fixed:n_elements(vsys)			, limited:[1,1], limits:[-1d,10000d]}		]	;systemic velocity
			
			
endif	
;;-------------------------------------------------------------------------------------------------------INITIAL GUESSES w/ VELOCITY DISPERSION	
	if n_elements(nostellar) eq 1 then begin
  
  vsys=0d
  kaxisratio = 1.0;0.733 ;1d;
  veldisp = 183d  ;initial guess
  kpa = self.parot*!dtor
  vel = 0d

	if (n_elements(vel) eq 0) then start = randomu(seed)*300d else start = vel
  if (n_elements(kpa) eq 0) then start = [start, randomu(seed)*!pi - !pi/2d] else start = [start, kpa]
  if (n_elements(kaxisratio) eq 0) then start = [start, 1d] else start = [start, kaxisratio]
  if (n_elements(vsys) eq 0) then start = [start, 0d] else start = [start, vsys]
  if (n_elements(veldisp) eq 0) then start = [start, 150d] else start = [start, veldisp]
  if (n_elements(veldisp_sys) eq 0) then start = [start, 0d] else start = [start, veldisp_sys]
  start[1] = kpa

  undefine,vel
  undefine,kpa
  undefine,vsys
  undefine,veldisp
;  undefine,kaxisratio

  pi = [  {fixed:n_elements(vel)      , limited:[1,1], limits:[-1000d,10000d]}, $ ;rotational velocity
;      {fixed:n_elements(kpa)      , limited:[1,1], limits:[-!pi/2d,!pi/2d]},  $ ;kinematic position angle
       {fixed:n_elements(kpa)      , limited:[1,1], limits:[-4*!pi/2d,4*!pi/2d]},  $ ;kinematic position angle
      {fixed:n_elements(kaxisratio) , limited:[1,1], limits:[0.1d,1.00001d]},   $ ;kinematic axisratio
      {fixed:n_elements(vsys)     , limited:[1,1], limits:[-1d,10000d]},  $    ;systemic velocity
      {fixed:n_elements(veldisp)  , limited:[1,1], limits:[-300d,300d]},  $   ;velocity dispersion
      {fixed:n_elements(veldisp_sys)  , limited:[1,1], limits:[-1d,10000d]}   ]  ;systemic vdisp
  
  endif
;;-------------------------------------------------------------------------------------------------------  
  
  
;	startveldisp = [-50,start[1:2],100,0]
;	piveldisp = [{fixed:0	, limited:[1,1], limits:[-1000d,1000d]},	$	;rotational velocity
;				{fixed:1	, limited:[1,1], limits:[-!pi/2d,!pi/2d]},	$	;kinematic position angle;1
;				{fixed:1	, limited:[1,1], limits:[0.1d,1.00001d]},		$	;kinematic axisratio
;				{fixed:0	, limited:[1,1], limits:[10d,1000d]},		$	;systemic velocity
;				{fixed:1	, limited:[1,1], limits:[-1d,1d]}			]	;signifies evenness
;		
;	starth3 = [-0.1,start[1:2],0.0]
;	pih3 = [	{fixed:0	, limited:[1,1], limits:[-1d,1d]*0.3},	$	;rotational velocity
;				{fixed:1	, limited:[1,1], limits:[-!pi/2d,!pi/2d]},	$	;kinematic position angle
;				{fixed:1	, limited:[1,1], limits:[0.1d,1.00001d]},		$	;kinematic axisratio
;				{fixed:1	, limited:[1,1], limits:[-1d,1d]}		]	;systemic velocity
;
;	
;	starth4 = [0.1,start[1:2],0.0,0]
;	pih4 = [	{fixed:0	, limited:[1,1], limits:[-1d,1d]*0.3},	$	;rotational velocity
;				{fixed:1	, limited:[1,1], limits:[-!pi/2d,!pi/2d]},	$	;kinematic position angle
;				{fixed:1	, limited:[1,1], limits:[0.1d,1.00001d]},		$	;kinematic axisratio
;				{fixed:0	, limited:[1,1], limits:[-1d,1d]},			$	;systemic velocity
;				{fixed:1	, limited:[1,1], limits:[-1d,1d]}			]	;signifies evenness
;	




	iters = 10
	icIter = 10
	

	mininc = 5
	inc = 33 < n_elements(data)  ;25
	ind = lindgen(n_elements(data))							;indices of the data


  if n_elements(nodiscrete) eq 1 then begin
    kpa=self.parot*!dtor
    kaxisratio=self.q
	endif
	if n_elements(nostellar) eq 1 then begin
	  kpa=self.parot*!dtor
	  kaxisratio=1d;0.733
	endif
	dummy = ellippar(data.x, data.y, kaxisratio, kpa - self.parot*!dtor, struc=data, /sort)	;install new a and b values for each object
  undefine,kpa
  undefine,kaxisratio


	undefine, bin
	
	
;------------------------------------------------------------------------set up inital inner bin
	j = 0
	loc = where(tag_names(data) eq 'VEL')
	for j=0,(inc - mininc - 1) do begin
		use = ind[0:mininc+j-1]
;		bin = struct_append(bin, kinemetry_run(data, loc, use, err_vel, pi, start))
    if n_elements(nodiscrete) eq 1 then begin
		  hold = kinemetry_monte_skim(data[use], iters, pi, start, piH3=piH3, startH3=startH3, piH4=piH4, startH4=startH4, $
						piveldisp=piveldisp, startveldisp=startveldisp, icIter=icIter, novel=novel)
;		endif
;		if n_elements(nostellar) eq 1 then begin
;		  hold = kinemetry_max(data[use].pa, data[use].vel, data[use].errvel, pi, start, lam=lam, verbose=verbose)
;;		  hold = kinemetry_monte_discrete(data[use], iters, pi, start, piH3=piH3, startH3=startH3, piH4=piH4, startH4=startH4, $
;;          piveldisp=0, startveldisp=startveldisp, icIter=icIter, novel=novel)
;stop
;    endif
		hold = struct_addtags(hold, {er:mean(data[use].er)})
		bin = struct_append(bin, hold)
		endif
;		print, j
	endfor	
	
	
	
	
;------------------------------------------------------------------------Rolling bins algorithm

	j = 0
	while 1 do begin		
	     use = ind[j:min([j+inc,n_elements(data)])-1]						;isolate the data for this elliptical bin
		if (n_elements(use) lt mininc) then break
		if n_elements(nodiscrete) eq 1 then begin
;    if max(data[use].er/self.reff) GE 1.4 then begin
;      kaxisratio = 0.677
;      vsys = 0.0
;        pi = [  {fixed:0      , limited:[1,1], limits:[-1000d,10000d]}, $ ;rotational velocity
;      {fixed:0      , limited:[1,1], limits:[-!pi/2d,!pi/2d]},  $ ;kinematic position angle
;      {fixed:1 , limited:[1,1], limits:[0.1d,1.00001d]},   $ ;kinematic axisratio
;      {fixed:1     , limited:[1,1], limits:[-1d,10000d]}   ] ;systemic velocity
;    endif
		hold = kinemetry_monte_skim(data[use], iters, pi, start, piH3=piH3, startH3=startH3, piH4=piH4, startH4=startH4, $
						piveldisp=piveldisp, startveldisp=startveldisp, icIter=icIter, novel=novel)
		endif
    if n_elements(nostellar) eq 1 then begin
      hold = kinemetry_max(data[use].pa, data[use].vel, data[use].errvel, pi, start, lam=lam, verbose=verbose)
    endif
		hold = struct_addtags(hold, {er:mean(data[use].er)})
;		print,strcompress('Finished bin with ER ='+string(hold.er/self.reff))
		bin = struct_append(bin, hold)
		++j
	endwhile
	
reff=self.reff
dd =replicate({er:0d, v:0d, err:0d, disp:0d, PA:0d, tag:''},n_elements(bin.er))
dd.tag = k
dd.er = bin.er/reff
if n_elements(nodiscrete) eq 1 then dd.PA = bin.VROT_PA/!dtor mod 180d
if n_elements(nostellar) eq 1 then dd.PA = bin.PA/!dtor mod 180d

if n_elements(nodiscrete) eq 1 then begin
  dd.v = (bin.vrot+bin.vrot_sys)/sige
;  dd.disp = bin.veldisp/sige
;  p=plot(bin.er/self.reff, $
;    (bin.vrot+bin.vrot_sys)/sige, linestyle=6,symbol='plus', xrange=[0,8],yrange=[-0.5,1.5])
endif
if n_elements(nostellar) eq 1 then begin
  dd.v = (bin.quan)/sige
  dd.disp = bin.disp/sige
;  p=plot(bin.er/self.reff, (bin.quan)/sige, $
;    linestyle=6,symbol='plus', color='red', xrange=[0,8],yrange=[-0.5,1.5], /overplot)
endif
dd.err = bin.err_vrot[0]


if k eq 0 then begin
  d1 = dd
endif else begin
  d1 = struct_append(d1,dd)
  endelse

;	cleanplot, /silent
;	!p.multi=[0,1,2]
;	plot, bin.er, bin.vrot, title=self.name
;	if n_elements(nodiscrete) eq 1 then plot, bin.er, bin.veldisp_sys
;  if n_elements(nodiscrete) eq 1 then p=plot(bin.er/self.reff, (bin.vrot+bin.vrot_sys)/193d, linestyle=6, symbol='plus', xrange=[0,8],yrange=[-0.5,2],title='N821 Stellar',xtitle='$R/R_e$',ytitle='$\itV/\sigma_e$')
;  if n_elements(nostellar) eq 1 then p=plot(bin.er/self.reff, (bin.vrot+bin.vrot_sys)/193d, linestyle=6, symbol='plus', xrange=[0,8],yrange=[-0.5,2],title='N821 Discrete',xtitle='$R/R_e$',ytitle='$\itV/\sigma_e$')

;	filename = 'skim_kinemetry_'+self.name+'.fits'
;	mwrfits, bin, filename, /create

  endfor; end bootstrap loop
  
  if n_elements(nodiscrete) eq 1 then dstellar=d1
  if n_elements(nostellar) eq 1 then ddiscrete=d1
  
  if n_elements(nodiscrete) eq 1 then begin
    innerring = dstellar[where(dstellar.er GT 0.75 AND dstellar.er LT 1.25)]
    innerVEL = median(innerring.v)
    innerPA = median(innerring.pa)
  endif
  
  if n_elements(nostellar) eq 1 then begin
    outerring = ddiscrete[where(ddiscrete.er GT 3.75 AND ddiscrete.er LT 4.25)]
    outerVEL = median(outerring.v)
    outerPA = median(outerring.pa)
  endif
  
  
  if n_elements(nostellar) then begin
  print,strcompress('Delta ='+string(outerVEL - innerVEL))
  print,strcompress('Gamma ='+string(abs(outerPA - innerPA))+'[degrees]')
  endif
  
;  shadyplot,bin,reff,sige, pp1=pp1, thisdevice=thisdevice,nodiscrete=nodiscrete,nostellar=nostellar, $
;  dstellar=dstellar, ddiscrete=ddiscrete, fullrange=fullrange

;  if n_elements(nostellar) eq 1 then p=plot(ddiscrete.er,ddiscrete.pa,linestyle=6,symbol='plus',yrange=[-140,200])
;p=plot(bin.er/self.reff,(bin.vrot+bin.vrot_sys)/sige,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='dodger blue',xrange=[0,8],yrange=[-0.5,3.5])


  profileplot,bin,reff,sige, pp1=pp1, thisdevice=thisdevice,nodiscrete=nodiscrete,nostellar=nostellar, $
  dstellar=dstellar, ddiscrete=ddiscrete, fullrange=fullrange, /dopa
  
  
  
end

;-----------------------------------------------------------------------------------------------

pro map_galaxy::calc_rms, d
	compile_opt idl2, hidden
	if (n_elements(d) gt 0) then begin
		print, 'fix errveldisp'
		if ~tag_exist(d, 'ERRVEL') then d = struct_addtags(d, replicate({errvel:20d},n_elements(d)))
		if ~tag_exist(d, 'ERRVELDISP') then d = struct_addtags(d, replicate({errveldisp:20d},n_elements(d)))
		ibad = where(~finite(d.errvel),nbad)
		if (nbad gt 0) then d[ibad].errvel = 20d
		ibad = where(~finite(d.errveldisp),nbad)
		if (nbad gt 0) then d[ibad].errveldisp = 20d
		d = struct_addtags(d, arr_struct({rms:sqrt(d.vel^2d + d.veldisp^2d), errrms:d.errvel + d.errveldisp}))
	endif
end

;-----------------------------------------------------------------------------------------------


pro map_galaxy__define
	compile_opt idl2, hidden

	void={map_galaxy,  $
			reff:0d, distance:0d, arc2kpc:0d, q:0d, parot:0d, galra:0d, galdec:0d, $
			vsys:0d, discrete:ptr_new(), stellar:ptr_new(), doreff:0, norm:0d, $
			obj:ptr_new(), plot:ptr_new(), xmid:0d, ymid:0d, image:ptr_new(), $
			pos:dblarr(4), dimensions:dblarr(2), current:ptr_new(), dss:ptr_new(), $
			name:'', maxDist:dblarr(2), maxPoints:0l, min_nrange:dblarr(2), $
			svar_info:ptr_new(), gsm:ptr_new(), variogram_info:ptr_new(), $
			kintype:'', vor:0, varplot:0}

end


