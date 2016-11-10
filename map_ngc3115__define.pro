;-----------------------------------------------------------------------------------------------

pro map_ngc3115::get_variogram_info,vs_fitting
	compile_opt idl2, hidden
	d = self->get_smeag(vs_fitting)
	d = struct_trimtags(d, select=['X','Y','VEL','VELDISP','H3','H4','RMS'])
	dummy = execute('d = struct_addtags(d, arr_struct({z:d.'+self.kinType+'}))')
	case self.kinType of 
		'VEL'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:300., ftype:'POLY', axisratio:0.5d, inc:3}
						END
		'VELDISP'	: 	BEGIN
							svar_info = {range:300., ftype:'POLY', axisratio:0.5d, inc:3}
						END
		'RMS'		: 	BEGIN
							svar_info = {range:190., ftype:'POLY', axisratio:0.5d, inc:3}
						END
		'H3'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:180., ftype:'POLY', axisratio:0.5d, inc:3}
						END
		'H4'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:180., ftype:'POLY', axisratio:0.5d, inc:3}
						END
		else		: stop
	endcase	
	stop
	if (n_elements(self.maxdist) eq 0) then self.maxDist = svar_info.range*0.5
	svar_info = struct_addtags(svar_info, {maxDist:self.maxDist})
	svar_info = struct_addtags(svar_info, {data:d})
	*(self.svar_info) = svar_info
end

;-----------------------------------------------------------------------------------------------

pro map_ngc3115::get_gal_info,vs_fitting
	compile_opt idl2, hidden
	self.name = 'NGC3115'
	self.distance = 10d6
	self.reff = 57.
	self.parot = 44d
	self.q = 0.29
	self.galra = 151.308250d  
	self.galdec = -7.718583d
	self.arc2kpc = self.distance/206265d/1000d
	self->get_variogram_info,vs_fitting
end

;-----------------------------------------------------------------------------------------------

function map_ngc3115::get_GCs, red=red, blue=blue
	compile_opt idl2, hidden
;	d = mrdfits('/Users/jaaarnol/kriging_data.fits',1,/silent)								;read in the demo data
	d = mrdfits('I:/Discrete Tracers/NGC3115/kriging_data.fits',1,/silent)                ;read in the demo data
	d = struct_addtags(d, arr_struct({vel:d.z}))
	d = struct_trimtags(d, except='Z')
  d = struct_addtags(d, arr_struct({pa:pacalc(d.x, d.y, self.parot,/radians)}))
  stop
	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc3115::get_PNe
	compile_opt idl2, hidden
  readcol, 'I:/Discrete Tracers/NGC3115/n3115_pndata.txt', ra, dec, vel, velerr, format='d,d,d,d', /silent
  d = arr_struct({ra:ra, dec:dec, vel:vel - 699.8, errvel:velerr})
  d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))  

	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc3115::get_Norris
	compile_opt idl2, hidden
;	d = mrdfits('/Users/jaaarnol/norris3115Kinematics.fits',1,/silent)
  d = mrdfits('I:/Discrete Tracers/NGC3115/norris3115Kinematics.fits',1,/silent)
	xyrot, d.x, d.y, (44d - 90d)*!dtor, xp, yp			;35d - self.parot
	d = struct_trimtags(d, except=['X','Y'])
	d = struct_addtags(d, arr_struct({ra:self.galra - xp/3600d, dec:self.galdec + yp/3600d}))
	d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))
	d.vel -= median(d.vel,/even)		
	self->calc_rms, d
	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc3115::get_smeag, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
;	d = self->get_skims('/Volumes/data4/skimspaper/NGC3115/NGC3115_skims.fits')
;	d = self->get_skims('/Volumes/data4/skimspaper/NGC3115/comb_deimos_ppxf_v0.fits')
;  d = self->get_skims('I:/Discrete Tracers/NGC3115/comb_deimos_ppxf_NGC3115.fits', medianing=medianing)
  d = self->get_skims('I:/Discrete Tracers/NGC3115/deimos_ppxf_NGC3115.fits',vs_fitting);, medianing=medianing)
	if ~tag_exist(d, 'ERRVEL') then d = struct_addtags(d, replicate({errvel:35d},n_elements(d)))
	if ~tag_exist(d, 'ERRVELDISP') then d = struct_addtags(d, replicate({errveldisp:35d},n_elements(d)))
	if ~tag_exist(d, 'ERRH3') then d = struct_addtags(d, replicate({errh3:0.05d},n_elements(d)))
	if ~tag_exist(d, 'ERRH4') then d = struct_addtags(d, replicate({errh4:0.05d},n_elements(d)))
;d.vel -= 700.	
	;d.vel -= median(d.vel,/even)
	d = d[where(d.y gt -300)]
	d = d[where(abs(d.h4) lt 0.2)]
;	d = d[where(d.sn gt 50)]
;stop	
	;d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))
	return, d
end
;-----------------------------------------------------------------------------------------------
PRO map_ngc3115::rolling_bins_kinemetry, d, nodiscrete=nodiscrete, nostellar=nostellar, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
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
;-----------------------------------------------------------------------------------------------


pro map_ngc3115::get_data, d, nopne=nopne, nogcs=nogcs, nodiscrete=nodiscrete, nostellar=nostellar, $
	nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, red=red, blue=blue, h3=h3, h4=h4, veldisp=veldisp, medianing=medianing,vs_fitting=vs_fitting
	compile_opt idl2, hidden

	if ((n_elements(nopne) eq 1) and (n_elements(nogcs) eq 1)) then nodiscrete = 1
	if (n_elements(nodiscrete) eq 0) then begin
		if (n_elements(nopne) eq 0) then discrete = struct_append(discrete, self->get_pne())
		if (n_elements(nogcs) eq 0) then discrete = struct_append(discrete, self->get_GCs(red=red, blue=blue))
		*(self.discrete) = discrete									;put the discrete data in a pointer
	endif
	
	if ((n_elements(nosauron) eq 1) and (n_elements(nolong) eq 1) and (n_elements(nosmeag) eq 1)) then nostellar = 1
	if (n_elements(nostellar) eq 0) then begin
		if (n_elements(nolong) eq 0) then stellar = struct_append(stellar, self->get_Norris())
		if (n_elements(nosmeag) eq 0) then stellar = struct_append(stellar, self->get_smeag(vs_fitting));medianing=medianing))
s = stellar
s = s[where(abs(s.h4) le 0.2)]
stellar = s

		if (n_elements(stellar) gt 0) then *(self.stellar) = stellar
	endif

    if n_elements(nosauron) eq 0 or n_elements(nosmeag) eq 0 or n_elements(nolong) eq 0 or n_elements(nostellar) eq 0 then d = struct_selecttags(s,select=['x','y','r','pa','vel','errvel']) ;subarray to hold stellar+discrete data (start by adding stellar data)
    if n_elements(nodiscrete) eq 0 then begin
      dummy = struct_selecttags(discrete,select=['x','y','r','pa','vel','errvel']) 
      d = struct_append(d,dummy)
    endif
end

;-----------------------------------------------------------------------------------------------

;pro map_ngc3115::isophote, p2
;	compile_opt idl2, hidden
;
;	s = *(self.plot)
;	pos = self.pos
;	return
;	
;	files = '/Users/jaaarnol/n3377_'+roundx(s.isophote_level,0)+'Re.dat'
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
;		con.x = smooth(con.x*s.xsiz/float(difference(minmax(s.xrange))) + self.xmid, smcon, /edge)
;		con.y = smooth(con.y*s.xsiz/float(difference(minmax(s.xrange))) + self.ymid, smcon, /edge)
;		con = struct_append(con, con[0])
;		p2 = plot(con.x, con.y, position=pos, /overplot, thick=s.c_thick)
;		++k
;	endforeach
;end

;-----------------------------------------------------------------------------------------------
pro map_ngc3115::isophote, p2
  compile_opt idl2, hidden

  s = *(self.plot)
  pos = self.pos
  
  
nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
q = 0.36 ; 2Re
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
q = 0.29 ; 1Re
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
pro map_ngc3115__define
	compile_opt idl2, hidden

	void={map_ngc3115, inherits map_galaxy}

end






 
 
