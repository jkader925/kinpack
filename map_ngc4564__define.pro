;-----------------------------------------------------------------------------------------------

pro map_ngc4564::get_variogram_info, medianing=medianing, vs_fitting
	compile_opt idl2, hidden
	d = self->get_smeag(vs_fitting)
	d = struct_trimtags(d, select=['X','Y','VEL','VELDISP','H3','H4','RMS'])
	dummy = execute('d = struct_addtags(d, arr_struct({z:d.'+self.kinType+'}))')
	case self.kinType of 
		'VEL'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:100., ftype:'POLY', axisratio:0.47d, inc:4}
						END
		'VELDISP'	: 	BEGIN
							svar_info = {range:300., ftype:'POLY', axisratio:0.43d, inc:3}
						END
		'RMS'		: 	BEGIN
							svar_info = {range:190., ftype:'POLY', axisratio:0.43d, inc:3}
						END
		'H3'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:180., ftype:'POLY', axisratio:0.43d, inc:3}
						END
		'H4'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:180., ftype:'POLY', axisratio:0.43d, inc:3}
						END
		else		: stop
	endcase	
	if (n_elements(self.maxdist) eq 0) then self.maxDist = svar_info.range*0.5
	svar_info = struct_addtags(svar_info, {maxDist:self.maxDist})
	svar_info = struct_addtags(svar_info, {data:d})
	*(self.svar_info) = svar_info



end

;-----------------------------------------------------------------------------------------------

pro map_ngc4564::get_gal_info,vs_fitting
	compile_opt idl2, hidden
	self.name = 'ngc4564'
	self.distance = 15.9d6;16d6
	self.reff = 20.4d  ;SLUGGS Survey Paper Brodie+2014
	self.parot = 48.5d; NED K_s=50, Brodie+2014=48.5d
	self.q = 0.47
	self.galra = 189.11248d 
	self.galdec = 11.4391667d;11.4393194d
	self.arc2kpc = self.distance/206265d/1000d
	self->get_variogram_info,vs_fitting
end

;-----------------------------------------------------------------------------------------------

;function map_ngc4564::get_GCs, red=red, blue=blue
;	compile_opt idl2, hidden
;;	d = mrdfits('/Users/jaaarnol/kriging_data.fits',1,/silent)								;read in the demo data
;	d = mrdfits('I:/Discrete Tracers/ngc4564/kriging_data.fits',1,/silent)                ;read in the demo data
;	d = struct_addtags(d, arr_struct({vel:d.z}))
;	d = struct_trimtags(d, except='Z')
;	return, d
;end

;-----------------------------------------------------------------------------------------------

function map_ngc4564::get_PNe
	compile_opt idl2, hidden
	
	readcol, 'I:/Discrete Tracers/NGC4564/N4564.PNdata', id, x, y, r, pa, a, a, vel, dv, a, format='d,d,d,d,d,d,d,d,d,d', /silent
	d = arr_struct({ra:x/3600d + self.galra, dec:y/3600d + self.galdec, vel:vel - 1155, errvel:dv})
	d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))	
	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc4564::get_sauron,vs_fitting
	compile_opt idl2, hidden
	d = self->map_galaxy::get_sauron(pxf = 'I:/Discrete Tracers/ngc4564/PXF_bin_MS_NGC4564_r3_idl.fits', $
					ms='I:/Discrete Tracers/ngc4564/MS_NGC4564_r3.fits',vs_fitting)
	return,d
end
	
;-----------------------------------------------------------------------------------------------

function map_ngc4564::get_halliday
   compile_opt idl2, hidden
   ;note errvel, veldisp, & errveldisp are guesses at the tag names... 3/25/2014
;readcol, 'I:/Discrete Tracers/NGC4564/Halliday01-kinem2.dat', b, loc, rr, vel, errvel, veldispp, errveldisp, format='a,a,d,d,d,d,d', /silent
;ind = where(strtrim(loc) eq '4564MN')
;;d_0 = arr_struct({xi:ind*0d, yi:-rr[ind], r:rr[ind], vel:vel[ind] - median(vel[ind])})
;d_0 = arr_struct({xi:ind*0d, yi:-rr[ind], r:-rr[ind], vel:vel[ind], errvel:errvel, veldisp:veldispp,errveldisp:errveldisp })
;xyrot,d_0.xi,d_0.yi,90d - self.parot,x,y
;ind = where(strtrim(loc) eq '4564MJ')
;;d_1 = arr_struct({xi:-rr[ind], yi:ind*0d, r:rr[ind], vel:vel[ind] - median(vel[ind])})
;d_1 = arr_struct({xi:-rr[ind], yi:ind*0d, r:-rr[ind], vel:vel[ind], errvel:errvel, veldisp:veldispp,errveldisp:errveldisp})
;d = struct_append(d_0, d_1)
;;xyrot, d.xi, d.yi, -0d*!dtor, xp, yp
;xyrot, d.xi, d.yi, self.parot, xp, yp
;
;d = struct_addtags(d, arr_struct({x:xp, y:yp}))
;d = struct_addtags(d, arr_struct({ra:d.x/3600d + self.galra, dec:d.y/3600d + self.galdec}))
;d1 = d
;  self->calc_rms, d
;;  ;self->get_Flux, d
;;readcol, 'I:/Discrete Tracers/NGC4564/Halliday_4564.txt', slit, loc, rr, vel, 
;;l_0 = arr_struct({slit:b,yi:-rr,vel:vel-1155.,errvel:errvel,veldisp:veldispp,errveldisp:errveldisp})
;stop
;return,d

readcol,'I:/Brutus_InternalDrive/Discrete Tracers/NGC4564/Halliday_4564_mj.dat', b, loc, r, vel, errvel, veldisp, errveldisp, format='a,a,d,d,d,d,d',/silent
relPA = 48d
d = arr_struct({pa:r*0d + relPA, r:r, dx:r*0d, dy:r*0d, vel:vel - 1155d, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp})
d.dx = d.r*cos((90d - d.pa)*!dtor)
d.dy = d.r*sin((90d - d.pa)*!dtor)

readcol,'I:/Brutus_InternalDrive/Discrete Tracers/NGC4564/Halliday_4564_mn.dat', b, loc, r, vel, errvel, veldisp, errveldisp, format='a,a,d,d,d,d,d',/silent
d1 = arr_struct({pa:r*0d + relPA + 90d, r:r, dx:r*0d, dy:r*0d, vel:vel, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp})
d1.dx = d.r*cos((90d - d1.pa)*!dtor)
d1.dx = d.r*sin((90d - d1.pa)*!dtor)

;d=struct_append(d,d1)

  mdm = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
  mdm = struct_addtags(mdm, struct_trimtags(d, except=['PA','R','DX','DY','SYM']))
  mdm = struct_addtags(mdm, arr_struct(projdist2(mdm.ra, mdm.dec, self.galra, self.galdec, self.parot)))

;  a=mdm[where(mdm.vel GE median(mdm.vel))]
;  b=mdm[where(mdm.vel LT median(mdm.vel))]
;  p=plot(a.x,a.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='red',xrange=[-60,60],yrange=[-60,60])
;  p=plot(b.x,b.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='dodger blue',/overplot)
  
  self->calc_rms, mdm 
  ;self->get_Flux, mdm
  return, mdm

end
;-----------------------------------------------------------------------------------------------

function map_ngc4564::get_smeag, medianing=medianing,vs_fitting
	compile_opt idl2, hidden

file='I:/Discrete Tracers/NGC4564/deimos_ppxf_NGC4564.fits'
d = self->get_skims(file,vs_fitting);,medianing=medianing)


  return, d
  
end

;-----------------------------------------------------------------------------------------------
PRO map_ngc4564::rolling_bins_kinemetry, d, nodiscrete=nodiscrete, nostellar=nostellar, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
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
pro map_ngc4564::get_data, d, nopne=nopne, nogcs=nogcs, nodiscrete=nodiscrete, nostellar=nostellar, $
	nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, red=red, blue=blue, h3=h3, h4=h4, veldisp=veldisp,vs_fitting=vs_fitting
	compile_opt idl2, hidden

	if ((n_elements(nopne) eq 1) and (n_elements(nogcs) eq 1)) then nodiscrete = 1
	if (n_elements(nodiscrete) eq 0) then begin
		if (n_elements(nopne) eq 0) then discrete = struct_append(discrete, self->get_PNe())
;		if (n_elements(nogcs) eq 0) then discrete = struct_append(discrete, self->get_GCs(red=red, blue=blue))
		*(self.discrete) = discrete									;put the discrete data in a pointer
		d=discrete
	endif
	if ((n_elements(nosauron) eq 1) and (n_elements(nolong) eq 1) and (n_elements(nosmeag) eq 1)) then nostellar = 1
	if (n_elements(nostellar) eq 0) then begin
    if (n_elements(nosauron) eq 0) then stellar = self->get_sauron(vs_fitting)
		if (n_elements(nolong) eq 0) then stellar = struct_append(stellar, self->get_halliday())
		if (n_elements(nosmeag) eq 0) then stellar = struct_append(stellar, self->get_smeag(vs_fitting))
    *(self.stellar) = stellar
		if (n_elements(stellar) gt 0) then *(self.stellar) = stellar
		s=stellar
	endif
    if n_elements(nosauron) eq 0 or n_elements(nosmeag) eq 0 or n_elements(nolong) eq 0 or n_elements(nostellar) eq 0 then d = struct_selecttags(s,select=['x','y','r','pa','vel','errvel']) ;subarray to hold stellar+discrete data (start by adding stellar data)
    if n_elements(nodiscrete) eq 0 then begin
      dummy = struct_selecttags(discrete,select=['x','y','r','pa','vel','errvel']) 
      d = struct_append(d,dummy)
    endif
end

;-----------------------------------------------------------------------------------------------
pro map_ngc4564::isophote, p2
 compile_opt idl2, hidden




  s = *(self.plot)
  pos = self.pos
;-----------------------------------------------------------------------------------------------
nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
;q = 0.68;0.8 ; outer isophote flattening
q = 0.67 ;made-up.
;maxa = 2*self.reff/sqrt(1-q)  ; DONT KNOW BUT PRETTY SURE THIS IS ACTUALLY 2 R_e
maxa = (4/7d)*(max(s.xsiz)/2)
;maxa = 1d/sqrt(q)*s.xsiz/float(difference(minmax(s.xrange)))
;maxa = self.reff*2.  ;this is just 2*Re [arcsec]
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
q = 0.47  ;inner isophote flattening
;maxa = 1d/sqrt(q)*s.xsiz/float(difference(minmax(s.xrange)))
;maxa = 0.5*self.reff/sqrt(1-q)
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

  p2 = plot(xhigh1, yhigh1, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0) ;outer?
  p2 = plot(xhigh2, yhigh2, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)  ;inner





















end
;-----------------------------------------------------------------------------------------------

;pro map_ngc4564::isophote, p2
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

pro map_ngc4564__define
	compile_opt idl2, hidden

	void={map_ngc4564, inherits map_galaxy}

end



