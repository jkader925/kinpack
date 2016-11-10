;-----------------------------------------------------------------------------------------------

pro map_ngc1023::get_variogram_info, vs_fitting
	compile_opt idl2, hidden
	d = self->get_smeag(vs_fitting)
	d = struct_trimtags(d, select=['X','Y','VEL','VELDISP','H3','H4','RMS'])
	dummy = execute('d = struct_addtags(d, arr_struct({z:d.'+self.kinType+'}))')	
	case self.kinType of 
		'VEL'		: 	BEGIN
							svar_info = {range:240., ftype:'GAUSSIAN', axisratio:0.41d, inc:50}
						END
		'VELDISP'	: 	BEGIN
							svar_info = {range:180., ftype:'POLY', axisratio:1.0d, inc:5}
						END
		'RMS'		: 	BEGIN
							svar_info = {range:180., ftype:'POLY', axisratio:1.0d, inc:5}
						END
		'H3'		: 	BEGIN
							self.varplot = 1
							svar_info = {range:120., ftype:'POLY', axisratio:1.0d, inc:4}
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

pro map_ngc1023::get_gal_info,vs_fitting
	compile_opt idl2, hidden
	self.name = 'NGC1023'
	self.reff = 48.
	self.distance = 11.8d6
	self.kintype = 'VEL' ;this wasnt here before -jk
	self.q = 0.41 ;1 R_e...              0.67 (old weird value)
	self.parot = 87d
	self.galra = 40.1000421d  
	self.galdec = 39.0632850d
	self.arc2kpc = self.distance/206265d/1000d
	self->get_variogram_info,vs_fitting
end

;-----------------------------------------------------------------------------------------------

function map_ngc1023::get_smeag, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
;	d = self->get_skims('I:/Discrete Tracers/NGC1023/NGC1023_skims.fits', medianing=medianing)
  d = self->get_skims('I:/Discrete Tracers/NGC1023/deimos_ppxf_NGC1023.fits',vs_fitting);, medianing=medianing)
d = d[where(d.vel gt -300)]	

a=d[where(d.vel lt median(d.vel))]
b=d[where(d.vel ge median(d.vel))]
;p=plot(a.x,a.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='red',title='SMEAGOL')
;p=plot(b.x,b.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='dodger blue',/overplot)
	return, d
end

;-----------------------------------------------------------------------------------------------

function map_ngc1023::get_debattista_mdm
  compile_opt idl2, hidden
  readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC1023/n1023_debattistalong.txt', r, pa, vel, errvel, veldisp, errveldisp
  relPA = 87d - self.parot
  d = arr_struct({r:r, pa:pa*0d, dx:r*0d, dy:r*0d, vel:vel-605d, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp})
  d.pa = 0d;relPA
;  d.dx = d.r*cos((d.pa + 90d - relPA)*!dtor)
;  d.dy = d.r*sin((d.pa + 90d - relPA)*!dtor)
  d.dx = d.r*cos((d.pa)*!dtor)
  d.dy = d.r*sin((d.pa)*!dtor)
;  xyrot,d.dx,d.dy,(90d + 80d)*:!dtor,x,y
;  d=struct_addtags(d,arr_struct({xx:x,yy:y}))

;  d.dx = -d.dx
;  d.dy = -d.dy

;  ind = where(d.pa ne 0)
;  d1 = d[where(d.pa ne 0)]
;  d1.dx = -d1.dx
;  d1.dy = -d1.dy
;  d = struct_append(d, d1)  

  readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC1023/Debdata1a_north12.txt', r, pa, vel, errvel, veldisp, errveldisp
  d2 = arr_struct({r:r, pa:pa*0d, dx:r*0d, dy:r*0d, vel:vel-605d, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp})
  relPA = 87d - self.parot
  d2.dx = d2.r*cos((d2.pa)*!dtor)
  d2.dy = (d2.r*sin((d2.pa)*!dtor))+12d
  ind = where(d2.pa ne 0)
  d3 = d2[where(d2.pa ne 0)]
  d3.dx = -d3.dx
  d3.dy = -d3.dy

  d = struct_append(d, d2)

  readcol, 'I:/Brutus_InternalDrive/Discrete Tracers/NGC1023/Debdata1a_south16.txt', r, pa, vel, errvel, veldisp, errveldisp
  d4 = arr_struct({r:r, pa:pa*0d, dx:r*0d, dy:r*0d, vel:vel-605d, errvel:errvel, veldisp:veldisp, errveldisp:errveldisp})
  relPA = 87d - self.parot
  d4.dx = d4.r*cos((d4.pa)*!dtor)
  d4.dy = (d4.r*sin((d4.pa)*!dtor))-16d
  ind = where(d4.pa ne 0)
  d5 = d4[where(d4.pa ne 0)]
  d5.dx = -d5.dx
  d5.dy = -d5.dy
  
  d = struct_append(d, d4)

  het = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
  het = struct_addtags(het, struct_trimtags(d, except=['PA','R','DX','DY']))
  het = struct_addtags(het, arr_struct(projdist2(het.ra, het.dec, self.galra, self.galdec, self.parot)))
  het.x = -het.x
  het=struct_trimtags(het,except=['PA'])
  het = struct_addtags(het, arr_struct({pa:pacalc(het.x, het.y, self.parot,/radians)}))
  
  self->calc_rms, het
;  self->get_Flux, het
;xyrot,d.dx,d.dy,(90d + 80d)*:!dtor,x,y
;  d=struct_addtags(d,arr_struct({xx:x,yy:y}))

  return, het
  
  
  ;  mdm = arr_struct({ra:d.dx/3600d + self.galra, dec:d.dy/3600d + self.galdec})
;  mdm = struct_addtags(mdm, struct_trimtags(d, except=['PA','R','SYM']))
;  mdm = struct_addtags(mdm, arr_struct({x:mdm.dx, y:mdm.dy}))
;  stop
;  mdm = struct_addtags(mdm, arr_struct(projdist2(mdm.ra, mdm.dec, self.galra, self.galdec, self.parot)))
;  self->calc_rms, mdm 
;  self->get_Flux, mdm
;  stop
;  return, mdm


end

;-----------------------------------------------------------------------------------------------


function map_ngc1023::get_GCs, red=red, blue=blue
  compile_opt idl2, hidden  
    readcol, 'I:/Discrete Tracers/NGC1023/N1023gcvels2.txt', id, x, y, a, a, vel, dvel, a, a, a, a, /silent
  d = arr_struct({ra:x/3600 + self.galra, dec:y/3600 + self.galdec, vel:vel-602.7, errvel:dvel})

  d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, self.parot)))  
  if (n_elements(red) gt 0) then d = d[where((g - i) gt 0.88)] else if (n_elements(blue) gt 0) then d = d[where((g - i) le 0.88)]
;  a=d[where(d.vel ge median(d.vel))]
;  b=d[where(d.vel lt median(d.vel))]
;  p=plot(a.x,a.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='red',title='N1023 GCs')
;  p=plot(b.x,b.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='dodger blue',/overplot)
  return, d
end

;-----------------------------------------------------------------------------------------------
function map_ngc1023::get_sauron, medianing=medianing,vs_fitting
	compile_opt idl2, hidden
;	d = self->map_galaxy::get_sauron(pxf = '/Users/jaaarnol/sauron/PXF_NGC1023_r2/PXF_bin_MS_NGC1023_r2.fits', $
;									ms='/Users/jaaarnol/sauron/PXF_NGC1023_r2/MS_NGC1023_r2.fits')
	d = self->map_galaxy::get_sauron(pxf = 'I:/Discrete Tracers/NGC1023/PXF_bin_MS_NGC1023_r2_idl.fits', $
                  ms='I:/Discrete Tracers/NGC1023/MS_NGC1023_r2.fits',vs_fitting);, medianing=medianing)
d = d[where(d.vel gt -300)]	
	return, d
end

;-----------------------------------------------------------------------------------------------
function map_ngc1023::get_PNe
  compile_opt idl2, hidden
  
;readcol, 'I:/Discrete Tracers/NGC1023/N1023.PNdata2a', r, x, y, pa, a, vel, dv, b , c, format='d,d,d,d,d,d,d,d,d,d', /silent
readcol, 'I:/Discrete Tracers/NGC1023/N1023.PNdata2a', x, y, r, pa, vel, dv, a , b, format='d,d,d,d,d,d,d,d', /silent
  d = arr_struct({xi:x, yi:y, vel:vel, errvel:dv, pa:pa, a:a, b:b})


  



  xyrot, d.xi, d.yi, -7d*!dtor, xp, yp
  d = struct_addtags(d, arr_struct({ra:xp, dec:yp}))
  d = arr_struct({ra:x/3600d + self.galra, dec:y/3600d + self.galdec, vel:vel-602.7, errvel:dv})
  d = struct_addtags(d, arr_struct(projdist2(d.ra, d.dec, self.galra, self.galdec, -self.parot)))
  d=struct_trimtags(d,except=['PA'])
  d = struct_addtags(d, arr_struct({pa:pacalc(d.x, d.y, -self.parot,/radians)}))
  
;  a = d[where(d.vel ge median(d.vel))]
;  b = d[where(d.vel lt median(d.vel))]
;  p=plot(a.x,a.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='red',title='N1023 PNe')
;  p=plot(b.x,b.y,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='dodger blue',/overplot)
  
  
  return,d
end

;-----------------------------------------------------------------------------------------------
PRO map_ngc1023::rolling_bins_kinemetry, d, nodiscrete=nodiscrete, nostellar=nostellar, innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
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

pro map_ngc1023::get_data, d, nopne=nopne, nogcs=nogcs, nodiscrete=nodiscrete, nostellar=nostellar, $
  nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, red=red, blue=blue, green=green, medianing=medianing, $
  vs_fitting=vs_fitting

  compile_opt idl2, hidden
  if ((n_elements(nopne) eq 1) and (n_elements(nogcs) eq 1)) then nodiscrete = 1
  if (n_elements(nodiscrete) eq 0) then begin
    if (n_elements(nopne) eq 0) then discrete = self->get_PNe()               ;retrieve the PN data
    if (n_elements(nogcs) eq 0) then discrete = struct_append(discrete, self->get_GCs(red=red, blue=blue))  ;retrieve the GC data
    *(self.discrete) = discrete                               ;put the discrete data in a pointer
    dis = discrete
  endif

  if ((n_elements(nosauron) eq 1) and (n_elements(nolong) eq 1) and (n_elements(nosmeag) eq 1)) then nostellar = 1
  if (n_elements(nostellar) eq 0) then begin
    if (n_elements(nosauron) eq 0) then stellar = self->get_sauron(vs_fitting);medianing=medianing)
    if (n_elements(nolong) eq 0) then begin
;     ls = struct_append(ls,  self->get_pinkney_mdm())
      ls = struct_append(ls, self->get_debattista_mdm())
      stellar = struct_append(stellar, ls)
    endif

    if (n_elements(nosmeag) eq 0) then stellar = struct_append(stellar, self->get_smeag(vs_fitting));medianing=medianing))
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
;-----------------------------------------------------------------------------------------------

pro map_ngc1023::isophote, p2
	compile_opt idl2, hidden

  s = *(self.plot)
  pos = self.pos

nbound = 100
xmid = ((s.xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((s.ysiz-1)/2.)                        ;y-midpoint of the map
q = 0.36 ; outer isophote flattening
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
q = 0.41  ;inner isophote flattening
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

pro map_ngc1023__define
	compile_opt idl2, hidden

	void={map_ngc1023, inherits map_galaxy}

end






 
 
