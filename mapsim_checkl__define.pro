;+
;PURPOSE
; To process simulation data into kinematic maps.
;
;TOOLS
; INIT, data, [file]
;   Read in simulation data and initialize default options
;   -Inputs:
;     data = simulation data in an array of structures with tags:
;       X, Y, Z, VX, VY, and VZ. Any additional tags are available in the
;       cull_data method, and will be averaged in the final voronoi bins.
;   -Keywords:
;     file = name of simulation data file for output naming purposes
; CULL_SATELLITES, [binsize=binsize]
;   Attempts to find and clip out satellite galaxies in cosmological simulations.
;   -Keywords:
;     binsize = binsize of the histogram of radii where satellites are found (default = 0.25 kpc)
; REALIGN_ANGULAR_MOMENTUM, [rot=rot, am_rad=am_rad]
;   Measures the direction of the net angular momentum of the inner galaxy and realigns
;   the galaxy to put this on the z-axis.
;   -Keywords:
;     rot = input 3x3 rotation matrix to rotate the coordinates by.
;     am_rad = radius within which the net angular momentum is calculated (default = 2 kpc)
; ROTATE_AROUND_AXIS, [xrot, yrot, zrot]
;   Rotate around any given axis by the specified amount
;   -Keywords:
;     XROT = angle, in degrees, to rotate around x-axis
;     YROT = angle, in degrees, to rotate around y-axis
;     ZROT = angle, in degrees, to rotate around z-axis
; CULL_DATA, [str]
;   CUll any data points via a string, which is executed in a where statement.
;   e.g., CULL_DATA, '(data.r lt 30) and (data.age lt 10)
;   Any tag name put in the input data structure is available.
; VORONOI_BIN, [restore=restore, quick=quick, targetN=targetN, pixelSize=pixelSize, eps=eps, nthreads=nthreads]
;   -Keywords:
;     restore = input save file to use instead of creating new voronoi bins
;     quick : set to use a default voronoi bin configuration
;     targetN = target number of objects in each voronoi bin (VORONOI_2D_BINNINGmod)
;     pixelsize, eps = see VORONOI_2D_BINNINGmod.pro
;     nthreads = number of sub-processes to use when accreting points into voronoi bins
;
;EXAMPLE
; function mapsim_checkL_sim_get_data, file
;   nRows = File_lines(file)
;   d = replicate({x:0d, y:0d, z:0d, VX:0d, VY:0d, VZ:0d}, nRows)
;   openr, lun, file, /get_lun
;   readf, lun, d
;   free_lun, lun
;   return, struct_addtags(d, arr_struct({r:sqrt(d.x^2d + d.y^2d + d.z^2d), mass:d.x*0 + 1.}))
; end
;
; pro mapsim_checkL
;   file = 'simulation.dat'         ;simulation data file
;   data = mapsim_checkL_sim_get_data(file)    ;get the simulation data
;   mapsim_checkL = obj_new('mapsim_checkL', data, file)  ;instantiate the mapsim_checkL object
;   mapsim_checkL->cull_satellites         ;cull the satellites
;   mapsim_checkL->realign_angular_momentum    ;realign the angular momentum with the z-axis
;   mapsim_checkL->rotate_around_axis        ;rotate around an axis
;   mapsim_checkL->cull_data, str=str        ;cull any data
;   mapsim_checkL->voronoi_bin           ;create and accrete into voronoi bins
;   mapsim_checkL->save              ;save the results
; end
;
;Written by Jacob A. Arnold, UCSC, January 6, 2012
;-

;-----------------------------------------------------------------------------------------------

FUNCTION mapsim_checkL::init, file, data=data
  compile_opt idl2, hidden

  if (n_elements(file) gt 0) then self.file = file else stop  ;self.file = '/NULL'
  self.am_rad = 2
  self.targetN = 400.
  self.pixelsize = 0.1
  self.eps = 0.1d
  self.nthreads = 6
  self.binsize = 0.25

  self.obj = ptr_new(/allocate_heap)
  self.data = ptr_new(/allocate_heap)
  self.mainGalaxy = ptr_new(/allocate_heap)
  ;self.rot = diag_matrix([1d,1d,1d])
  ;self.axisrot = diag_matrix([1d,1d,1d])

  if (n_elements(data) eq 0) then *(self.data) = self->getData(file) $
  else *(self.data) = data

  if 0*(n_elements(data) gt 0) then begin
    d = data
    ;req_tags = ['X','Y','Z','VX','VY','VZ']
    ;foreach tag, req_tags do if ~tag_exist(d, tag) then stop
    if ~tag_exist(d, 'R') then d = struct_addtags(d, {r:sqrt(d.x^2d + d.y^2d + d.z^2d)})
    stop
    if ~tag_exist(d, 'MASS') then d = struct_addtags(d, {mass:d.x*0 + 1d})
    stop
  endif
  ;*(self.data) = d

  return, 1
END

;-----------------------------------------------------------------------------------------------

PRO mapsim_checkL::cleanup
  compile_opt idl2, hidden
  ptr_free, self.data
  ptr_free, self.obj
  ptr_free, self.mainGalaxy
END

;-----------------------------------------------------------------------------------------------

FUNCTION mapsim_checkL::getdata, file
  compile_opt idl2, hidden
  outFits = (strsplit(file,'/',/extract))[-1]+'.fits'
  fitsfile = file_search(outFits, count=nfits)
  if ~nFits then begin
    nRows = File_lines(file)
    case (strsplit(file,'.',/extract))[-1] of
      'dat' : begin
        hdr = replicate({hdr:''},8)
        d = replicate({num:0l, x:0d, y:0d, z:0d, VX:0d, VY:0d, VZ:0d, mass:0d, id:0l}, nRows - 8 - 2)
        openr, lun, file, /get_lun
        readf, lun, hdr
        readf, lun, d
        free_lun, lun
      end
      'dat0'  : begin
        hdr = replicate({hdr:''},6)
        d = replicate({id:0l, x:0d, y:0d, z:0d, VX:0d, VY:0d, VZ:0d, mass:0d, age:0d, feh:0d, afe:0d}, nRows - 7)
        openr, lun, file, /get_lun
        readf, lun, hdr
        readf, lun, d
        free_lun, lun
      end
      'dat1'  : begin
        d = replicate({id:0l, x:0d, y:0d, z:0d, Rx:0d, PA:0d, VX:0d, VY:0d, VZ:0d, mass:0d, age:0d, feh:0d, afe:0d, lr:0d}, nRows)
        openr, lun, file, /get_lun
        readf, lun, hdr
        readf, lun, d
        free_lun, lun
      end
      'dat2' : begin
        hdr = replicate({hdr:''},8)
        d = replicate({id:0l, x:0d, y:0d, z:0d, R:0d, PA:0d, VX:0d, VY:0d, VZ:0d, mass:0d, rr:0d, Rm:0l}, nRows - 8 - 2)
        openr, lun, file, /get_lun
        readf, lun, hdr
        readf, lun, d
        free_lun, lun
      end
      'ascii' : begin
        hdr = {num:0l, dimen:0d, time:0d}
        d0 = replicate({mass:0d},(nRows - 1l)/3.)
        d = replicate({x:0d, y:0d, z:0d}, (nRows - 1l)/3.)
        d1 = replicate({VX:0d, VY:0d, VZ:0d}, (nRows - 1l)/3.)
        openr, lun, file, /get_lun
        readf, lun, hdr
        readf, lun, d0
        readf, lun, d
        readf, lun, d1
        free_lun, lun
        d0.mass *= 10^10d
        d = struct_addtags(d, d1)
        d = struct_addtags(d, d0)
        d = struct_addtags(d, arr_struct({id:lindgen(n_Elements(d)), r:sqrt(d.x^2d + d.y^2d + d.z^2d)}))
      end
      else  : stop
    endcase
    ;   d = struct_trimtags(d, select=['X','Y','Z','VX','VY','VZ','MASS','FEH'])
    d = struct_trimtags(d, select=['num','X','Y','Z','VX','VY','VZ','MASS','ID'])

    ;d = struct_addtags(d, arr_struct({mass:d.x*0 + 1.}))
    d = struct_addtags(d, arr_struct({r:sqrt(d.x^2d + d.y^2d + d.z^2d)}))
    ;    d = struct_addtags(d, arr_struct({r:sqrt(d.x^2d + d.y^2d + d.z^2d)}))
    ;mwrfits, d, outFits, /create
  endif else d = mrdfits(outFits, 1, /silent)
  return, d

END

;-----------------------------------------------------------------------------------------------

PRO mapsim_checkL::cull_satellites, binsize=binsize
  compile_opt idl2, hidden

  if (n_elements(binsize) gt 0) then self.binsize = binsize
  d = *(self.data)
  k = 0l
  ;  locs = [100,50,40,30,20,10,9,8,7,6,5,4,3]
  ;  locs = [10,5,3]
  locs = reverse((indgen(19)+2)*5)
  ; window, 1, xsize=600,ysize=600
  ; window, 0, xsize=600,ysize=600

  for j=0,(n_elements(locs)-1) do begin
    while 1 do begin
      h = histogram(d.r, binsize=self.binsize, locations=loc, reverse_indices=R)

      h = h - smooth(h,51,/edge_truncate)
      ;     if ((k eq 0) and (j eq 0)) then window, 2, xsize=600,ysize=600 else wset, 2
      ;     plot, loc, h, yrange=[0,100]

      hold = where(loc gt locs[j])                          ;radius bins that are farther out than locs[j]
      maxh = max(h[hold],imax)                            ;find the highest bin
      ;     oplot, intarr(10)+loc[hold[imax]],findgen(10)*1d6,color=fsc_color('red')    ;mark this bin
      dsat = d[where( abs(d.r - loc[hold[imax]]) lt 2.0 )]                ;isolate the particles with 2 kpc of this radius
      ;      dsat = d[where( abs(d.r - loc[hold[imax]]) lt 0.25 )]                ;isolate the particles with 0.25 kpc of this radius




      xh = histogram(dsat.x,binsize=self.binsize,locations=xloc)            ;make a histogram of x-coords for the satellite
      maxh = max(xh,imax)
      dsat1 = dsat[where( abs(dsat.x - xloc[imax]) lt 2.0,ndsat)]
      ;            dsat1 = dsat[where( abs(dsat.x - xloc[imax]) lt 2 ,ndsat)]

      ;     if ((k eq 0) and (j eq 0)) then window, 0, xsize=600,ysize=600 else wset, 0
      ;      if (ndsat lt 10) then break
      if (ndsat lt 2) then break
      ;     plothist, dsat.y, bin=0.01

      sat = [median(dsat.x),median(dsat.y),median(dsat.z)]
      data = d[where(sqrt((d.x - sat[0])^2d + (d.y - sat[1])^2d + (d.z - sat[2])^2d) lt 2.0)]
      if (n_elements(data) lt 10) then break

      crad = 10.*robust_sigma(sqrt((data.x - sat[0])^2d + (data.y - sat[1])^2d + (data.z - sat[2])^2d))   ;boundary of clump is defined as the radius where the point spread dispersion is ten times the central clump value.

      ;      crad = 1.*robust_sigma(sqrt((data.x - sat[0])^2d + (data.y - sat[1])^2d + (data.z - sat[2])^2d))   ;boundary of clump is defined as the radius where the point spread dispersion is ten times the central clump value.

      isat = where(total(d[where(d.r lt (rad = sqrt(sat[0]^2d + sat[1]^2d + sat[2]^2d)))].mass)/rad^2d lt $
        total(d[(isat1 = where((srad = sqrt((d.x - sat[0])^2d + (d.y - sat[1])^2d + (d.z - sat[2])^2d)) lt crad,complement=ikeep1))].mass)/srad^2d,nsat,complement=ikeep)

      ;        isat = d[where((srad = sqrt((d.x - sat[0])^2d + (d.y - sat[1])^2d + (d.z - sat[2])^2d)) lt crad, nsat,complement=ikeep)]
      ;the clump qualifies as a satellite if its mass per unit radius is greater than that of the galaxy, out to the position of the clump. Mass per unit radius of the
      ;clump is calculated within 'crad' of the clump center.

      ;maybe if the clump boundary is extended to something like 20 times the central dispersion, it will relax the condition on what qualifies as a satellite. Maybe the
      ;inner satellites are not qualifying because the mass-per-unit-radius of the galaxy at small radii is very large!

      if (nsat lt 2) then break

      ;     if ((k eq 0) and (j eq 0)) then p=plot(d.x,d.y,linestyle=6,symbol='dot',title='before',layout=[2,1,1])
      ;     p=plot(d[ikeep].x,d[ikeep].y,linestyle=6,symbol='dot',title='after',layout=[2,1,2],/CURRENT)
      ;      d1 = isat


      d1 = d[isat]
      d = d[ikeep]
      ;     if ((k eq 0) and (j eq 0)) then window, 1, xsize=600,ysize=600 else wset, 1
      ;     plot, d1.x, d1.y, psym=3, /iso
      ;      print, strcompress('satellite size ='+string(crad)+'kpc '+' nparticles:'+string(n_elements(d1)))
      ;      wset, 0
      ++k
      ;     pause

    endwhile
  endfor

  h = histogram(d.r, binsize=self.binsize, locations=loc, reverse_indices=R)
  h = h - smooth(h,11,/edge_truncate)
  ; if ((k eq 0) and (j eq 0)) then window, 2 else wset, 2
  ; plot, loc, h, yrange=[0,100]

  *(self.data) = d

  ;  window,0,xsize=1000,ysize=1000
  ;  plot,d.x,d.y,psym=3
  ;stop
END

;-----------------------------------------------------------------------------------------------
PRO mapsim_checkL::clump_destroyer, data, reff, sige, q, cbinsize=cbinsize, diagnostics=diagnostics
  compile_opt idl2, hidden

data = *(self.data)
odat = data
if diagnostics eq 1 then window,2,xsize=500,ysize=500

odat = odat[sort(odat.y)]
binsize = 1d4
ndbins=round((n_elements(odat.y))/binsize) ;number of "depth dimension" bins which have 1d4 stars (140 bins)



for ii = 0, ndbins-2 do begin
  
  data = odat[ii*binsize:(ii+1)*binsize]
  
  
  if diagnostics eq 1 then wset,2
  if diagnostics eq 1 then begin
    if ii eq 0 then plot,odat.y,odat.z,psym=3,xrange=[-30,30],yrange=[-30,30] 
    oplot,data.y,data.z,psym=3,color=fsc_color('green')
  endif
    ;----------Calculate elliptical radii of data, given an ellipticity and position angle-------\
    nbins = 12
    binsize = binsize
    kaxisratio = 1.0 ; Shape of rolling bins
    kpa = 0.*!dtor ; Position angle of elliptical rolling bins and value of fixed kPA
    data = struct_addtags(data, arr_struct({pa:pacalc(data.x, data.y, 0d)*!dtor}))  ;calculate position angle
    dummy = ellippar(data.x, data.y, kaxisratio, kpa, struc=data, /sort) ;Install new a and b values for each object. Also sort by elliptical radius.
    ;--------------------------------------------------------------------------------------------

    bins = replicate({binnum:0d,rad:0d,PA:0d}, nbins)
    locs = [5,4,3,2,1,0]
    i=0
    j=0
    for i=1, nbins-1 do begin ;for each of the 12 annular radius bins
     for j=0,n_elements(locs)-1 do begin ;step backwards in PA from [286deg - 0deg] to [230deg - 0deg] ... [57deg - 0deg] (save time by not doing this?)
      k=0
      while 1 do begin ;isolate and destroy clumps in PA range (PA < locs[j]), until largest clump has density less than 7 sigma from central density
        if diagnostics eq 1 then wset,0
        d = data[where(data.er GE i*reff and data.er LT (i+1)*reff)] ;annulus data
        if n_elements(d) LT 100 then break
        if diagnostics GT 0 and j eq 0 then window,0,xsize=500,ysize=500 else wset,0
        if diagnostics GT 0 and j eq 0 then plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3  
        if diagnostics GT 0 then oplot,d.x,d.y,psym=3,color=fsc_color('red')
        if n_elements(d.pa) LT 2 then break
        h=histogram(d.pa, binsize=0.1, locations=loc);in each annulus, search for clump as overabundance of points near a given PA (loc is vector holding PA-value of histogram bins)
        if diagnostics eq 1 then begin
          if i eq 1 and j eq 0 then window,3,xsize=500,ysize=500 else wset,3
        endif
        if diagnostics eq 1 then plot, loc, h,xrange=[0,8],yrange=[0,500]

        hold = where(loc gt locs[j])                          ;look for overdensities in PA bins further out than locs[j]
        maxh = max(h[hold],imax)                            ;find the highest bin, call it imax
        if k eq 0 then sigma = robust_sigma(h)              ;before destroying any clump, we need to define a characteristic density, we'll use the dispersion of the distribution of PA values in the current annulus
        
        if maxh LT 7d*sigma then break ;it's a clump if its spread-out-ness is LESS than 7x the initial overall spread-out-ness -- i.e. if it has at least 7x the compactness of the overall initial compactness of the annulus
        
        if diagnostics eq 1 then oplot, intarr(10)+loc[hold[imax]],findgen(10)*1d6,color=fsc_color('red')    ;mark this bin

        dataclump = where(abs(data.pa - loc[hold[imax]]) LT 10d*!dtor AND data.er GE i*reff and data.er LT (i+1)*reff, complement = dkeep) ; isolate clump as data with [PA_bin - 10deg, PA_bin + 10deg] and [R_bin - 0.5R_e, R_bin + 0.5R_e]
       
       
        if diagnostics eq 1 then wset,0
        if diagnostics eq 1 then oplot,data[dataclump].x,data[dataclump].y,psym=3,color=fsc_color('green')
;        if 
;        data = struct_append(data,data[dkeep])
        data=data[dkeep]
        k = k+1
      endwhile
;       if n_elements(ikeep) GT 2 and diagnostics eq 1 then  oplot,d[ikeep].x,d[ikeep].y,psym=3,color=fsc_color('white') else oplot,d.x,d.y,psym=3,color=fsc_color('white')
     endfor
    endfor
    if ii eq 0 then datta=data ;this will collect cleaned data from each depth bin
    if ii ne 0 then datta=struct_append(datta,data)
endfor
    
;test that datta really holds clump-less data
undefine,data
undefine,d
stop
*(self.data) = datta


END
;-----------------------------------------------------------------------------------------------

PRO mapsim_checkl::angular_momentum_sweep, L_rad
  compile_opt idl2, hidden

d = *(self.data)
i = 0
holdAM = replicate({bincenter:0d,AM:0d},1) ;structure to hold angular momentum for each bin

WHILE 1 DO BEGIN
;binsize = 0.15 Re
;increment size = 0.1 R_e
  if i eq 0 then begin
    d_bin = d[where(d.r GT 0 and d.r LT (i/10.)*self.reff + 0.15*self.reff)] ; initial bin center at 0.1Re, r_min=0, r_max=0.25R_e
    holdAM1 = holdAM
  endif else begin 
    d_bin = d[where(d.r GT (i/10.)*self.reff - 0.25 and d.r LT (i/10.)*self.reff + 0.25)]
  endelse
  
  mass0 = d_bin.mass
  mass = transpose([[mass0],[mass0],[mass0]])         ;mass of each particle (identical mass particles)
  coord0 = transpose([[d_bin.x],[d_bin.y],[d_bin.z]])     ;create a matrix to hold the coordinates
  vCoord0 = transpose([[d_bin.vx],[d_bin.vy],[d_bin.vz]])     ;create a matrix to hold the velocity coordinates


  aM = crosspn(coord0, vcoord0*mass)            ;calculate each component (x,y,z) of the angular momentum for each particle (mv x r)
  totAM = total(aM,2,/double)               ;total angular momentum in each direction (x,y,z)
  magAM = sqrt(totAM[0]^2d + totAM[1]^2d + totAM[2]^2d) ;magnitude of the net angular momentum vector
  
  holdAM1.bincenter = (i/10.)*self.reff
  holdAM1.AM = magAM

  holdAM = struct_append(holdAM,holdAM1)

  if (i/10.)*self.reff + 0.15*self.reff GT 29d then break

  i++

ENDWHILE

p=plot(holdam.bincenter,holdam.am,linestyle=6,symbol='star',sym_size=0.5,sym_filled=1,color='blue',xrange=[0,30],xtitle='r',ytitle='L')

undefine,i
  ;now step through the radial |L| profile in a moving window of width approx. 5 kpc
  
i=0
holdmeanL = replicate({windowcenter:0d, meanL:0d},1)
WHILE 1 DO BEGIN

  if i eq 0 then begin
    meanL = mean(holdAM[where(holdAM.bincenter GT 0d and holdAM.bincenter LT 5d)].am)
    holdmeanL1 = holdmeanL
    holdmeanL1.windowcenter = 2.5d
    holdmeanL1.meanL = meanL
  endif else begin
    meanL = mean(holdAM[where(holdAM.bincenter GT i*5d and holdAM.bincenter LT (i+1)*5d)].am)
    holdmeanL1.windowcenter = median(holdAM[where(holdAM.bincenter GT i*5d and holdAM.bincenter LT (i+1)*5d)].bincenter)
    holdmeanL1.meanL = meanL
  endelse

  
  holdmeanL = struct_append(holdmeanL,holdmeanL1)

  if (i+1)*5d GT 29d then break
  
  i++
  
ENDWHILE

L_rad = holdmeanL[where(holdmeanL.meanl eq max(holdmeanL.meanl))].windowcenter


END
;-----------------------------------------------------------------------------------------------

PRO mapsim_checkL::realign_angular_momentum, magAM, rot=rot, am_rad=am_rad, kk=kk, outercalc=outercalc
  compile_opt idl2, hidden

  d = *(self.data)
  
  
    
    
  if outercalc eq 1 then r_e = am_rad/4d

  if (n_elements(rot) gt 0) then self.rot = double(rot)
  if (n_elements(am_rad) gt 0) then self.am_rad = double(am_rad)
  ;align angular momentum
  if ((min(self.rot) eq 0d) and (max(self.rot) eq 0d)) then begin
    if outercalc eq 1 then begin 
      d_inner = d[where(d.r gt (am_rad - r_e/4d) and d.r lt (am_rad + r_e/4d))] ;binwidth of 0.5 R_e.
    endif else begin
    d_inner = d[where(d.r lt self.am_rad)]                ;isolate the data within radius am_rad
    endelse
    coord0 = transpose([[d_inner.x],[d_inner.y],[d_inner.z]])     ;create a matrix to hold the coordinates
    vCoord0 = transpose([[d_inner.vx],[d_inner.vy],[d_inner.vz]])     ;create a matrix to hold the velocity coordinates
    newd = am2z_angmom(coord0, vcoord0, d_inner.mass, magAM, struc=d, rot=xrot, /verbose)   ;re-orient the data to have all the angular momentum in the z-direction
    self.rot = xrot
  endif else begin
    d_coord = transpose([[d.x],[d.y],[d.z]]) ## self.rot
    num = n_elements(d_coord[0,*]) - 1l
    d.x = (d_coord[0,*])[0:num]
    d.y = (d_coord[1,*])[0:num]
    d.z = (d_coord[2,*])[0:num]
    d_vcoord  = transpose([[d.vx],[d.vy],[d.vz]]) ## self.rot
    d.vx = (d_vcoord[0,*])[0:num]
    d.vy = (d_vcoord[1,*])[0:num]
    d.vz = (d_vcoord[2,*])[0:num]
  endelse
  *(self.data) = d

;    window,4,xsize=500,ysize=500
;    plot,d.x,d.z,psym=3,xrange=[-10,10],yrange=[-10,10],title='Principle Axis Frame'
END

;-----------------------------------------------------------------------------------------------
PRO mapsim_checkL::coord_transform
  compile_opt idl2, hidden

  d = *(self.data)
  d = arr_struct({x:d.x,y:d.z, z:d.y, vx:d.vx, vz:d.vy, vy:d.vz, mass:d.mass})  ;change of variables. xhori,yvert,zoutofpage, vz=radial velocity
  *(self.data) = d

END
;-----------------------------------------------------------------------------------------------

PRO mapsim_checkL::rotate_around_axis, reff, sige, xrot=xrot, yrot=yrot, zrot=zrot, diagnostics=diagnostics
  compile_opt idl2, hidden

data = *(self.data)

  


data = struct_addtags(data, arr_struct({pa:pacalc(data.x, data.z)*!dtor}))  ;calculate position angle
a=data[where(data.PA GT 0d and data.PA LT 90d*!dtor)]
;  if diagnostics eq 1 then begin
;    window,0,xsize=500,ysize=500
;    plot,data.x,data.z,psym=3,xrange=[-10,10],yrange=[-10,10]
;    oplot,a.x,a.z,psym=3,color=fsc_color('blue')
;  endif
  
  
  if (n_elements(xrot) gt 0) then rot = xrot else rot = 0d
  if (n_elements(yrot) gt 0) then rot = [rot, yrot] else rot = [rot, 0d]
  if (n_elements(zrot) gt 0) then rot = [rot, zrot] else rot = [rot, 0d]
  rot = rot*!dtor

  if (max(abs(rot)) gt 0d) then begin
    d = *(self.data)
    foreach ax, ['x','y','z'] do begin
      if (n_elements(lambda) gt 0) then undefine, lambda
      case ax of
        'x' : if (rot[0] ne 0d) then lambda = double([ [1, 0, 0], [0, cos(rot[0]), -sin(rot[0])], [0, sin(rot[0]), cos(rot[0])] ])  ;rotation matrix about the x-axis (y?)
        'y' : if (rot[1] ne 0d) then lambda = double([ [cos(rot[1]), 0, sin(rot[1])], [0, 1, 0], [-sin(rot[1]), 0, cos(rot[1])] ])  ;rotation matrix about the z-axis (x?)
        'z' : if (rot[2] ne 0d) then lambda = double([ [cos(rot[2]), -sin(rot[2]), 0], [sin(rot[2]), cos(rot[2]), 0], [0, 0, 1] ])  ;rotation matrix about the z-axis (x?)
      endcase

      if (n_elements(lambda) gt 0) then begin                 ;if yrot is set, then pivot around the x-axis
        ;lambda = double([ [1, 0, 0], [0, cos(yrot*!dtor), -sin(yrot*!dtor)], [0, sin(yrot*!dtor), cos(yrot*!dtor)] ])  ;rotation matrix about the x-axis (y?)
        ;lambda = double([ [cos(yrot*!dtor), -sin(yrot*!dtor), 0], [sin(yrot*!dtor), cos(yrot*!dtor), 0], [0, 0, 1] ])  ;rotation matrix about the z-axis (x?)
        d_coord = transpose([[d.x],[d.y],[d.z]]) ## lambda
        num = n_elements(d_coord[0,*]) - 1l
        d.x = (d_coord[0,*])[0:num]
        d.y = (d_coord[1,*])[0:num]
        d.z = (d_coord[2,*])[0:num]
        d_vcoord  = transpose([[d.vx],[d.vy],[d.vz]]) ## lambda
        d.vx = (d_vcoord[0,*])[0:num]
        d.vy = (d_vcoord[1,*])[0:num]
        d.vz = (d_vcoord[2,*])[0:num]
        self.rot = self.rot ## lambda
      endif
    endforeach
    *(self.data) = d
    self.axisrot = rot
  endif

  data=*(self.data)
    data = struct_addtags(data, arr_struct({pa:pacalc(data.x, data.z)*!dtor}))  ;calculate position angle
;
match,data.id,a.id,suba,subb,/SORT
b=data[suba]
;  if diagnostics eq 1 then begin
;    window,1,xsize=500,ysize=500
;    plot,data.x,data.z,psym=3,xrange=[-10,10],yrange=[-10,10]
;    oplot,b.x,b.z,psym=3,color=fsc_color('blue')
;;    stop
;  endif
;  x,y,z,vx,vy,vz coordinates attached to each particle prior to rotation. After rotation of data with coordinate axes x,y,z, the data are assigned NEW x,y,z,vx,vy,vz coords!

END

;-----------------------------------------------------------------------------------------------

PRO mapsim_checkL::cull_data, str, data
  compile_opt idl2, hidden

  if (n_elements(str) gt 0) then begin
    data = *(self.data)
    dummy = execute('ind = where('+str+',count)')
    if (count eq 0) then begin
      message, 'All of the objects were culled...', /info
      stop
    endif else *(self.data) = data[ind]
  endif
END

;-----------------------------------------------------------------------------------------------

FUNCTION mapsim_checkL::kMeans, d0, nClusters, var=var, clusterLocs=clusterLocs, numPnts=numPnts
  d = d0

  arr = float(transpose([[d.x],[d.y],[d.z]]))                               ;create an array of object positions

  clusterLocs = clust_wts(arr, n_clusters=nClusters)                            ;determine the cluster centers
  result = cluster(arr, clusterLocs, n_clusters=nClusters)                        ;determine which cluster each object "belongs" to

  ; device, decomposed=0
  ; loadct, 39, /silent
  ; maxcolor = 200
  ; colors = round(findgen(nClusters)*maxColor/(float(nClusters) - 1) + (255 - maxColor)/2d)
  ; window, 0
  ; plot, d.x, d.z, psym=3
  ; for i=0,(nClusters - 1) do oplot, d[(ind = where(result eq i))].x, d[ind].z, psym=3, color=colors[i]

  numPnts = fltarr(nClusters)                                       ;the number of points belonging to each sluster
  var = 0d
  for i=0,(nClusters - 1) do begin                                    ;for each cluster
    ind = where(result eq i, count)                                   ;retrieve the indices of objects belonging to this cluster
    var += total(( arr[*,[ind]] - cmreplicate(clusterLocs[*,i],count) )^2d)/float(n_elements(d))    ;measure the total variance of the distance residuals
    numPnts[i] = count
  endfor

  return, result
END

;-----------------------------------------------------------------------------------------------

FUNCTION mapsim_checkL::getNClusters, d0, frac
  d = d0
  if (n_elements(frac) eq 0) then frac = 0.1                                ;fraction of objects to use
  if (n_elements(itime) gt 0) then undefine, itime
  ;time, itime
  num = round(n_elements(d)*frac)                                     ;the number of objects to use in determining the optimum number of clusters
  d = d[((lindgen(n_elements(d)))[sort(randomu(seed,n_elements(d)))])[0:num-1l]]              ;select a random sample of NUM objects

  nIters = 10                                               ;number of cluster sizes to try
  nClusterArr = indgen(nIters)+1                                      ;create the cluster sizes array
  varArr = fltarr(nIters)                                         ;create an array to hold the total variance of cluster-point distance residuals
  minClustDiff = fltarr(nIters)                                     ;create an array to hold the minimum distance between any two clusters

  k = 0                                                 ;counting variable to record iterations
  foreach nClusters, nClusterArr do begin                                 ;for each specified number of clusters
    result = self->kMeans(d, nClusters, var=var, clusterLocs=clusterLocs, numPnts=numPnts)
    varArr[k] = var

    ;   window, 2
    ;   plot, nClusterArr, sqrt(varArr), psym=-6

    if (nClusters gt 1) then begin                                    ;if there is more than one cluster
      minClustDiff[k] = min(distance_measure(clusterLocs))                      ;measure the minimum distance between clusters
      ;     window, 1
      validity = varArr/minClustDiff                                  ;validity value from Siddheswar Ray and Rose H. Turi 1999, 'Determination of Number of Clusters in K-Means Clustering and Application in Colour Image Segmentation'
      ;     plot, nClusterArr, validity, psym=-6
      if ((validity[k] - validity[k-1]) gt 0) then break                        ;if this is greater than 0, then the last number of clusters is the best
    endif
    ;time, itime
    ++k
  endforeach
  return, nClusters - 1
END

;-----------------------------------------------------------------------------------------------

FUNCTION mapsim_checkL::getReff
  compile_opt idl2, hidden

  d0 = *(self.data)
  d = d0[where( (sqrt(d0.x^2d + d0.y^2d + d0.z^2d) lt 80) )]                        ;cull data inside 50 kpc

  ; nClusters = self->getNClusters(d, frac)                                 ;determine the optimum number of clusters using k-means clustering
  ;
  ; result = self->kMeans(d, nClusters, var=var, clusterLocs=clusterLocs, numPnts=numPnts)          ;run kMeans using the determined number of clusters
  ;
  ; dummy = max(numPnts,imax)                                       ;find the biggest cluster (this is the main galaxy)
  ; d = d[where(result eq imax)]                                      ;isolate the biggest cluster

;    d = arr_struct({x:d.x, y:d.z, r:sqrt(d.x^2d + d.z^2d), mass:d.mass, vel:d.vy})              ;change of variables
  d = arr_struct({r:sqrt(d.x^2d + d.z^2d), mass:d.mass, vel:d.vx})              ;change of variables

  ;  d = arr_struct({x:d.x, y:d.y, z:d.z, r:sqrt(d.x^2d + d.y^2d + d.z^2d), mass:d.mass, vel:d.vy})
  d = d[sort(d.r)]                                            ;sort by projected radius
  prof = total(d.mass, /cumulative)/total(d.mass)                             ;measure the normalized cumulative mass profile
  self.reff = interpol(d.r, prof, 0.5)                                  ;measure the radius at the half-mass mark

;    window, 2
;    plot, d.r, prof, xtitle='Projected Radius [kpc]', ytitle='Normalized Cumulative Mass', $
;      title='R_e = '+roundx(self.reff,2)+' kpc'
;      stop
  ; wset, 0
  ;  print, 'R_e = '+roundx(self.reff,2)+' kpc'
  *(self.mainGalaxy) = d                                          ;define a pointer to the object data for the largest kmeans cluster
  return, self.reff
END
;-----------------------------------------------------------------------------------------------
FUNCTION mapsim_checkL::getReff2
  compile_opt idl2, hidden

  d0 = *(self.data)
  d = d0[where( (sqrt(d0.x^2d + d0.y^2d + d0.z^2d) lt 50) )]                        ;cull data inside 50 kpc

  ;  nClusters = self->getNClusters(d, frac)                                 ;determine the optimum number of clusters using k-means clustering
  ;
  ;  result = self->kMeans(d, nClusters, var=var, clusterLocs=clusterLocs, numPnts=numPnts)          ;run kMeans using the determined number of clusters
  ;
  ;  dummy = max(numPnts,imax)                                       ;find the biggest cluster (this is the main galaxy)
  ;  d = d[where(result eq imax)]                                      ;isolate the biggest cluster

  d = arr_struct({x:d.x, z:d.z, r:sqrt(d.x^2d + d.y^2d + d.z^2d), mass:d.mass})              ;E's are to a good approximation transparent -- so we use full 3D r

  d=d[sort(abs(d.x))]
  major = d[where(d.z GT -0.01*max(d.z) and d.z LT 0.01*max(d.z) and d.x GT -max(d.x) and d.x LT max(d.x))]
  ;  major = d[where(d.z GT -0.01*max(d.z) and d.z LT 0.01*max(d.z) and d.x GT 0d and d.x LT max(d.x))]
  d=d[sort(abs(d.z))]
  minor = d[where(d.x GT -0.01*max(d.x) and d.x LT 0.01*max(d.x) and d.z GT -max(d.z) and d.z LT max(d.z))]
  ;  minor = d[where(d.x GT -0.01*max(d.x) and d.x LT 0.01*max(d.x) and d.z GT 0d and d.z LT max(d.z))]

  profa = total(major.x, /cumulative)/total(major.x)
  profb = total(minor.z, /cumulative)/total(minor.z)
  reffa = interpol(abs(major.x), profa, 0.5)
  reffb = interpol(abs(minor.z), profb, 0.5)
  self.reff = sqrt(reffa*reffb) ; Binney & Merrifield p186


  window, 2
  plot, abs(major.x), profa, xtitle='Major (solid) and Minor (dotted) Axis [kpc]',ytitle='Normalized Cumulative Number Density'
  oplot, abs(minor.z), profb, linestyle=3
  print, 'R_e = '+roundx(self.reff,2)+' kpc'
  stop
  *(self.mainGalaxy) = d                                          ;define a pointer to the object data for the largest kmeans cluster
  return, self.reff
END
;-----------------------------------------------------------------------------------------------

FUNCTION mapsim_checkL::getSige
  compile_opt idl2, hidden

  if (self.reff eq 0d) then reff = self->getReff()                            ;get the effective radius if it hasn't been determined
  d = *(self.mainGalaxy)                                          ;get the largest kmeans cluster
  prof = sqrt(total(d.mass*d.vel^2d,/cumulative)/total(d.mass,/cumulative))               ;measure the sigma_e profile
  self.sige = interpol(prof, d.r, self.reff)                                ;record the sigma_e value at the half-mass mark
;    print, 'sig_e = '+roundx(self.sige,2)+' km/s'
;   plot, d.r, prof, xtitle='Projected Radius [kpc]', xrange=[0,30], $
;     ytitle=textoidl('\sigma_e = (\Sigma M(R)*V(R) / \Sigma M(R))^{1/2}'), $
;     title=textoidl('\sigma_e = '+roundx(self.sige,2))
;   oplot, intarr(10)+self.reff, findgen(10)/9.*self.sige, linestyle=2, color=fsc_color('red')
;   xyouts, self.reff+1, 10, textoidl('R_e = '+roundx(self.reff,2)+' kpc')

  return, self.sige
END
;-----------------------------------------------------------------------------------------------
PRO mapsim_checkL::kinemetry, xrot, yrot, zrot, nrots, l, kk, f1, q, vfromfixPA, cosmo, diagnostics=diagnostics, dovel=dovel, dopa=dopa, twobins=twobins, onebin=onebin, adbin=adbin, $
  fullprofile=fullprofile, data, outx, outy, peakvQ, peakplot=peakplot
  compile_opt idl2, hidden
  
  
  

  reff = self.reff
  sige = self.sige
  d = *(self.data)



  d = arr_struct({x:d.x, y:d.z, z:d.y, vx:d.vx, vy:d.vz, vz:d.vy, r:sqrt(d.x^2d + d.z^2d), mass:d.mass})              ;change of variables (x,z,y, -> x,y,z and vx,vz,vy -> vx,vy,vz),
                                                                                                                      ;otherwise "pacalc" doesn't work properly
  


  data=d  ;x horizontal, z vertical (from am2z function)
  data = struct_addtags(data, arr_struct({pa:pacalc(data.x, data.y, 0d)*!dtor}))  ;calculate position angle
  data = struct_addtags(data, arr_struct({errvel:data.vz}))
  data = struct_addtags(data, arr_struct({vel:data.vz}))
  data.errvel = 0d

  data.r = sqrt(data.x^2 + data.y^2)
  data = data[where(finite(data.pa))]
  
;window,5,xsize=500,ysize=500
;plot,data.x,data.y,psym=3,xrange=[-10,10],yrange=[-10,10],title='Principle Axis Frame'
 diagnostics = 0
  Axisratio_estimate, data, reff, sige, q, diagnostics=diagnostics, cosmo
;  q=0.8
diagnostics = 1
;  data = struct_addtags(data, arr_struct({pa:pacalc(data.x, data.y, 0d)*!dtor}))  ;calculate position angle
;  data = struct_addtags(data, arr_struct({errvel:data.vz}))
;  data = struct_addtags(data, arr_struct({vel:data.vz}))
;  data.errvel = 0d


;--------------------------------------------------------------------------<<<<<<<<Adaptive Binning>>>>>>>>
IF adbin EQ 1 THEN BEGIN


obs=0

reff=self.reff
sige=self.sige


IF fullprofile EQ 1 THEN BEGIN
  inc = 200d < n_elements(data)
  mininc = 150d
  incr = 500d
  IF cosmo EQ 1 THEN BEGIN
   inc = 200d < n_elements(data)
   mininc = 150d
   incr = 1000d
  ENDIF
  rrange = 15*reff
  profplot = 1
ENDIF


IF twobins EQ 1 THEN BEGIN
innerbinwidth = 0.2*self.reff
innerbinQ = q
innerbinPA = 0d*!dtor
Nbootstraps_in = 1
Nbootstraps_out = 1
if cosmo eq 1 then innerbinINC = 100d
if cosmo eq 1 then innerbinINCR = 500d

if cosmo eq 0 then innerbinINC = 500;250 ;apparently this is for both annuli
if cosmo eq 0 then innerbinINCR = 500;250
profplot = 0
ENDIF

obs = 0
RollingBins_Kinemetry, data, q, bin, twobins=twobins, fullprofile=fullprofile, diagnostics=diagnostics, cosmo=cosmo, reff, sige, $
                            innerbinwidth,innerbinQ,innerbinPA,Nbootstraps_in,Nbootstraps_out,innerbinINC,innerbinINCR,profplot,F1,obs=obs, $
                            mininc,inc,incr,rrange, outx, outy, l, peakvQ, peakplot=peakplot

  ENDIF

  ;-----------------------------------------------------------------------------------------------

;  d = arr_struct({x:d.x, y:d.z, z:d.y, vx:d.vx, vy:d.vz, vz:d.vy, r:sqrt(d.y^2d + d.z^2d), mass:d.mass})              ;change of variables (swap x,y and vx,vy)
  undefine,bin
  undefine,plotdat
;  *(self.data) = data

END


;-----------------------------------------------------------------------------------------------
PRO mapsim_checkL::vmax
  compile_opt idl2, hidden


  d = *(self.data)

  reff = self.reff
  d.x = abs(d.x)
  d.vy = abs(d.vy)
  d.y = abs(d.y)
  d.x = d.x/self.reff
  d.z = d.z/self.reff

  xbounded = d[where(d.x lt 7)]
  xybounded = xbounded[where((xbounded.y gt -0.1) and (xbounded.y lt 0.1))]

  p = plot(d.x,d.z,linestyle=6,symbol='dot')
  p = plot(xybounded.x,xybounded.y,linestyle=6,symbol='plus',color='red',/overplot)
  stop
  plt = plot(xybounded.x,xybounded.vy,symbol='o',xrange=[0,7],yrange=[0,200],linestyle=6,thick=2)

  xfit = (indgen(700))/100.

  stop
END

;-----------------------------------------------------------------------------------------------


PRO mapsim_checkL::voronoi_bin, restore=restore, quick=quick, targetN=targetN, pixelSize=pixelSize, eps=eps, nthreads=nthreads
  compile_opt idl2, hidden

  if (n_elements(targetN) gt 0) then self.targetN = targetN
  if (n_elements(pixelsize) gt 0) then self.pixelsize = pixelsize
  if (n_elements(eps) gt 0) then self.eps = eps
  if (n_elements(nthreads) gt 0) then self.nthreads = nthreads

  if (n_elements(quick) gt 0) then begin
    ;load voronoi default save file
  endif

  data = *(self.data)
  ;only use the following if you have NOT already run the kinemetry routine
  ;-------------------------------------------------------------------
;  data = arr_struct({x:d.x, y:d.z, vel:d.vy})                         ;change of variables
;  d = struct_trimtags(d, except=['ID','X','Y','Z','VX','VY','VZ','R'])
;  if (type(d,/string) eq 'structure') then data = struct_addtags(data, d)
  ;-------------------------------------------------------------------
  
  *(self.obj) =  vormap_fun(data, self.targetN, eps=self.eps, pixelSize=self.pixelSize, $   ;create and accrete into voronoi bins
    nthreads=self.nthreads, restore=restore, binNumber=binNumber)

  newData = *(self.data)
  *(self.data) = struct_addtags(newData, arr_struct({binNumber:binNumber}))


END
;-----------------------------------------------------------------------------------------------

PRO mapsim_checkL::save, outfile
  compile_opt idl2, hidden

  pre = 'map_N='  + strtrim(round(self.targetN),1) + '_' + (strsplit(self.file,'/',/extract))[-1]

  all = create_struct(name=obj_class(self))
  struct_assign, self, all
  summary = struct_trimtags(all, select=['FILE','TARGETN','ROT','REFF','SIGE','AXISROT'])
  summary = struct_addtags(summary, {nBins:n_elements(*(self.obj))})
  summary = struct_addtags(summary, {total_mass:total((*(self.data)).mass)})

  k = 0
  while 1 do begin
    ;   outfile = pre + 'DATA_v' +strtrim(k,2)+'.fits'
    outfile = 'DATA_v' +strtrim(k,2)+'.fits'
    dummy = file_search(outfile,count=nfile)
    if ~nfile then begin
      mwrfits, *(self.obj), outfile
      mwrfits, summary, outfile
      mwrfits, *(self.data), outfile
      break
    endif
    ++k
  endwhile
END

;-----------------------------------------------------------------------------------------------

FUNCTION mapsim_checkL::getObj
  compile_opt idl2, hidden
  return, *(self.obj)
END

;-----------------------------------------------------------------------------------------------

PRO mapsim_checkL__define
  compile_opt idl2, hidden

  void={mapsim_checkL, $
    file:'', am_rad:0d, targetN:0., pixelsize:0., eps:0d, nthreads:0, $
    data:ptr_new(), rot:dblarr(3,3), binsize:0d, reff:0d, sige:0d, $
    axisrot:dblarr(3), obj:ptr_new(), mainGalaxy:ptr_new()}

END




