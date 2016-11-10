;It's possible that the small values I'm getting for Delta (fig2B) are from too low of an estimate of sigma at 4R_e!
pro mapsim_checkL, file, nthreads=nthreads, rot=rot, xrot=xrot, yrot=yrot, zrot=zrot, restore=restore, am_rad=am_rad, targetN=targetN
  compile_opt idl2, hidden

  ;restore='vormap_29.sav'
  ngals = 16 ;num sims
  nrots = 1 ;number of random projections
  twobins = 0 ;two elliptical bins mode
  diagnostics = 1
  dovel = 1 ;turn this on to calculate delta(V_rot) from fixed PA
  dopa = 0
  fig2a = 0
  fig2b = 0
  onebin = 0
  adbin = 1
  adbinprof = 0
  adbin2bin = 1
  vfromfixPA = 1
  chklplot = 0
  outercalc = 0 ;calculate angular momentum in annulus around 4R_e. This is best for dry 1to1 mergers.
  Lsweep = 0
  voronoi = 0
  cosmo = 1
  fullprofile = 1
  peakplot = 1
  
  
  if chklplot GT 0 then begin
    a=indgen(17)/10.
    b=indgen(17)*0
      p=plot(a,b,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='medium slate blue',xtitle='Net Inner Angular Momentum', ytitle='$\Gamma_{2bin}$', /nodata)
  endif
  
  
  
  ;plotting variables (fig.2a)
  if fig2a GT 0 then begin
    a=indgen(17)/10.;x-coordinates of 1-to-1 relation line and hori. line
    b=indgen(17)*0. ;hori. line y-coords
    c=a ;1-to-1 line
    p=plot(a,b,linestyle=6,xrange=[0,1.6],$
      yrange=[-0.3,1.6],xtitle='$V/\sigma$ ($1R_e$)',$
      sym_filled=1, ytitle='$V/\sigma$ ($4R_e$)', /nodata)
    p=plot(a,b,linestyle=2,color='black',/overplot)
    p=plot(a,c,linestyle=2,color='black',/overplot)
  endif


  ;plotting variables (fig.2b)
 
  if fig2b GT 0 then begin
    ;    if twobins eq 1 then begin
    vlx = indgen(200)*0
    vly = indgen(200)
    p=plot(vlx,vly,linestyle=6,symbol='star',sym_filled=1,sym_size=0.5,$
      color='dodger blue',xrange=[-1,1],yrange=[0,180],xtitle='Normalized change in rotation $\Deltav_{4-1} [km/s]$',$
      ytitle='$\Gamma_{kin} [deg]$',/nodata)
    p=plot(vlx,vly,linestyle=1,/overplot)
  endif

;  openw,1,'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\wet1to3_07_05_2016.dat'
  ;endif


  for l=0,ngals-1 do begin  ;iterate the randomLOS kinemetry routine over a set of sims
    undefine,a
    undefine,am_rad
    undefine,b
    undefine,c
    undefine,mapsim_checkL
    undefine,reff
    undefine,sige
    undefine,xrot
    undefine,rot
    undefine,seed


    print,strcompress('galaxy:'+string(l+1))

    files=FILE_SEARCH('I:\Cosmo_data\LatestZ\*.dat0')
;    files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\Dissipationless\1to1\7kpc\noBulge\snap_orb_*.dat')
;        files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\all3to1\snap_orb_*.dat')
;files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\withGas\1to1\snap_orb_*.dat')
;    files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\Fig1c\*.dat0')
;    files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\all3to1\snap_orb_*.dat')
;     files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\withGas\3to1\*.dat')
;     files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\Dissipationless\3to1\7kpc\noBulge\*.dat')
;    files=FILE_SEARCH('F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\dry1to1\*.dat')
;files=FILE_SEARCH('H:\Cosmo_data\LatestZ\*.dat0')
    file = files[l]
;    file='H:\Cosmo_data\LatestZ\MW1_SZ602.a0.420.dat0'
;file = 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\HoffSims\Dissipationless\3to1\7kpc\noBulge\snap_orb_p32.dat'
    print,file

    mapsim_checkL = obj_new('mapsim_checkL', file)                ;instantiate the mapsim_checkL object

    if strmid(file,strlen(file)-4,4) eq 'dat0' then begin
      mapsim_checkL->cull_satellites, binsize=0.25           ;cull the satellites
    endif


    mapsim_checkL->cull_data, '(data.r lt 30)'               ;cull any data, set str to a string executed in a where statement



    if (n_elements(amrad) eq 0) then am_rad = mapsim_checkL->getreff() $ ;calculate the effective radius
    else am_rad = amrad

    if outercalc eq 1 then     am_rad = am_rad*5d ;calculate L from outer bin, instead of noisy inner bin

    if Lsweep eq 1 then begin
      mapsim_checkL->angular_momentum_sweep, L_rad
      am_rad = L_rad
      print,strcompress('r_Lmax ='+string(L_rad)+'[kpc]')
    endif
    

    mapsim_checkL->realign_angular_momentum, magAM, rot=rot, am_rad=am_rad, kk=kk, outercalc=outercalc     ;realign the angular momentum with the z-axis -- seems to make y horizontal and x out of page.

;  print,'Net angular momentum ='+string(abs(meantotAM))

    sige = mapsim_checkL->getsige()                    ;calculate the sigma_e value
    reff = mapsim_checkL->getreff()

    print, strcompress('pre-rotation R_e ='+string(reff)+'kpc')
    print, strcompress('pre-rotation Sigma_e ='+string(sige)+'km/s')

    a2=systime() ; for recording total runtime

;      cbinsize=reff
;      if n_elements(cosmo) eq 1 then mapsim_checkL->clump_destroyer, data, reff, sige, q, cbinsize=cbinsize, diagnostics=diagnostics

    for kk = 0,nrots-1 do begin

    ;Note on rotations and coordinates. "rotate_around_axis" uses data from *(self.data), where x->horizontal, z->vertical, y->depth
      a1=systime() ; time elapsed to analyze one viewing angle

      print,strcompress('rotation'+string(kk+1))
      xrot=0d ; Horizontal (inclinations)
      yrot=0d ; LOS (PA rotations)
      zrot=0d ; Vertical (azimuthal rotations)
      if kk ne 0 then zrot = randomu(seed,1)*360d ;changed to 360, was 0 to 180.w
;      if kk ne 0 and cosmo eq 1 then yrot = randomu(seed,1)*360d ;cosmo sims: z-axis is into the screen, so azimuthal rotations are about Y
      rotz = zrot
;      if cosmo eq 1 then rotz = yrot

      print,strcompress('z-rotation (azimuthal) = '+string(round(zrot))+'deg')
      mapsim_checkL->rotate_around_axis, reff, sige, xrot=xrot, yrot=yrot, zrot=zrot, diagnostics=diagnostics    ;set x/y/zrot to rotate around that axis
  xrot=0
  yrot=0
  zrot=0

      if kk ne 0 then begin
        while 1 do begin
          val = randomn(seed,1)
          if abs(val) LT 0.34d then break
        endwhile

        xrot = asin(val)/!dtor  ;changed from randomU to randomN
        
        print,strcompress('Ang.Momentum ='+string(val))
        endif
      rotx=xrot
      print,strcompress('x-rotation (inclination) = '+string(round(xrot))+'deg')
      mapsim_checkL->rotate_around_axis, reff, sige, xrot=xrot, yrot=yrot, zrot=zrot, diagnostics=diagnostics    ;set x/y/zrot to rotate around that axis



      reff = mapsim_checkL->getreff()                    ;calculate the effective radius ;***Warning*** V. axis = d.z, H. axis = d.x!! This is reassigned to d.x,d.y in mapsim_checkL::kinemetry module only!! ****
      sige = mapsim_checkL->getsige()                    ;calculate the sigma_e value
      print, strcompress('post-rotation R_e ='+string(reff)+'kpc')
      print, strcompress('post-rotation Sigma_e ='+string(sige)+'km/s')
      mapsim_checkL->kinemetry, xrot, yrot, zrot, nrots, l, kk, f1, q, vfromfixPA, cosmo, diagnostics=diagnostics, dovel=dovel, dopa=dopa, twobins=twobins, onebin=onebin, adbin=adbin, $
        fullprofile=fullprofile, data, outx, outy, peakvQ, peakplot=peakplot;, doaxisratio=doaxisratio



;  if diagnostics GT 0 then window,0,xsize=500,ysize=500
;  if diagnostics GT 0 then plot,data.x,data.y,xrange=[-7,7]*reff,yrange=[-7,7]*reff,psym=3

    ;      Voronoi Map
      if voronoi GT 0 then begin
          restore='C:\Users\Justin\Dropbox\justin\IDL Programs\vormap_29.sav'
          targetN=500d
          mapsim_checkL->voronoi_bin, nthreads=nthreads, restore=restore, targetN=targetN      ;create and accrete into voronoi bins
          mapsim_checkL->save, outfile                          ;save the results
          simfile='C:\Users\Justin\DATA_v17.fits';outfile
          xmap_justin, simfile
      endif


     if fullprofile eq 0 then begin
      if twobins eq 1 or adbin eq 1 then begin
        if kk eq 0 then begin
          f=f1
        endif else begin
          f = struct_append(f,f1)
          print,'delta:',f1.delta
          print,'2-bin gamma:',f1.gamma2bin
        endelse
      endif
     endif
      ;    if l eq 0 then color = 'red'
      ;    if l eq 1 then color = 'orange'
      ;    if l eq 2 then color = 'yellow'
      ;    if l eq 3 then color = 'green'
      ;    if l eq 4 then color = 'blue'
      ;    if l eq 5 then color = 'medium slate blue'

      if fig2a GT 0 then begin
        if kk ne 0 then begin
          p=errorplot(f.innerv,f.outerv,f.div,f.dov,linestyle=6,$
            symbol='o',color='medium slate blue',xrange=[0,1.6],yrange=[-0.3,1.6],$
            title='1:1 merger sims (gas-poor)',xtitle='$v_{rot}/\sigma_e$ ($1R_e$)',$
            sym_size=0.5,sym_filled=1, ytitle='$v_{rot}/\sigma_e$ ($4R_e$)', /overplot)
        endif
      endif



;      if fig2b GT 0 then begin
;        if kk ne 0 then begin
;          p=errorplot(f.delta,f.gamma2bin,f.derror,f.g2binerror,linestyle=6,symbol='o',sym_filled=1,sym_size=0.5,$
;            color='firebrick',xrange=[-1,1],yrange=[0,180],xtitle='Normalized change in rotation $\Deltav_{4-1}$ [km/s]',$
;            ytitle='$\Gamma_{kin} [deg]$',/overplot)
;        endif
;      endif
      
;      out=dblarr(8,1)
;      out[0,*]=f1.innerv
;      out[1,*]=f1.outerv
;      out[2,*]=f1.div
;      out[3,*]=f1.dov
;      out[4,*]=f1.delta
;      out[5,*]=f1.derror
;      out[6,*]=f1.gamma2bin
;      out[7,*]=f1.g2binerror



      if fig2b GT 0 then begin
        if kk ne 0 then begin
          p=plot(f.delta,f.gamma2bin,linestyle=6,symbol='o',sym_filled=1,sym_size=0.5,$
            color='dodger blue',xrange=[-1,1],yrange=[0,180],xtitle='Normalized change in rotation $\Deltav_{4-1}$ [km/s]',$
            ytitle='$\Gamma_{kin} [deg]$',/overplot)
        endif
      endif


  if fullprofile eq 0 then begin
    out=dblarr(4,1)
      out[0,*]=f1.innerv
      out[1,*]=f1.outerv
      out[2,*]=f1.delta
      out[3,*]=f1.gamma2bin


;      openu,2+kk+l,'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\wet1to3_07_05_2016.dat',/append
;      printf,2+kk+l,FORMAT='(4F)',out
;      close,2+kk+l
  endif


  
    if chklplot GT 0 then begin
      if kk GT 0 then begin
        p=plot(abs(meantotAM),f.gamma2bin,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,color='medium slate blue',/overplot)
      endif
    endif
      
      
      


      b1=systime()

      ;------calc time elapsed
      c1=strmid(a1,11,8)
      c1=strjoin(strsplit(c1,':',/extract),' ')
      h1 = strmid(c1,0,2)    ;hours
      m1 = strmid(c1,3,2)    ;minutes
      s1 = strmid(c1,6,2)    ;seconds
      c2=strmid(b1,11,8)
      c2=strjoin(strsplit(c2,':',/extract),' ')
      h2 = strmid(c2,0,2)    ;hours
      m2 = strmid(c2,3,2)    ;minutes
      s2 = strmid(c2,6,2)    ;seconds

      dh = float(h2) - float(h1)
      if dh LT 0 then dh = 24 + dh
      dm = float(m2) - float(m1)
      if dm LT 0 then dm = 60 + dm
      ds = float(s2) - float(s1)
      if ds LT 0 then ds = 60 + ds

      print,strcompress(string(dh)+' hour(s)'+string(dm)+' minute(s)'+string(ds)+' second(s)'+' elapsed for rotation'+string(kk+1))
      ;------





















      ;;undo rotations
      IF kk GT 0 THEN BEGIN

        xrot=-rotx
        yrot=0
        zrot=0
        print,strcompress('undoing (inclination) x-rotation:'+string(xrot))
        mapsim_checkL->rotate_around_axis, reff, sige, xrot=xrot, yrot=yrot, zrot=zrot, diagnostics=diagnostics    ;set x/y/zrot to rotate around that axis

        xrot=0
        yrot=0
        zrot=-rotz
        print,strcompress('undoing (azimuthal) z-rotation:'+string(zrot))
        mapsim_checkL->rotate_around_axis, reff, sige, xrot=xrot, yrot=yrot, zrot=zrot, diagnostics=diagnostics    ;set x/y/zrot to rotate around that axis
        
      ENDIF
      undefine,xrot
      undefine,yrot
      undefine,zrot




    endfor ;end loop over multiple viewing angles


if fullprofile eq 0 then begin
    if l eq 0 then begin
      ff=f
    endif else begin
      ff = struct_append(ff,f)
    endelse
endif



  if l eq 0 then begin
;    profs = fltarr(ngals+1,29)
;    profs[0,*] = outx[0:28]
;    pg1 = gaussfit(outx[0:20],outy[0:20])
;    pg2 = gaussfit(outx[20:40],outy[20:40])
;    pg3 = gaussfit(outx[40:60],outy[40:60])
;    pg = fltarr(61)
;    pg[0:20] = pg1
;    pg[20:40] = pg2
;    pg[40:60] = pg3
;    p=plot(outx,-outy,linestyle=6,symbol='o',sym_size=0.5,sym_filled=1,xrange=[0,6],yrange=[-0.3,1.4],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$',title='Dry 1:1 orbit j4')
;    p=plot(outx,-profs[1,*],thick=2,color='medium slate blue',/overplot)
    profs = replicate({x:0d, y:0d, tag:0l},n_elements(outx))
    profs.x = outx
    profs.y = outy
    profs.tag = l
  endif else begin
    tmp = replicate({x:0d, y:0d, tag:0l},n_elements(outx))
    tmp.x = outx
    tmp.y = outy
    tmp.tag = l
    profs = struct_append(profs,tmp)
  endelse



  endfor  ;end loop over multiple sims
  
  if fullprofile eq 1 then begin
     linex=indgen(7)
     liney=indgen(7)*0
;     color = ['indian red', 'medium slate blue', 'violet', 'green', 'lime green']
     color = ['dark blue','navy','midnight blue','light sky blue','deep sky blue','sky blue','dodger blue','cornflower','cornflower','steel blue','royal blue','blue','medium blue','teal','aqua','cyan','cadet blue']
;     names = ['Wet 1:1 i', 'Dry 1:1 j', 'Dry 1:1 k', 'Wet 1:3 l23', 'Dry 1:3 l23']
     names = ['MW1 z=1.38','MW3 z=1.38','MW4 z=1.38','MW7 z=0.52','MW8 z=0.35','MW9 z=1','SFG5 z=1.27']
     for jj = 0,ngals-1 do begin
      if jj eq 0 then begin 
        p=plot(profs[where(profs.tag eq 0)].x,-profs[where(profs.tag eq 0)].y,thick=2,color=color[jj],xrange=[0,6],yrange=[-0.3,1.4],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$')
;        p=plot(profs[0,*],-profs[1,*],thick=2,color=color[0],xrange=[0,6.5],yrange=[-0.3,1.4],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$')
        p=plot(linex,liney,linestyle=1,color='gray',/OVERPLOT)
      endif
       if jj ne 0 then p=plot(profs[where(profs.tag eq jj)].x,-profs[where(profs.tag eq jj)].y,thick=2,color=color[jj],/overplot)
;      if jj ne 0 then p=plot(profs[0,*],-profs[jj+1,*],thick=2,color=color[jj+1],/OVERPLOT)
;      t=text(4.5,1.4 - jj/8.,names[jj],color=color[jj],target=p,/DATA)
     endfor
  endif
  stop
  
  
  
  
  b2=systime()


  c1=strmid(a2,11,8)
  c1=strjoin(strsplit(c1,':',/extract),' ')
  h1 = strmid(c1,0,2)    ;hours
  m1 = strmid(c1,3,2)    ;minutes
  s1 = strmid(c1,6,2)    ;seconds
  c2=strmid(b2,11,8)
  c2=strjoin(strsplit(c2,':',/extract),' ')
  h2 = strmid(c2,0,2)    ;hours
  m2 = strmid(c2,3,2)    ;minutes
  s2 = strmid(c2,6,2)    ;seconds

  dh = float(h2) - float(h1)
  if dh LT 0 then dh = 24 + dh
  dm = float(m2) - float(m1)
  if dm LT 0 then dm = 60 + dm
  ds = float(s2) - float(s1)
  if ds LT 0 then ds = 60 + ds

  print,strcompress(string(dh)+' hour(s)'+string(dm)+' minute(s)'+string(ds)+' second(s)'+' elapsed for simulation set')

  ;undo rotations








  ;---------------------------------------------------------------------make plot
;  if n_elements(fig2a) GT 0 then begin
;    a=indgen(17)/10.;x-coordinates of 1-to-1 relation line and hori. line
;    b=indgen(17)*0. ;hori. line y-coords
;    c=a ;1-to-1 line

;    p=plot(f.inner,f.outer,linestyle=6,$
;      symbol='star',color='medium slate blue',$
;      xrange=[0,1.6],yrange=[-0.3,1.6],$
;      xtitle='$V/\sigma$ ($1R_e$)',$
;      sym_filled=1, ytitle='$V/\sigma$ ($4R_e$)', /overplot)

;    p=plot(a,b,linestyle=2,color='black',/overplot)
;    p=plot(a,c,linestyle=2,color='black',/overplot)
;  endif
  ;---------------------------------------------------------------------
  ;  mapsim_checkL->extract_vmax

  ;  mapsim_checkL->vmax



;  mapsim_checkL->save                          ;save the results



  ; obj_destroy, mapsim_checkL
end


