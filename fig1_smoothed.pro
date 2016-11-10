PRO Fig1_smoothed, observations, rarr, file=file, struct=struct, rawdat=rawdat, obs=obs, qplot=qplot, dblsmooth=dblsmooth


struct = 0
file = 1
qplot = 1
dblsmooth = 1
rawdat = 1


obs = 1
cosmo = 0
w1to1 = 0
d1to1 = 0
w1to3 = 0
d1to3 = 0

;secondary boxcar average smoothing with constant increment
;dblsmooth = 0
IF qplot eq 0 then dblsmooth=1
IF dblsmooth EQ 1 THEN width = 5 ELSE width = 0

;k=0
;FOR k = 0,3 DO BEGIN
;
;if k eq 0 then begin
;  w1to1 = 0
;  d1to1 = 1
;  w1to3 = 0
;  d1to3 = 0
;  cosmo = 0
;  obs = 0
;endif
;if k eq 1 then begin
;  w1to1 = 0
;  d1to1 = 0
;  w1to3 = 0
;  d1to3 = 1
;  cosmo = 0
;  obs = 0
;endif
;if k eq 2 then begin
;  w1to1 = 1
;  d1to1 = 0
;  w1to3 = 0
;  d1to3 = 0
;  cosmo = 0
;  obs = 0
;endif
;if k eq 3 then begin
;  w1to1 = 0
;  d1to1 = 0
;  w1to3 = 1
;  d1to3 = 0
;  cosmo = 0
;  obs = 0
;endif
;observations


if obs eq 1 then begin
  if file eq 1 then begin
    readcol, 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\N4564_08_11_2016_3.dat', er, v, id, /silent
;   readcol, 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\obs_n821_07_27_2016.dat', er, v, id, /silent
  endif
  
  if struct eq 1 then begin
    er = reform(observations[0,*])
    v = reform(observations[1,*])
    id = reform(observations[2,*])  
  endif
    
  ob = replicate({er:0d, v:0d, id:0d},n_elements(er))
  ob.er = er
  ob.v = v
  ob.id = id
  masterbin = ob
  title = 'NGC 4564';'Observations'
  colors = ['dodger blue', 'dark orange', 'lime green', 'royal blue', 'red', 'grey', 'orange']
;  colors = ['dark goldenrod', 'gold', 'goldenrod', 'dark orange', 'yellow', 'orange', 'khaki']
;  colors = 'gold'
endif

;cosmo sims
if cosmo eq 1 then begin
;  readcol, 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\cosmo_kinresults_07_26_2016.dat', er, v, id, /silent
  readcol, 'C:\Users\Justin\Dropbox\justin\kinemetry\2016_07_28\cosmo_profs.dat',er,v,sig,id,/silent

  cmo = replicate({er:0d, v:0d, id:0d},n_elements(er))
  cmo.er = er
  cmo.v = v
  cmo.id = id
  masterbin = cmo
  title = 'Major Merger Simulations'
  colors = ['cadet blue', 'steel blue', 'cornflower', 'dodger blue', 'royal blue', 'medium blue', 'blue', 'dark blue', 'navy', 'midnight blue', 'sky blue', 'turquoise', 'aqua', 'teal', 'dark cyan', 'aquamarine']
;  colors = 'lime green'
endif

;wet 1:1 mergers
if w1to1 eq 1 then begin
;  readcol, 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\w1to1_kinresults_07_26_2016.dat',er, v, id, /silent
  readcol, 'C:\Users\Justin\Dropbox\justin\kinemetry\2016_07_28\w1to1_profs.dat',er,v,sig,id,/silent
  w1 = replicate({er:0d, v:0d, id:0d},n_elements(er))
  w1.er = er
  w1.v = v
  w1.id = id
  masterbin = w1
  title = 'Major Merger Simulations'
  colors = ['red', 'dark red', 'maroon', 'firebrick', 'brown', 'indian red', 'crimson', 'tomato']
;  colors='crimson'
endif

;dry 1:1 mergers
if d1to1 eq 1 then begin
  readcol, 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\d1to1_kinresults_07_26_2016.dat',er, v, id, /silent
  d1 = replicate({er:0d, v:0d, id:0d},n_elements(er))
  d1.er = er
  d1.v = v
  d1.id = id
  masterbin = d1
  title = 'Hoffman Dry 1:1 Mergers'
  colors = ['deep pink', 'fuchsia', 'magenta', 'medium purple', 'medium orchid', 'orchid', 'violet', 'plum', 'thistle', 'pink', 'light pink', 'pale violet red', 'medium violet red', $
    'hot pink', 'indigo', 'medium slate blue', 'slate blue', 'purple', 'dark magenta', 'blue violet', 'dark violet', 'dark orchid', 'cornflower', 'royal blue', $
    'deep pink', 'fuchsia', 'magenta', 'medium purple', 'indigo']
;colors = 'plum'
endif

;wet 1:3 mergers
if w1to3 eq 1 then begin
  readcol, 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\w1to3_kinresults_07_26_2016.dat',er, v, id, /silent
;  readcol, 'C:\Users\Justin\Dropbox\justin\kinemetry\2016_07_28\w1to3_profs.dat',er,v,sig,id,/silent
  w3 = replicate({er:0d, v:0d, id:0d},n_elements(er))
  w3.er = er
  w3.v = v
  w3.id = id
  masterbin = w3
  title = 'Hoffman Wet 1:3 Mergers'
  colors = ['forest green', 'dark green', 'olive', 'olive drab', 'dark olive green', 'sea green', 'dark sea green', 'pale green', 'light green', 'spring green', 'medium spring green', $
    'medium sea green', 'lime', 'lime green', 'green']
;colors = 'medium blue'
endif

;dry 1:3 mergers
if d1to3 eq 1 then begin
  readcol, 'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\d1to3_kinresults_07_26_2016.dat',er, v, id, /silent
;  readcol, 'C:\Users\Justin\Dropbox\justin\kinemetry\2016_07_28\d1to3_profs.dat',er,v,sig,id,/silent

  d3 = replicate({er:0d, v:0d, id:0d},n_elements(er))
  d3.er = er
  d3.v = v
  d3.id = id
  masterbin = d3
  title = 'Hoffman Dry 1:3 Mergers'
  colors = ['forest green', 'dark green', 'olive', 'olive drab', 'dark olive green', 'sea green', 'dark sea green', 'pale green', 'light green', 'spring green', 'medium spring green', $
    'medium sea green', 'lime', 'lime green', 'green']
;colors='silver'
endif


;==================================MOVING WINDOW VARIABLES
if d1to1 eq 1 or w1to1 eq 1 or w1to3 eq 1 or d1to3 eq 1 then begin
  xsize = 8 ;xrange of profile plot
  wsize = 0.2;moving window width
  maxdx = 0.05d;0.2d;moving window increment for outermost data
  len = 0.005d;0.005d;smoothing length scale
endif
  
if obs eq 1 then begin
  xsize = 8 ;xrange of profile plot
  wsize = 1d;moving window width
  maxdx = 0.2d;0.09d;moving window increment for outermost data
  len = 0.005d;0.001d;smoothing length scale
endif  
  
if cosmo eq 1 then begin
  xsize = 8
  wsize = 0.3
  maxdx = 0.1d
  len = 0.005d
endif

plotdat = replicate({er:0d, v:0d, sig:0d},1)
tmp = plotdat
;==========================================

unq = uniq(masterbin.id)
ngals = n_elements(uniq(masterbin.id))

if ngals eq 1 then begin
  nrows = 1
  ncols = 1
endif else begin
  if ngals mod 2 eq 0 then nrows = ngals/2d
  if ngals mod 2 ne 0 then nrows = (ngals+1)/2d
  if ngals eq 7 then nrows = 4
  ncols = 2
endelse 


FOR j = 0,ngals-1 DO BEGIN;
  bin = masterbin[where(masterbin.id eq j)]
;  bin = masterbin ;all data at once -- for figure 1 panels D&F shady regions
  window,0,xsize=500,ysize=500
  plot,bin.er,bin.v,psym=4,xrange=[0,6],yrange=[-0.3,2]



  i = 0
  WHILE 1 DO BEGIN

   if i eq 0 then begin
    currentdat = bin[where(bin.er GT 0d and bin.er LT 0.05d)]
    mn = 0d
    mx = max(currentdat.er)
    dx = 0.01
   endif

   if i gt 0 then dx = maxdx - len*alog(n_elements(currentdat))
;   if i gt 0 then dx = 0.05
  
  
  if j eq 0 then wsize = 1d else wsize = 0.5d
   if max(currentdat.er) lt wsize then begin ;regime where you still have not incremented over to window of width wsize
    currentdat = bin[where(bin.er LT mx+dx)]
    mn = min(currentdat.er)
    mx = max(currentdat.er)
   endif else begin
    currentdat = bin[where(bin.er GT mn+dx AND bin.er LT (mn+dx)+wsize)]
    mn = min(currentdat.er)
   endelse

   if n_elements(currentdat.v) LT 2 then break


   tmp.er = median(currentdat.er)
   tmp.v = median(currentdat.v)
   tmp.sig = stddev(currentdat.v)


   if i eq 0 then begin
    plotdat = tmp
   endif else begin
    plotdat = struct_append(plotdat,tmp)
   endelse
      oplot,bin.er,bin.v,psym=4,color=fsc_color('white')
      if n_elements(currentdat.er) gt 2 then oplot,currentdat.er,currentdat.v,psym=4,color=fsc_color('green')
      discrete = 1
      stellar = 0
      
   if n_elements(currentdat.er) lt 2 then begin
    print,strcompress('break at'+string(tmp.er)+'R_e')
    break    
   endif
   
   if stellar eq 1 and plotdat[i].er gt 3d then break
   if discrete eq 1 and plotdat[i].er gt 8d then break
   
   print,strcompress('N='+string(n_elements(currentdat.er))+'  dx='+string(maxdx - len*alog(n_elements(currentdat))))

   i=i+1
 ENDWHILE

 hold = replicate({er:0d, v:0d, sig:0d, tag:0d},n_elements(plotdat))
 hold.er = plotdat.er
 hold.v = plotdat.v
 hold.sig = plotdat.sig
 hold.tag = j

 ;secondary smoothing
 IF dblsmooth EQ 1 THEN BEGIN
  s = smooth(hold.v,width)
;  s = smooth(s,width)
;  s = smooth(s,width)
  hold.v = s
 ENDIF


 ;-------------------------------------------------------------------find peak of rotation 
 IF qplot EQ 1 THEN BEGIN
  if j eq 0 then tmp2 = fltarr(1)
  maxr = where(hold.v eq max(hold[where(hold.er LT 5.5d)].v))
  print,strcompress('peak rotation at'+string(plotdat[maxr].er)+'R_e')
  maxr = plotdat[maxr].er
  tmp2 = maxr
  if j eq 0 then begin
    rarr = tmp2
  endif else begin
    rarr = [rarr, tmp2]
  endelse
 ENDIF
 ;-------------------------------------------------------------------




 lx = indgen(7)
 ly = indgen(7)*0
 IF rawdat EQ 1 THEN BEGIN
  p1=plot(bin.er,bin.v,linestyle=6,symbol='dot',xrange=[0,6],yrange=[-0.3,1.4],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$',layout=[ncols,nrows,j+1],/CURRENT)
  p1=plot(hold.er,hold.v,thick=2,color=colors[j],/overplot)  
  IF qplot EQ 1 THEN BEGIN
    vx = indgen(60)*0 + maxr
    vy = indgen(60)-30d
    p1=plot(vx,vy,linestyle=1,color='red',thick=2,/overplot)
    tx=text(2,1,strcompress('Peak Rot. Radius ='+string(maxr, format='(F0.4)')+'$R_e$'),font_size=8,target=p1,/DATA)
  ENDIF  
 ENDIF

 IF qplot EQ 1 and rawdat EQ 0 THEN BEGIN
  p1=plot(hold.er,hold.v,linestyle=6,symbol='dot',xrange=[0,6],yrange=[-0.3,1.4],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$',layout=[ncols,nrows,j+1],/CURRENT)
  vx = indgen(60)*0 + maxr
  vy = indgen(60) - 30d
  p1=plot(vx,vy,linestyle=1,color='red',thick=2,/overplot)
  tx=text(2,1,strcompress('Peak Rot. Radius ='+string(maxr, format='(F0.4)')+'$R_e$'),font_size=8,target=p1,/DATA)
 ENDIF

  
 IF j EQ 0 THEN BEGIN
  profs = hold
 ENDIF ELSE BEGIN
  profs = struct_append(profs, hold)
 ENDELSE


  


ENDFOR  ;close loop over N galaxies

IF qplot EQ 0 THEN BEGIN
  FOR j = 0,ngals-1 DO BEGIN
    if j eq 0 then begin
      p=plot(profs[where(profs.tag eq 0)].er, profs[where(profs.tag eq 0)].v, thick=2, color=colors[j], xrange=[0,6], yrange=[-0.3,1.4], xtitle='$R/R_e$', ytitle='$v_{rot}/\sigma_e$',font_size=18, title=title)
      p=plot(lx, ly, linestyle=1, color='grey',/overplot)
    endif else begin
      p=plot(profs[where(profs.tag eq j)].er, profs[where(profs.tag eq j)].v, thick=2, color=colors[j], xrange=[0,6], yrange=[-0.3,1.4], /overplot)
    endelse
  ENDFOR
ENDIF
    
    
    
    
;  p=plot(plotdat.er,plotdat.v,thick=2,color='grey',xrange=[0,6],yrange=[-0.3,1.4], xtitle='$R/R_e$', ytitle='$v_{rot}/\sigma_e$',font_size=18, title=title, /NODATA)
;  p=plot(plotdat.er,plotdat.v+plotdat.sig,thick=1,color='grey',/overplot)
;  p=plot(plotdat.er,plotdat.v-plotdat.sig,thick=1,color='grey',/overplot)


;134


END