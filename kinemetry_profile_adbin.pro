PRO Kinemetry_profile_adbin, reff, sige, q, data, rrange, bin44, bin55, bin1, bin2, twobin=twobin, binfullprof, fullprofile=fullprofile, name, stellar, discrete, inplotdat, inmaxr, outx, $
 outy, obs, cosmo, l, peakvQ, peakplot=peakplot, maxr, bin


  IF fullprofile EQ 1 THEN BEGIN
  bin=binfullprof

if n_elements(bin2) gt 1 then begin
  bin=bin1[where(bin1.er/reff LT 2d)]
  bin2=bin2[where(bin2.er/reff GT 2d)]
  bin=struct_append(bin,bin2)
;  p=plot(bin.er/reff,bin.quan/sige,linestyle=6,symbol='dot',xrange=[0,6],yrange=[-0.3,2.5])
  histdata = histogram(bin.er)
  num = n_elements(histdata)
  depdat = findgen(num)*((max(bin.er) - min(bin.er))/(num-1)) + min(bin.er)
  plot,bin.er/reff,bin.quan/sige,psym=5,yrange=[-1,2.5]
  oplot,depdat/reff,histdata/250d,color=fsc_color('green')
endif

  ;=============MOVING WINDOW VARIABLES
  if obs eq 1 then begin
    xsize = 8d ;xrange of profile plot
    wsize = 1d;moving window width
    maxdx = 0.2d;0.09d;moving window increment for outermost data
    len = 0.005d;0.001d;smoothing length scale
  endif
  
  if obs eq 0 and cosmo eq 0 then begin
    xsize = 8d
    wsize = 2d
    maxdx = 0.2d
    len = 0.005d
  endif
  
  if cosmo eq 1 then begin
    xsize = 8
    wsize = 0.8
    maxdx = 0.2d
    len = 0.005d
  endif
  
    plotdat = replicate({er:0d, v:0d, sig:0d},1)
    tmp = plotdat
    i=0
  ;=============MOVING WINDOW ALGORITHM
window,0,xsize=500,ysize=500
plot,bin.er/reff,bin.quan/sige,psym=4,xrange=[0,6],yrange=[-0.3,1.4]
WHILE 1 DO BEGIN


   if i eq 0 then begin
    currentdat = bin[where(bin.er/reff GT 0d and bin.er/reff LT 0.1d)]
    mn = 0d
    mx = max(currentdat.er/reff)
   endif
   
   
   dx = maxdx - len*alog(n_elements(currentdat))
;    dx = 0.2d
  
  if i gt 0 then begin
   if max(currentdat.er/reff) lt wsize then begin ;regime where you still have not incremented over to window of width wsize
    currentdat = bin[where(bin.er/reff LT mx+dx)]
    mn = min(currentdat.er/reff)
    mx = max(currentdat.er/reff)
   endif else begin
;    currentdat = bin[where(bin.er/reff GT mn+dx - wsize/2d AND bin.er/reff LT (tmp.er + dx) + wsize/2d)]
    currentdat = bin[where(bin.er/reff GT mn+dx AND bin.er/reff LT (mn+dx)+wsize)]
    mn = min(currentdat.er/reff)
   endelse
  endif

   tmp.er = mean(currentdat.er/reff)
   tmp.v = mean(currentdat.quan/sige)
   tmp.sig = stddev(currentdat.quan/sige)
   

   if i eq 0 then begin
    plotdat = tmp
   endif else begin
    plotdat = struct_append(plotdat,tmp)
   endelse
      oplot,bin.er/reff,bin.quan/sige,psym=4,color=fsc_color('white')
      if n_elements(currentdat.er) gt 2 then oplot,currentdat.er/reff,currentdat.quan/sige,psym=4,color=fsc_color('green')
      
;      discrete = 1
;      stellar = 0
   if n_elements(currentdat.er) lt 2 then  begin
    print,strcompress('break at'+string(tmp.er)+'R_e')
    break    
   endif
   if stellar eq 1 and plotdat[i].er gt 3d then break
   if discrete eq 1 and plotdat[i].er gt 8d then break
   print,strcompress('N='+string(n_elements(currentdat.er))+'  dx='+string(maxdx - len*alog(n_elements(currentdat))))
;if i gt 30 then stop
   i=i+1
ENDWHILE
  ;=============
  

  IF stellar EQ 1 OR discrete EQ 1 THEN BEGIN
  plotdat = plotdat[where(plotdat.er*sige LT 1000d)]
;  plotdat = plotdat[where(abs(plotdat[4,*]) NE 180d)]
  ENDIF

;  if discrete eq 1 then plotdat = plotdat[*,where(plotdat[0,*] GT 1)]
  
  if stellar eq 1 then begin 
    maxr = i-1
    inmaxr = maxr ;save for overplot
    inplotdat = plotdat ;save for overplot
  endif
  
  if discrete eq 1 then maxr = (n_elements(plotdat.er)-2)
  ;=============PLOT
   linex=indgen(9)
   liney=indgen(9)*0
   if stellar eq 1 then begin
    title = strcompress(string(name)+' stellar')
    colorv = 'medium slate blue'
    colorpa = 'orange'
   endif
   if discrete eq 1 and obs eq 1 then begin
    title = strcompress(string(name)+' discrete')
    colorv = 'medium slate blue'
    colorpa = 'green'
   endif
   if discrete eq 0 and stellar eq 0 then begin
    title = 'some simulated galaxy rv profile'
    colorv = 'lime green'
    colorpa = 'indigo'
    maxr = n_elements(plotdat[0].er)-1
   endif
   
   
   
   
   
   IF discrete eq 0 and stellar eq 0 then begin 
;    p=plot(plotdat.er,-plotdat.v,color='medium slate blue',thick=2,xrange=[0,6],yrange=[-0.3,1.4],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$')
   ENDIF ELSE BEGIN
;   p=plot(plotdat.er[0:maxr],plotdat.v[0:maxr],linestyle=0,thick=2,color=colorv,xrange=[0,8],yrange=[-1,2],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$', title=title)
;   p=plot(linex,liney,linestyle=1,color='gray',/OVERPLOT)
   ENDELSE
;   p=plot(plotdat[0,0:maxr],plotdat[1,0:maxr],linestyle=0,thick=2,color=colorv,xrange=[0,8],yrange=[-0.5,2],xtitle='R/R_e',ytitle='$v_{rot}/\sigma_e$', title=title,layout=[2,1,1])
;   p=plot(linex,liney,linestyle=1,color='gray',/OVERPLOT)
;   IF twobin EQ 0 THEN BEGIN
;    p=plot(plotdat[0,0:maxr],plotdat[4,0:maxr],linestyle=0,thick=2,color=colorpa,xrange=[0,8],xtitle='R/R_e',ytitle='kPA',yrange=[-180,180],layout=[2,1,2],/current)
;    p=plot(linex,liney,linestyle=1,color='grey',/OVERPLOT)
;   ENDIF
  outx = plotdat.er
  outy = plotdat.v
;  p=plot(bin.er/reff,-bin.quan/sige,linestyle=6,symbol='dot',xrange=[0,6],yrange=[-0.3,1.5])
;  p=plot(plotdat.er,-plotdat.v,thick=2,color='orange',/overplot)
  ;=============
  IF peakplot EQ 1 THEN BEGIN
  plotdat.v = plotdat.v
  maxr = where(plotdat.v eq max(plotdat[where(plotdat.er LT 5.5d)].v))
  lx=indgen(17)*0 + plotdat[maxr].er
  ly=indgen(17)/10d - 0.3
  print,strcompress('peak rotation at'+string(plotdat[maxr].er)+'R_e')
  oplot,plotdat.er,plotdat.v,thick=2,color=fsc_color('lime green')
  oplot,lx,ly,linestyle=1,color=fsc_color('red')
  maxr = plotdat[maxr].er
  ENDIF
  oplot,lx,ly,linestyle=1,color=fsc_color('red')

;if discrete eq 1 then stop
     if discrete eq 1 then begin
      if l eq 0 then begin
        openw,lun,'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\N4564_08_11_2016_3.dat',/GET_LUN
        out=dblarr(3,n_elements(bin.er))
        out[0,*]=bin.er/reff
        out[1,*]=bin.quan/sige
        out[2,*]=l
        printf,lun,FORMAT='(3F)',out
      endif
      
      if l gt 0 then begin
        openw,lun,'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\N4564_08_11_2016_3.dat',/append, /GET_LUN
        out=dblarr(3,n_elements(bin.er))
        out[0,*]=bin.er/reff
        out[1,*]=bin.quan/sige
        out[2,*]=l
        printf,lun,FORMAT='(3F)',out
        close,lun
      endif
     endif
  ENDIF
  
;**********************************************Inner and Outer Annuli*************************************************
;seems like i get way better vrot fits when kpa is left as a free parameter


  IF twobin EQ 1 THEN BEGIN
    
    
    IF fullprofile EQ 0 THEN BEGIN
    
  ;=============MOVING WINDOW VARIABLES
    xsize = 8 
    wsize = 0.75;
    dx = 0.1
    incs = (xsize-wsize)/dx
    plotdat = fltarr(5,incs)

    ENDIF
    
    ;inner bin name
    bin4 = bin44
    ;outer bin name
    bin5 = bin55

    plotdatinv = fltarr(5,incs)

    i=0
  ;=============MOVING WINDOW ALGORITHM
for i=0,incs-1 do begin
   currentdatinv = bin4[where(bin4.er/reff GT i*dx and bin4.er/reff LT i*dx+wsize)]
   plotdatinv[0,i] = mean([0,wsize])+i*dx
   if n_elements(currentdatinv) GT 2 then begin
    plotdatinv[1,i] = mean(currentdatinv.quan/sige)
    plotdatinv[2,i] = median(currentdatinv.quan/sige)
    plotdatinv[3,i] = stddev(currentdatinv.quan/sige)
    plotdatinv[4,i] = median(currentdatinv.PA/!dtor)
   endif
endfor

    plotdatoutv = fltarr(5,incs)


i=0

for i=0,incs-1 do begin
   currentdatoutv = bin5[where(bin5.er/reff GT i*dx and bin5.er/reff LT i*dx+wsize)]
   plotdatoutv[0,i] = mean([0,wsize])+i*dx
   if n_elements(currentdatoutv) GT 2 then begin
    plotdatoutv[1,i] = mean(currentdatoutv.quan/sige)
    plotdatoutv[2,i] = median(currentdatoutv.quan/sige)
    plotdatoutv[3,i] = stddev(currentdatoutv.quan/sige)
    plotdatoutv[4,i] = median(currentdatoutv.PA/!dtor)
   endif
endfor
  ;=============
    
    

  endif

;*************************************************inner and outer bins (kPA fits)*********************************************************************



  if n_elements(bin4) GT 0 and n_elements(bin5) GT 0 then begin

    plotdatinp = fltarr(5,incs)

    i=0
  ;=============MOVING WINDOW ALGORITHM
for i=0,incs-1 do begin
   currentdatinp = bin4[where(bin4.er/reff GT i*dx and bin4.er/reff LT i*dx+wsize)]
   plotdatinp[0,i] = mean([0,wsize])+i*dx
   if n_elements(currentdatinp) GT 2 then begin
    plotdatinp[1,i] = mean(currentdatinp.quan/sige)
    plotdatinp[2,i] = median(currentdatinp.quan/sige)
    plotdatinp[3,i] = stddev(currentdatinp.quan/sige)
    plotdatinp[4,i] = median(currentdatinp.PA/!dtor)
   endif
endfor

    plotdatoutp = fltarr(5,incs)


i=0

for i=0,incs-1 do begin
   currentdatoutp = bin5[where(bin5.er/reff GT i*dx and bin5.er/reff LT i*dx+wsize)]
   plotdatoutp[0,i] = mean([0,wsize])+i*dx
   if n_elements(currentdatoutp) GT 2 then begin
    plotdatoutp[1,i] = mean(currentdatoutp.quan/sige)
    plotdatoutp[2,i] = median(currentdatoutp.quan/sige)
    plotdatoutp[3,i] = stddev(currentdatoutp.quan/sige)
    plotdatoutp[4,i] = median(currentdatoutp.PA/!dtor)
   endif
endfor
  ;=============
    
  ;=============Vrot PLOTS
  IF fullprofile EQ 0 THEN BEGIN
    linex=indgen(9)
    liney=indgen(9)*0
    p=plot(plotdatinv[0,*],plotdatinv[2,*],linestyle=0,thick=2,color='blue',xrange=[0,8],yrange=[-1.5,1.5],xtitle='R/R_e',ytitle='$v_{rot}/\sigma_e$')
    p=plot(linex,liney,linestyle=1,color='gray',/OVERPLOT)
    p=plot(plotdatoutv[0,*],plotdatoutv[2,*],linestyle=0,thick=2,color='blue',/overplot)
  ENDIF
  
  IF fullprofile EQ 1 AND twobin EQ 1 THEN BEGIN
   p=plot(plotdatinv[0,*],plotdatinv[2,*],linestyle=0,thick=2,color='blue',/overplot)
   p=plot(plotdatoutv[0,*],plotdatoutv[2,*],linestyle=0,thick=2,color='blue',/overplot)
  ENDIF
  ;=============

    
    
  ;=============kPA PLOTS
   IF fullprofile EQ 1 AND twobin EQ 1 THEN BEGIN

    p=plot(plotdat[0,*],plotdat[4,*],linestyle=0,thick=2,color='orange',xrange=[0,8],xtitle='R/R_e',ytitle='kPA',yrange=[-180,180])
    p=plot(plotdatinp[0,*],plotdatinp[4,*],linestyle=0,thick=2,color='red',/overplot)
    p=plot(plotdatoutp[0,*],plotdatoutp[4,*],linestyle=0,thick=2,color='red',/overplot)
   ENDIF
   
   IF fullprofile EQ 0 THEN BEGIN
    linex=indgen(9)
    liney=indgen(9)*0
    p=plot(plotdatinp[0,*],plotdatinp[4,*],linestyle=0,thick=2,color='red',xrange=[0,8],xtitle='R/R_e',ytitle='kPA',yrange=[-180,180],layout=[2,1,2],/current)
    p=plot(linex,liney,linestyle=1,color='grey',/OVERPLOT)
    p=plot(plotdatoutp[0,*],plotdatoutp[4,*],linestyle=0,thick=2,color='red',/overplot)
   ENDIF
  ;=============

endif
;**********************************************************************************************************************
;undefine,bin
END







