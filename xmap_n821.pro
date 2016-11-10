

pro xmap_n821, minerr, veldisp=veldisp, h3=h3, h4=h4, rms=rms, ps=ps, nosauron=nosauron, colorbar=colorbar, compass=compass

	xticknames_latex = strtrim([-6,-4,-2,0,2,4,6],1)
	yticknames_latex = strtrim([-4,-2,0,2,4],1)

  ngals = 2
  xsiz = 50
  range = 183

  xmap_prep, s, colorbar=colorbar, compass=compass, xsiz, ngals, range, nrows, ncols  ;for Brutus, lives in \cygwin directory

	s.smooth = round(0/100.*s.xsiz)
	s.smooth = 0
;	s.range = [-1,1]
	loc = 0

	kinType = 'VEL'
	if (n_elements(rms) gt 0) then kinType = 'RMS'
	if (n_elements(veldisp) gt 0) then kinType = 'VELDISP'
	if (n_elements(h3) gt 0) then kinType = 'H3'
	if (n_elements(h4) gt 0) then kinType = 'H4'


	yoff = 0.006


  vs_fitting=0

;  galname = 'ngc1344'
;  s.title = 'NGC '
;  s.range = [-1,1]*100
;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,0], nox=0, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, vs_fitting=vs_fitting, minerr=35d, /noextrap, /nogcs, /nosmeag, /nosauron




;  **For Two-Phase Assembly paper: SAURON,SKiMS,long-slit,PNe,RGCs. (EVERYTHING BUT BLUE GC's)**
;  **For SLUGGS survey paper, N821 should have SAURON, +long/SMEAGOL, +GCs
; ** V_sys is applied to Sauron/Skims/long in map_galaxy__define.pro. Discrete sets have offset applied in map_ngc****__define.pro!!**


;   galname = 'NGC4697'
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, $
;    doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nosauron, /nolong, /nosmeag;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne

;    galnames = ['NGC0821', 'NGC1023', 'NGC1344', 'NGC2768', 'NGC3115', 'NGC3377', 'NGC4564', 'NGC4697']
    galnames = ['NGC0821', 'NGC1023', 'NGC2768', 'NGC3115', 'NGC3377', 'NGC4564', 'NGC4697']
    
    hold = replicate({iv:0d, div:0d, ov:0d, dov:0d, ipa:0d, dipa:0d, opa:0d, dopa:0d, d:0d, dd:0d, g:0d, dg:0d, vdisp:0d},1)
    
    N0821 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N1023 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N1344 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N2768 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N3115 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N3377 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N4473 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N4564 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    N4697 = replicate({plotx:fltarr(75),ploty:fltarr(75),vradi:fltarr(1000),pai:fltarr(1000),dvradi:fltarr(1000),vroti:0d,vsysi:0d,kpai:0d,kqi:0d,vrado:fltarr(1000),pao:fltarr(1000),dvrado:fltarr(1000),vroto:0d, vdisp:0d,vsyso:0d,kpao:0d,kqo:0d,sige:0d,reff:0d},1)
    
    ;Whether or not to use v_systemic from kinemetry
    vs_fitting = 0
    
    
;    openw,1,'F:\obs_10_29_2016.dat'
    
    FOR i = 0,0 DO BEGIN;n_elements(galnames)-1 DO BEGIN
      i=5
      galname = galnames[i]
      
      l=i
      
      print,galname
      
;      s.range = [-1,1]*100
;      map_galaxy, galname, s, p2=p2, pos=s.pos[*,0], nox=0, noy=0, doreff=s.doreff, isophote=1, norm=1., /current, /new, maxPoints=1d2, maxDist=maxDist, minerr=15d,vs_fitting=0, /noextrap;,/nolong;, /nodiscrete;, /nosauron, /nolong
;      stop
;      IF galname EQ 'NGC0821' OR galname EQ 'NGC1023' THEN BEGIN
;       kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, $
;        doreff=s.doreff, norm=1., innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
;        maxPoints=1d5, maxDist=maxDist, minerr=15d, N0821, N1023, N2768, N3115, N3377, N4564, N4697, $ 
;        /noextrap, /isophote, /current, /new, /nodiscrete
;      ENDIF
;      IF galname EQ 'NGC4697' THEN BEGIN
;       kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, $
;        doreff=s.doreff, norm=1., innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
;        maxPoints=1d5, maxDist=maxDist, minerr=15d, N0821, N1023, N2768, N3115, N3377, N4564, N4697, $ 
;        /noextrap, /isophote, /current, /new, /nolong, /nodiscrete, /nosmeag
;      ENDIF ELSE BEGIN

nbootstraps = 1
peakplot = 0
twobin =0
fullprofile = 1
        kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, $
        doreff=s.doreff, norm=1., innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
        maxPoints=1d5, maxDist=maxDist, minerr=15d, N0821, N1023, N1344, N2768, N3115, N3377, N4564, N4697, N4473, inplotdat, inmaxr, $
        bin1, i, /noextrap, /isophote, /current, /new, /nodiscrete, nbootstraps, twobin=twobin, obser, peakplot=peakplot, qin, dqin, nbootstrapsq
        
      hold.iv = innerv
      hold.div = dinnerv  
      hold.ipa = innerpa
      hold.dipa = dinnerpa
      
        kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, $
        doreff=s.doreff, norm=1., innerv, outerv, dinnerv, douterv, innerpa, outerpa, dinnerpa, douterpa, vdisp, dvdisp, vs_fitting, $
        maxPoints=1d5, maxDist=maxDist, minerr=15d, N0821, N1023, N1344, N2768, N3115, N3377, N4564, N4697, N4473, inplotdat, inmaxr, $
        bin1, i, /noextrap, /isophote, /current, /new, /nostellar, /nolong, /nosauron, /nosmeag, nbootstraps, twobin=twobin, obser, $
        peakplot=peakplot, qin, dqin, nbootstrapsq
        
      hold.ov = outerv
      hold.dov = douterv
      hold.opa = outerpa
      hold.dopa = douterpa
      
      hold.d = hold.ov - hold.iv
      hold.dd = sqrt(hold.div^2d + hold.dov^2d)
;      if hold.ipa gt !pi/2d and hold.opa lt -!pi/2d then hold.g = difference([-!pi,hold.opa]) - difference([!pi,hold.ipa])
;      if hold.ipa lt -!pi/2d and hold.opa gt !pi/2d then hold.g = difference([!pi,hold.opa]) - difference([-!pi,hold.ipa])
      hold.g = hold.opa - hold.ipa
      hold.g = abs(hold.g)
      if hold.g GT 180d then hold.g = 360d - hold.g
      hold.dg = sqrt(hold.dopa^2d)

      hold.vdisp = vdisp
      IF i eq 0 THEN BEGIN  
      kinresult = hold
      ENDIF ELSE BEGIN
      kinresult = struct_append(kinresult, hold)
      ENDELSE
      
      out = dblarr(8,1)
      out[0,*] = hold.iv
      out[1,*] = hold.div
      out[2,*] = hold.ov
      out[3,*] = hold.dov
      out[4,*] = hold.d
      out[5,*] = hold.dd
      out[6,*] = hold.g
      out[7,*] = hold.dg
      
      out2 = dblarr(2,75)
      out2[0,*] = N1023.plotx
      out2[1,*] = N1023.ploty
      
      
;      openu,2+i,'F:\cygwin64\home\Justin\Exelis\IDL82\lib\SIMDATA\obs_07_19_2016_3.dat',/append
;      printf,2+i,FORMAT='(8F)',out
;      close,2+i
      
;      openu,2+1,'F:\obs_10_29_2016.dat',/append
;      printf,2+1,out2
;      close,2+1
    ENDFOR

    linex=indgen(7)
    liney=indgen(7)*0
;;    pf1 = poly_fit(N0821.plotx,N0821.ploty,5,yfit=yfit1)
;;    pf2 = poly_fit(N1023.plotx,N1023.ploty,5,yfit=yfit2)
;;    pf3 = poly_fit(N2768.plotx,N2768.ploty,5,yfit=yfit3)
;;    pf4 = poly_fit(N3115.plotx,N3115.ploty,5,yfit=yfit4)
;    pf5 = poly_fit(N3377.plotx,N3377.ploty,5,yfit=yfit5)
; stop
;    pf6 = poly_fit(N4564.plotx,N4564.ploty,5,yfit=yfit6)
;    pf7 = poly_fit(N4697.plotx,N4697.ploty,5,yfit=yfit7)
;
;    p=plot(N0821.plotx[where(N0821.plotx NE 0)],N0821.ploty[where(N0821.plotx NE 0)],thick=2,color='olive',xrange=[0,6],yrange=[-0.3,1.4],xtitle='$R/R_e$',ytitle='$v_{rot}/\sigma_e$')
;    p=plot(linex,liney,linestyle=1,color='grey',/overplot)
;;
;    p=plot(N1023.plotx[where(N1023.plotx NE 0)],N1023.ploty[where(N1023.plotx NE 0)],thick=2,color='dark goldenrod',/overplot)
;    p=plot(N2768.plotx[where(N2768.plotx NE 0)],N2768.ploty[where(N2768.plotx NE 0)],thick=2,color='gold',/overplot)
;    p=plot(N3115.plotx[where(N3115.plotx NE 0)],N3115.ploty[where(N3115.plotx NE 0)],thick=2,color='goldenrod',/overplot)
;    p=plot(N3377.plotx[where(N3377.plotx NE 0)],N3377.ploty[where(N3377.plotx NE 0)],thick=2,color='dark orange',/overplot)
;    p=plot(N4697.plotx[where(N4697.plotx NE 0)],N4697.ploty[where(N4697.plotx NE 0)],thick=2,color='orange',/overplot)
;    p=plot(N4564.plotx[where(N4564.plotx NE 0)],N4564.ploty[where(N4564.plotx NE 0)],thick=2,color='yellow green',/overplot)
stop

;kinemetry_plots, N0821, N1023, N2768, N3115, N3377, N4473, N4564, N4697, kinresult
;
;a=indgen(17)/10.;x-coordinates of 1-to-1 relation line and hori. line
;b=indgen(17)*0. ;hori. line y-coords
;c=a ;1-to-1 line
;p=errorplot(kinresult.d, abs(kinresult.g), kinresult.dd, kinresult.dg, linestyle=6,symbol='o',sym_filled=1,sym_size=0.9,xrange=[-1,1],yrange=[0,180],'black')
;p=errorplot(kinresult.iv,kinresult.ov,kinresult.div,kinresult.dov,linestyle=6,symbol='o',sym_filled=1,sym_size=0.9,xrange=[0,1.6],yrange=[-0.3,1.6],color='orange')
;p=plot(a,c,linestyle=3,color='black',/overplot)





;   galname = 'ngc0821'
;  s.title = 'NGC 821'
;;  s.title = ''
;;; s.range = [-69,66];[-1,1]*69.;*71.65
;  s.range = [-1,1]*90.
; map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap,/nodiscrete,/nosmeag,/nolong;, /nodiscrete;, /medianing;,/vmax
;;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, $
;;    doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nosauron, /nolong, /nosmeag;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;stop
;;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar;, /red;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;;++loc
;;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap,/nodiscrete;, /nodiscrete;, /medianing;,/vmax
;
;;++loc
;;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=0, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap, /nopne;, /nosauron, /nolong, /nosmeag;, /nopne;,/nodiscrete,/nosmeag,/nolong;, /nodiscrete;, /medianing;,/vmax
;;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=0, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap,/nodiscrete,/nosmeag,/nolong;, /nodiscrete;, /medianing;,/vmax
;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,0], nox=1, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap;, /nosauron, /nolong, /nodiscrete
;  p2 = text(s.pos[2,0]-0.395, s.pos[1,0] + 0.19, '$\bfNGC 821$', /overplot, FONT_SIZE=11)
;
;;  xoff = -0.07
;	++loc
;  galname = 'ngc1023'
;;;;	
;;  s.title = 'NGC 1023'
;; s.title = ''
;  s.range = [-1,1]*200.
;    map_galaxy, galname, s, p2=p2, pos=s.pos[*,2], nox=1, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap;, /nosauron, /nolong, /nodiscrete
;  p2 = text(s.pos[2,2]-0.395, s.pos[1,2] + 0.19, '$\bfNGC 1023$', /overplot, FONT_SIZE=11)
;
;
;;	map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote;,/nosmeag,/nolong,/nodiscrete;, /nodiscrete
;;	xoff = -0.07
;	;p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, $
;    doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nodiscrete;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;stop



;;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar;, /red;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;
;;;;  xoff = -0.07
;;;;  ++loc
;  ++loc
;  galname = 'ngc2768'
;;;  s.title = 'NGC 2768'
;; s.title = ''
;;  s.range = [-1,1]*150.
;s.range = [-1,1]*200.
;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,4], nox=1, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap;, /nosauron, /nolong, /nodiscrete
;  p2 = text(s.pos[2,4]-0.395, s.pos[1,4] + 0.19, '$\bfNGC 2768$', /overplot, FONT_SIZE=11)
;
;;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote;,/nosmeag,/nolong,/nodiscrete
;;  xoff = -0.07
; ; p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote;, /nodiscrete, /nosauron, /nolong;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar;, /red;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;
;  ++loc
;  galname = 'ngc3115' 
;;;  
;;;  s.title = 'NGC 3115
;;   s.title = ''a
;;;;;  s.range = [-1,1]*200
;  s.range = [-1,1]*280
;    map_galaxy, galname, s, p2=p2, pos=s.pos[*,6], nox=0, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap, /nosauron;, /nolong, /nodiscrete
;    p2 = text(s.pos[2,6]-0.395, s.pos[1,6] + 0.19, '$\bfNGC 3115$', /overplot, FONT_SIZE=11)
  
;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nosauron;,/nolong,/nodiscrete
;  xoff = -0.07
;  p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
; kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote;, /nolong, /nosauron, /nodiscrete;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar;, /red;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne

;  ++loc
;  galname = 'ngc3377'  
;  
;  s.title = 'NGC 3377'
;;   s.title = ''
;;;  s.range = [-1,1]*90.
;;
;  s.range = [-1,1]*110
;    map_galaxy, galname, s, p2=p2, pos=s.pos[*,1], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap;, /nosmeag,/nosauron,/nodiscrete
;    p2 = text(s.pos[2,1]-0.395, s.pos[1,1] + 0.19, '$\bfNGC 3377$', /overplot, FONT_SIZE=11)
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote;, /nodiscrete;, /nosauron, /nolong;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne

;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=0, noy=0, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote;,/nosmeag,/nolong,/nodiscrete
;  xoff = -0.07
 ; p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nolong, /nosauron, /nodiscrete;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar, /red;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne


;++loc
;  galname = 'ngc4473'  
;  s.title = 'NGC 4473'
;  s.range = [-1,1]*100;75
;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nopne,  /green
;  xoff = -0.07
;  p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nolong, /nosauron, /nodiscrete;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar;, /red;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne

;++loc
;
;  galname = 'NGC4564' 
;;  
;;  s.title = 'NGC 4564'
;; s.title = ''
;  s.range = [-150,180] ;[-1,1]*100.
;    map_galaxy, galname, s, p2=p2, pos=s.pos[*,3], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap; /nodiscrete
;    p2 = text(s.pos[2,3]-0.395, s.pos[1,3] + 0.19, '$\bfNGC 4564$', /overplot, FONT_SIZE=11)

;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=0, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nogcs;,/nosmeag,/nolong,/nodiscrete
;  xoff = -0.07
;  p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote;, /nolong, /nosauron, /nodiscrete;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, p2=p2, innervel=innervel, innerpa=innerpa, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar;, /red;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne

;stop
;++loc
;
;  galname = 'ngc4697' 
;;  
;;  s.title = 'NGC 4697'
; s.title = ''
;  s.range = [-1
;  ,1]*140.
;;;;  s.range = [-1,1]*89.1
;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,5], nox=0, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d10, maxDist=maxDist, minerr=15d, /noextrap, /nolong;, /nosauron, /nolong
;  p2 = text(s.pos[2,5]-0.395, s.pos[1,5] + 0.19, '$\bfNGC 4697$', /overplot, FONT_SIZE=11)

;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=0, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /nolong, /isophote, /nogcs;,/nosmeag,/nodiscrete
;  xoff = -0.07
 ; p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, innervel=innervel, innerpa=innerpa, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nodiscrete, /nosauron, /nolong;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne
;   kinemetry_skims, galname, s, pp1=pp1, thisdevice=thisdevice, innervel=innervel, innerpa=innerpa, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /isophote, /nostellar;, /nogcs, /nodiscrete;, /diagnostics;,/nolong;, /nodiscrete;, /vmax;/nopne; /nodiscrete, /nolong, /nosmeag; /red, /nopne

;;galname = 'ngc4526'
;;
;;   s.title = 'NGC 4526'
;;  s.range = [-1.2,0.5]
;;  map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=15d, /noextrap, /nolong, /nodiscrete
;;  xoff = -0.07
;;  p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10, /fill_background, fill_color='white')
;  
;
;  
;
;;do the discrete tracers not have vel dispersions recorded?
;;	++loc
;;	s.range = [100,200]
;;	map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=10d, /noextrap, /veldisp, /nodiscrete, nosauron=nosauron
;;	xoff = -0.08
;;	p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0])+'/'+roundx(s.range[1]), /overplot, FONT_SIZE=10)
;;
;;	++loc
;;	s.title = ''
;;	map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /isophote, /dataloc, /fill, /current, nosauron=nosauron , /nodiscrete
;;
;;
;;	++loc
;;	minerr_h3 = 0.05d
;;	s.range = [-1,1]*0.1
;;	map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=[0,1], noy=[0,1], doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=minerr_h3, /noextrap, /h3, nosauron=nosauron, /nodiscrete
;;	xoff = -0.10
;;	p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0],2)+'/'+roundx(s.range[1],2), /overplot, FONT_SIZE=10)
;;
;;
;;	++loc
;;	minerr_h4 = 0.05d	
;;	s.range = [-0.10,0.10]
;
;	map_galaxy, galname, s, p2=p2, pos=s.pos[*,loc], nox=1, noy=1, doreff=s.doreff, norm=1., /current, /new, maxPoints=1d5, maxDist=maxDist, minerr=minerr_h4, /noextrap, /h4, /nodiscrete ; nosauron=nosauron, 
;	xoff = -0.10
;	p2 = text(s.pos[2,loc] + xoff, s.pos[1,loc] + yoff, roundx(s.range[0],2)+'/'+roundx(s.range[1],2), /overplot, FONT_SIZE=10, target = p2)

;	if (n_elements(compass) gt 0) then begin
;		++loc
;		xmap_compass, p2, s, s.pos[*,loc], head_indent=0.3, head_angle=30, minor=0.15, thick=thick, rot=s.parot
;	endif
;	
;
;s = struct_addtags(s, {colorbar_height:300d, colorbar_thick:15d})
;
;	if (n_elements(colorbar) eq 0) then begin
;		++loc
;		xmap_colorbar, p2, s;, /hori;, /horizontal, /notags ;commented out "Sharpcorner" from xmap_colorbar.pro
;	endif
;
;
;	outfile = 'smeagol_'+kinType+'_'+galname+'_maps_v'
;		file = file_search(outfile+'*',count=count)
;		num = 0
;		if (count gt 0) then for i=0,(count - 1) do num = [num, fix((strsplit((strsplit(file[i],outfile,/extract,/regex))[-1],'.eps',/extract))[0])]
;		p2.save, outFile + strtrim(max(num)+1,2) + '.eps';, /cmyk
;	
;	
;	sc = 1.0
;	latexify, outFile + strtrim(max(num)+1,2) + '.eps', [s.xtitle,s.ytitle,s.xticknames, s.yticknames], $
;		['$\Sigma$',s.ytitle_latex,s.xticknames_latex, s.yticknames_latex], $
;		[[1,1]*sc,replicate(sc,n_elements(s.xticknames)+n_elements(s.yticknames))]

;	p2.close
;	obj_destroy, p2

end


