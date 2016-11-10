PRO kinemetry_plots, N0821, N1023, N2768, N3115, N3377, N4473, N4564, N4697, kinresult, ngals, qplot=qplot

if ngals eq 1 then begin
  nrows = 1
  ncols = 2
endif else begin
  if ngals mod 2 eq 0 then nrows = ngals/2d
  if ngals mod 2 ne 0 then nrows = (ngals+1)/2d
  if ngals eq 7 then nrows = 4
  ncols = 4
endelse 

if qplot eq 1 then begin
  nrows = 4
  ncols = 2
endif



j=0      
if n_elements(N0821) gt 0 then begin
      j=j+1
      ;inner kinemetry plot
      N0821.vroti = kinresult[0].iv*N0821.sige
      N0821.kpai = kinresult[0].ipa*!dtor
      N0821.vroto = kinresult[0].ov*N0821.sige
      N0821.kpao = kinresult[0].opa*!dtor
      N0821.vdisp = kinresult[0].vdisp*N0821.sige
      result821_in = [N0821.vroti,N0821.kpai,N0821.kqi,N0821.vsysi]
      result821_out = [N0821.vroto,N0821.kpao,N0821.kqo,N0821.vsyso]
      vmod821_in = kinem_func(N0821.pai, [result821_in[0],result821_in[1],result821_in[2],result821_in[3]])
      vmod821_out = kinem_func(N0821.pao, [result821_out[0],result821_out[1],result821_in[2],result821_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
      p=errorplot(N0821.pai,N0821.vradi,N0821.dvradi,linestyle=6,symbol='o',sym_size=0.25,sym_filled=1,xrange=[0,2*!pi],title='NGC 821',yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi])
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N0821.pai[sort(N0821.pai)],vmod821_in[sort(N0821.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[0].iv*N0821.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[0].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[0].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data) 
      if qplot eq 1 then t=text(3,200,strcompress('$kQ_{in} =$'+string(N0821.kqi)),target=pp,font_size=5,color='dodger blue',/data)    
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin  
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N0821.pao,N0821.vrado,N0821.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N0821.pao[sort(N0821.pao)],vmod821_out[sort(N0821.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N0821.pao[sort(N0821.pao)],vmod821_out[sort(N0821.pao)]+N0821.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N0821.pao[sort(N0821.pao)],vmod821_out[sort(N0821.pao)]-N0821.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[0].ov*N0821.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[0].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[0].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[0].vdisp*N0821.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      ;------------------------------
    endif
    
    
    
endif

if n_elements(N1023) gt 0 then begin
      j = j+1
      
      
      ;inner kinemetry plot
      N1023.vroti = kinresult[1].iv*N1023.sige
      N1023.kpai = kinresult[1].ipa*!dtor
      N1023.vroto = kinresult[1].ov*N1023.sige
      N1023.kpao = kinresult[1].opa*!dtor
      N1023.vdisp = kinresult[1].vdisp*N1023.sige
      result1023_in = [N1023.vroti,N1023.kpai,N1023.kqi,N1023.vsysi]
      result1023_out = [N1023.vroto,N1023.kpao,N1023.kqo,N1023.vsyso]
      vmod1023_in = kinem_func(N1023.pai, [result1023_in[0],result1023_in[1],result1023_in[2],result1023_in[3]])
      vmod1023_out = kinem_func(N1023.pao, [result1023_out[0],result1023_out[1],result1023_in[2],result1023_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
      p=errorplot(N1023.pai,N1023.vradi,N1023.dvradi,linestyle=6,symbol='o',sym_size=0.25,sym_filled=1,title='NGC 1023',xrange=[0,2*!pi],yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi],/current)
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N1023.pai[sort(N1023.pai)],vmod1023_in[sort(N1023.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[1].iv*N1023.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[1].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[1].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data)    
      if qplot eq 1 then t=text(3,200,strcompress('$kQ_{in} =$'+string(N1023.kqi)),target=pp,font_size=5,color='dodger blue',/data)    
       
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N1023.pao,N1023.vrado,N1023.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N1023.pao[sort(N1023.pao)],vmod1023_out[sort(N1023.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N1023.pao[sort(N1023.pao)],vmod1023_out[sort(N1023.pao)]+N1023.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N1023.pao[sort(N1023.pao)],vmod1023_out[sort(N1023.pao)]-N1023.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[1].ov*N1023.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[1].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[1].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[1].vdisp*N1023.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      ;------------------------------
    endif
    
    
endif

if n_elements(N2768) gt 0 then begin
      j=j+1
    ;inner kinemetry plot
      N2768.vroti = kinresult[2].iv*N2768.sige
      N2768.kpai = kinresult[2].ipa*!dtor
      N2768.vroto = kinresult[2].ov*N2768.sige
      N2768.kpao = kinresult[2].opa*!dtor
      N2768.vdisp = kinresult[2].vdisp*N2768.sige
      result2768_in = [N2768.vroti,N2768.kpai,N2768.kqi,N2768.vsysi]
      result2768_out = [N2768.vroto,N2768.kpao,N2768.kqo,N2768.vsyso]
      vmod2768_in = kinem_func(N2768.pai, [result2768_in[0],result2768_in[1],result2768_in[2],result2768_in[3]])
      vmod2768_out = kinem_func(N2768.pao, [result2768_out[0],result2768_out[1],result2768_in[2],result2768_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
      p=errorplot(N2768.pai,N2768.vradi,N2768.dvradi,linestyle=6,symbol='o',sym_size=0.25,title='NGC 2768',sym_filled=1,xrange=[0,2*!pi],yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi],/current)
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N2768.pai[sort(N2768.pai)],vmod2768_in[sort(N2768.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[2].iv*N2768.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[2].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[2].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data)
      if qplot eq 1 then t=text(3,200,strcompress('$kQ_{in} =$'+string(N2768.kqi)),target=pp,font_size=5,color='dodger blue',/data)      
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin      
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N2768.pao,N2768.vrado,N2768.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N2768.pao[sort(N2768.pao)],vmod2768_out[sort(N2768.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N2768.pao[sort(N2768.pao)],vmod2768_out[sort(N2768.pao)]+N2768.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N2768.pao[sort(N2768.pao)],vmod2768_out[sort(N2768.pao)]-N2768.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[2].ov*N2768.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[2].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[2].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[2].vdisp*N2768.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      ;------------------------------
    endif
    
    
endif

if n_elements(N3115) gt 0 then begin
      
      j=j+1
      ;inner kinemetry plot
      N3115.vroti = kinresult[3].iv*N3115.sige
      N3115.kpai = kinresult[3].ipa*!dtor
      N3115.vroto = kinresult[3].ov*N3115.sige
      N3115.kpao = kinresult[3].opa*!dtor
      N3115.vdisp = kinresult[3].vdisp*N3115.sige
      result3115_in = [N3115.vroti,N3115.kpai,N3115.kqi,N3115.vsysi]
      result3115_out = [N3115.vroto,N3115.kpao,N3115.kqo,N3115.vsyso]
      vmod3115_in = kinem_func(N3115.pai, [result3115_in[0],result3115_in[1],result3115_in[2],result3115_in[3]])
      vmod3115_out = kinem_func(N3115.pao, [result3115_out[0],result3115_out[1],result3115_in[2],result3115_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
      
      p=errorplot(N3115.pai,N3115.vradi,N3115.dvradi,linestyle=6,symbol='o',sym_size=0.25,title='NGC 3115',sym_filled=1,xrange=[0,2*!pi],yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi],/current)
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N3115.pai[sort(N3115.pai)],vmod3115_in[sort(N3115.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[3].iv*N3115.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[3].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[3].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data)   
      if qplot eq 1 then t=text(3,200,strcompress('$kQ_{in} =$'+string(N3115.kqi)),target=pp,font_size=5,color='dodger blue',/data)   
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin      
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N3115.pao,N3115.vrado,N3115.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N3115.pao[sort(N3115.pao)],vmod3115_out[sort(N3115.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N3115.pao[sort(N3115.pao)],vmod3115_out[sort(N3115.pao)]+N3115.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N3115.pao[sort(N3115.pao)],vmod3115_out[sort(N3115.pao)]-N3115.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[3].ov*N3115.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[3].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[3].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[3].vdisp*N3115.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      ;------------------------------
    endif
    
    
    
endif

if n_elements(N3377) gt 0 then begin
      
      j=j+1
      ;inner kinemetry plot
      N3377.vroti = kinresult[4].iv*N3377.sige
      N3377.kpai = kinresult[4].ipa*!dtor
      N3377.vroto = kinresult[4].ov*N3377.sige
      N3377.kpao = kinresult[4].opa*!dtor
      N3377.vdisp = kinresult[4].vdisp*N3377.sige
      result3377_in = [N3377.vroti,N3377.kpai,N3377.kqi,N3377.vsysi]
      result3377_out = [N3377.vroto,N3377.kpao,N3377.kqo,N3377.vsyso]
      vmod3377_in = kinem_func(N3377.pai, [result3377_in[0],result3377_in[1],result3377_in[2],result3377_in[3]])
      vmod3377_out = kinem_func(N3377.pao, [result3377_out[0],result3377_out[1],result3377_in[2],result3377_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
            
      p=errorplot(N3377.pai,N3377.vradi,N3377.dvradi,linestyle=6,symbol='o',sym_size=0.25,title='NGC 3377',sym_filled=1,xrange=[0,2*!pi],yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi],/current)
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N3377.pai[sort(N3377.pai)],vmod3377_in[sort(N3377.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[4].iv*N3377.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[4].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[4].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data)     
      if qplot eq 1 then t=text(3,200,strcompress('$kQ_{in} =$'+string(N3377.kqi)),target=pp,font_size=5,color='dodger blue',/data) 
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin      
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N3377.pao,N3377.vrado,N3377.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N3377.pao[sort(N3377.pao)],vmod3377_out[sort(N3377.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N3377.pao[sort(N3377.pao)],vmod3377_out[sort(N3377.pao)]+N3377.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N3377.pao[sort(N3377.pao)],vmod3377_out[sort(N3377.pao)]-N3377.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[4].ov*N3377.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[4].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[4].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[4].vdisp*N3377.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      ;------------------------------
    endif
    
    
endif

undefine,N4473
if n_elements(N4473) gt 0 then begin
  
      j=j+1
      ;inner kinemetry plot
      N4473.vroti = kinresult[4].iv*N4473.sige
      N4473.kpai = kinresult[4].ipa*!dtor
      N4473.vroto = kinresult[4].ov*N4473.sige
      N4473.kpao = kinresult[4].opa*!dtor
      N4473.vdisp = kinresult[4].vdisp*N4473.sige
      result4473_in = [N4473.vroti,N4473.kpai,N4473.kqi,N4473.vsysi]
      result4473_out = [N4473.vroto,N4473.kpao,N4473.kqo,N4473.vsyso]
      vmod4473_in = kinem_func(N4473.pai, [result4473_in[0],result4473_in[1],result4473_in[2],result4473_in[3]])
      vmod4473_out = kinem_func(N4473.pao, [result4473_out[0],result4473_out[1],result4473_in[2],result4473_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
            
      p=errorplot(N4473.pai,N4473.vradi,N4473.dvradi,linestyle=6,symbol='o',sym_size=0.25,title='NGC 4473',sym_filled=1,xrange=[0,2*!pi],yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi],/current)
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N4473.pai[sort(N4473.pai)],vmod4473_in[sort(N4473.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[4].iv*N4473.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[4].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[4].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data)    
      if qplot eq 1 then t=text(3,200,strcompress('$kQ_{in} =$'+string(N4564.kqi)),target=pp,font_size=5,color='dodger blue',/data)  
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin      
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N4473.pao,N4473.vrado,N4473.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N4473.pao[sort(N4473.pao)],vmod4473_out[sort(N4473.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N4473.pao[sort(N4473.pao)],vmod4473_out[sort(N4473.pao)]+N4473.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N4473.pao[sort(N4473.pao)],vmod4473_out[sort(N4473.pao)]-N4473.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[4].ov*N4473.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[4].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[4].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[4].vdisp*N4473.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      ;------------------------------
    endif
    
    
endif

if n_elements(N4564) gt 0 then begin
      
      j=j+1
      ;inner kinemetry plot
      N4564.vroti = kinresult[5].iv*N4564.sige
      N4564.kpai = kinresult[5].ipa*!dtor
      N4564.vroto = kinresult[5].ov*N4564.sige
      N4564.kpao = kinresult[5].opa*!dtor
      N4564.vdisp = kinresult[5].vdisp*N4564.sige
      result4564_in = [N4564.vroti,N4564.kpai,N4564.kqi,N4564.vsysi]
      result4564_out = [N4564.vroto,N4564.kpao,N4564.kqo,N4564.vsyso]
      vmod4564_in = kinem_func(N4564.pai, [result4564_in[0],result4564_in[1],result4564_in[2],result4564_in[3]])
      vmod4564_out = kinem_func(N4564.pao, [result4564_out[0],result4564_out[1],result4564_in[2],result4564_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
            
      p=errorplot(N4564.pai,N4564.vradi,N4564.dvradi,linestyle=6,symbol='o',sym_size=0.25,title='NGC 4564',sym_filled=1,xrange=[0,2*!pi],yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi],/current)
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N4564.pai[sort(N4564.pai)],vmod4564_in[sort(N4564.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[5].iv*N4564.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[5].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[5].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data)     
      if qplot eq 1 then t=text(3,200,strcompress('$kQ_{in} =$'+string(N4564.kqi)),target=pp,font_size=5,color='dodger blue',/data) 
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin      
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N4564.pao,N4564.vrado,N4564.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N4564.pao[sort(N4564.pao)],vmod4564_out[sort(N4564.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N4564.pao[sort(N4564.pao)],vmod4564_out[sort(N4564.pao)]+N4564.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N4564.pao[sort(N4564.pao)],vmod4564_out[sort(N4564.pao)]-N4564.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[5].ov*N4564.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[5].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[5].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[5].vdisp*N4564.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      ;------------------------------
    endif
    
    
endif

if n_elements(N4697) gt 0 then begin
        
      
      j=j+1
      ;inner kinemetry plot
      N4697.vroti = kinresult[6].iv*N4697.sige
      N4697.kpai = kinresult[6].ipa*!dtor
      N4697.vroto = kinresult[6].ov*N4697.sige
      N4697.kpao = kinresult[6].opa*!dtor
      N4697.vdisp = kinresult[6].vdisp*N4697.sige
      result4697_in = [N4697.vroti,N4697.kpai,N4697.kqi,N4697.vsysi]
      result4697_out = [N4697.vroto,N4697.kpao,N4697.kqo,N4697.vsyso]
      vmod4697_in = kinem_func(N4697.pai, [result4697_in[0],result4697_in[1],result4697_in[2],result4697_in[3]])
      vmod4697_out = kinem_func(N4697.pao, [result4697_out[0],result4697_out[1],result4697_in[2],result4697_out[3]])
      
      linex=indgen(7)
      liney=indgen(7)*0.
      if qplot eq 0 then posi = 2*j-1
      if qplot eq 1 then posi = j
      
            
      p=errorplot(N4697.pai,N4697.vradi,N4697.dvradi,linestyle=6,symbol='o',sym_size=0.25,title='NGC 4697',sym_filled=1,xrange=[0,2*!pi],yrange=[-400,400],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',LAYOUT=[ncols,nrows,posi],/current)
      pp=plot(linex,liney,linestyle=2,color='black',/overplot)
      pp=plot(N4697.pai[sort(N4697.pai)],vmod4697_in[sort(N4697.pai)],color='dodger blue',linestyle=0,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[6].iv*N4697.sige))+'km/s'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[6].iv)),color='dodger blue',target=pp,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[6].ipa))+'degrees'),color='dodger blue',target=pp,font_size=5,/data)
      t=text(0.25,300,'~1 R_e',target=pp,font_size=9,/data)     
      if qplot eq 1 then t=text(0.25,200,strcompress('$kQ_{in} =$'+string(N4697.kqi)),target=pp,font_size=5,color='dodger blue',/data) 
      ;------------------------------
      
      ;outer kinemetry plot
    if qplot eq 0 then begin      
      linex=indgen(7)
      liney=indgen(7)*0.
      p=errorplot(N4697.pao,N4697.vrado,N4697.dvrado,linestyle=6,symbol='o',sym_filled=1,yrange=[-400,400],xrange=[0,6.3],xtitle='P.A. [rad]',ytitle='$V_{LOS} [km/s]$',sym_size=0.25,LAYOUT=[ncols,nrows,2*j],/current)
      p=plot(linex,liney,linestyle=2,color='black',/overplot)
      p=plot(N4697.pao[sort(N4697.pao)],vmod4697_out[sort(N4697.pao)],color='red',linestyle=0,thick=2,/overplot)
      p=plot(N4697.pao[sort(N4697.pao)],vmod4697_out[sort(N4697.pao)]+N4697.vdisp,color='red',linestyle=1,thick=2,/overplot)
      p=plot(N4697.pao[sort(N4697.pao)],vmod4697_out[sort(N4697.pao)]-N4697.vdisp,color='red',linestyle=1,thick=2,/overplot)
      t=text(3,350.,strcompress('Vrot='+string(round(kinresult[6].ov*N4697.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(3,300.,strcompress('V/sig='+string(kinresult[6].ov)),color='orange red',target=p,font_size=5,/data)
      t=text(3,250.,strcompress('kPA='+string(round(kinresult[6].opa))+'degrees'),color='orange red',target=p,font_size=5,/data)
      t=text(3,200.,strcompress('Vdisp='+string(round(kinresult[6].vdisp*N4697.sige))+'km/s'),color='orange red',target=p,font_size=5,/data)
      t=text(0.25,300,'~4 R_e',target=p,font_size=9,/data)     
      
      ;------------------------------
    endif 
     
endif    

END