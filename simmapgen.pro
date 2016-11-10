


;pro simmapgen, obj, p2=p2, min_points=min_points, method=method, c_thick=c_thick, tickfont_size=tickfont_size, pos=pos, xsiz=xsiz, $
;	current=current, xrange=xrange, yrange=yrange, nox=nox, noy=noy, xticknames=xticknames, yticknames=yticknames, range=range, $
;	smooth=smooth, cont=cont, c_value=c_value, xtitle=xtitle, ytitle=ytitle, rnorm=rnorm, dimensions=dimensions, reff=reff, $
;	xticklen=xticklen, smcon=smcon, yticklen=yticklen

pro simmapgen, obj, s, p2=p2, pos=pos, nox=nox, noy=noy, reff=reff, current=current, cont=cont

;sauron_colormap, rgb_table=rgb_table
;device, decomposed=0






vals = create_struct( $
		'xrange'		, max(abs(obj.x))*[-1,1], $
		'yrange'		, max(abs(obj.y))*[-1,1], $
		'xsiz'			, 300		, $
		'rnorm'			, 1.		, $
		'xticknames'	, ''		, $
		'yticknames'	, ''		, $
		'min_points'	, 15		, $
		'method'		, 'kriging'	, $
		'smooth'		, 0			, $
		'yticklen'		, 0.015		, $
		'xticklen'		, 0.015		, $
		'smcon'			, 0.		, $
		'xminor'		, 1			, $
		'yminor'		, 1			, $
		'isophote_level', [1,4]		, $
		'tickfont_size'	, 15		  $
	)

if (n_elements(reff) eq 0) then reff = 1.
if tag_exist(s, 'dimensions') then dimensions = s.dimensions
valtags = tag_names(vals)
for i=0,(n_tags(vals)-1) do if ~tag_exist(s, valtags[i]) then s = create_struct(valtags[i], vals.(i), s)

if ~tag_exist(s,'range') then begin
	range = max(abs(obj.val))*[-1,1]
	if (min(obj.val) gt 0) then range[0] = 0
	s = struct_addtags(s, {range:range})
	undefine, range
endif



ysiz = fix(difference(s.yrange)/float(difference(s.xrange))*s.xsiz)
if tag_exist(s, 'ysiz') then s.ysiz = ysiz else s = struct_addtags(s, {ysiz:ysiz})
undefine, ysiz

xmid = round((s.xsiz-1)/2.)
ymid = round((s.ysiz-1)/2.)
holdim = fltarr(s.xsiz, s.ysiz) + s.range[1]
dim = round(difference(minmax(obj.x)/float(s.rnorm)/float(difference(s.xrange))*s.xsiz))
dim = [dim, round(dim * (max(obj.y)-min(obj.y))/float(max(obj.x)-min(obj.x)))]

;obj.x += randomn(seed, n_elements(obj))*0.1
;obj.y += randomn(seed, n_elements(obj))*0.1

triangulate, obj.x, obj.y, tr
;if (n_elements(method) eq 0) then method = 'linear'
grid = griddata(obj.x, obj.y, obj.val, dimension=dim, method=s.method, triangles=tr, min_points=s.min_points)



dim = [n_elements(grid[*,0]), n_elements(grid[0,*])]

if (n_elements(smooth) gt 0) then grid = smooth(grid,s.smooth)

stop
triangulate, obj.x, obj.y, tr, b
xfake = round((obj[b].x - min(obj.x))/float(max(obj.x)-min(obj.x))*(dim[0]-1.))
yfake = round((obj[b].y - min(obj.y))/float(max(obj.y)-min(obj.y))*(dim[1]-1.))

xcoord = rebin(indgen(dim[0]), dim)
ycoord =  transpose(rebin(indgen(dim[1]), reverse(dim)))


out = where(box(xcoord, ycoord, reverse(xfake), reverse(yfake)) eq 1, count)
if ((n_elements(xcoord) - count) lt 10) then begin
	out = where(box(xcoord, ycoord, (xfake), (yfake)) eq 1)
	;plot, xcoord, ycoord, psym=3
	;oplot, xcoord[out], ycoord[out], psym=3, color=fsc_color('red')
endif
if (n_elements(out) lt n_elements(grid)) then grid[out] = s.range[1]		;turn this on to turn off extrapolation



gonts = holdim

xov = round((min(obj.x)/float(s.rnorm) - s.xrange[0])*float(s.xsiz)/float(difference(s.xrange))) 
yov = round((min(obj.y)/float(s.rnorm) - s.yrange[0])*float(s.ysiz)/float(difference(s.yrange)))


if ((xov gt 0) or (yov gt 0)) then begin
	xshift = xov
	yshift = yov
	gonts[xshift,yshift] = grid[0:-xshift-1,0:-yshift-1]
endif else begin
	gonts[0,0] = grid[abs(xov):min([n_elements(grid[*,0])-1,s.xsiz+abs(xov)-1]),abs(yov):min([n_elements(grid[0,*])-1,s.ysiz+abs(yov)-1])]
	;gonts[xov-1:*,*] = range[1]
	;gonts[*,yov-1:*] = range[1]

endelse




;ov = n_elements(grid[*,0]) - n_elements(holdim[*,0])
;if (ov gt 0) then grid = grid[round(ov/2.):-(ov - round(ov/2.)) - 1,*]
;ov = n_elements(grid[0,*]) - n_elements(holdim[0,*])
;if (ov gt 0) then grid = grid[*,round(ov/2.):-(ov - round(ov/2.)) - 1]




;p2 = image(gonts, min_value=range[0], max_value=range[1], yticks=1, xticks=1, xminor=1, xtickname=[' ',' '], yminor=1, ytickname=[' ',' '], xsubticklen=0, xticklen=0.02, rgb_table=rgb_table, axis_style=0, position=pos, current=current)

;;p2 = image(holdim, min_value=range[0], max_value=range[1], yticks=1, xticks=1, xminor=1, xtickname=[' ',' '], yminor=1, ytickname=[' ',' '], xsubticklen=0, xticklen=0.02, rgb_table=rgb_table, axis_style=0, position=pos, current=current)
;p2 = image(reverse(grid,1), min_value=range[0], max_value=range[1], xticks=1, xminor=1, yticks=1, yminor=1, ytickname=[' ',' '], xtickname=[' ',' '], rgb_table=rgb_table, axis_style=0, image_location=[n_elements(holdim[*,0]) - n_elements(grid[*,0]),n_elements(holdim[0,*]) - n_elements(grid[0,*])]/2, /current, /overplot, position=pos)
;;p2 = image(reverse(gonts,1), min_value=range[0], max_value=range[1], xticks=1, xminor=1, yticks=1, yminor=1, ytickname=[' ',' '], xtickname=[' ',' '], rgb_table=rgb_table, axis_style=0, image_location=[0,0], /current, /overplot, position=pos)



;why do i have the reverse here?
;stop
;p2 = image(reverse(gonts,1), min_value=s.range[0], max_value=s.range[1], yticks=1, xticks=1, xminor=1, xtickname=[' ',' '], $
;	yminor=1, ytickname=[' ',' '], xsubticklen=0, xticklen=s.xticklen, rgb_table=s.rgb_table, axis_style=0, position=pos, $
;	current=current, dimensions=dimensions, yticklen=s.yticklen, tickfont_size=s.tickfont_size)

;p2 = image(gonts, min_value=s.range[0], max_value=s.range[1], yticks=1, xticks=1, xminor=1, xtickname=[' ',' '], $
p2 = image(gonts, min_value=s.range[0], max_value=s.range[1], xminor=1, xtickname=[' ',' '], yminor=1, ytickname=[' ',' '], xsubticklen=0, xticklen=s.xticklen, rgb_table=s.rgb_table, axis_style=0, position=pos, $
	current=current, dimensions=s.dimensions, yticklen=s.yticklen, xtickfont_size=s.tickfont_size, ytickfont_size=s.tickfont_size) ;The actual voronoi binned plot


;p2 = image(reverse(grid,1), min_value=range[0], max_value=range[1], xticks=1, xminor=1, yticks=1, yminor=1, ytickname=[' ',' '], xtickname=[' ',' '], rgb_table=rgb_table, axis_style=0, image_location=[n_elements(holdim[*,0]) - n_elements(grid[*,0]),n_elements(holdim[0,*]) - n_elements(grid[0,*])]/2, /current, /overplot, position=pos)
;p2 = image(reverse(gonts,1), min_value=range[0], max_value=range[1], xticks=1, xminor=1, yticks=1, yminor=1, ytickname=[' ',' '], xtickname=[' ',' '], rgb_table=rgb_table, axis_style=0, image_location=[0,0], /current, /overplot, position=pos)

stop




if (n_elements(cont) gt 0) then begin
	
	if ~tag_exist(cont, 'vor_count') then begin
		cont = struct_addtags(cont, arr_struct({vor_count:cont.vorcount, vor_area:cont.area, vor_mu:(cont.vorcount/double(cont.area))}))
		
	endif
	
	ind = where( (cont.vor_Count ne 0) and (cont.vor_Area ne 0) and (cont.y/s.rnorm gt s.yrange[0]) and (cont.y/s.rnorm lt s.yrange[1]) $
		and (cont.x/s.rnorm gt s.xrange[0]) and (cont.x/s.rnorm lt s.xrange[1])) 
	cont = cont[ind]
	
	
	
	
	cdim = round(difference(minmax(cont.x))/float(s.rnorm)/float(difference(s.xrange))*s.xsiz)
	cdim = [cdim, round(cdim * (max(cont.y)-min(cont.y))/float(max(cont.x)-min(cont.x)))]
	
	triangulate, cont.x, cont.y, ctr, b
	cgrid = griddata(cont.x, cont.y, cont.vor_mu, dimension=cdim, /linear, triangles=ctr)
	cdim = [n_elements(cgrid[*,0]), n_elements(cgrid[0,*])]
	
	stop
	
	
	
	conts = holdim
	
	xov = round((min(cont.x)/float(s.rnorm) - s.xrange[0])*float(s.xsiz)/float(difference(s.xrange))) 
	yov = round((min(cont.y)/float(s.rnorm) - s.yrange[0])*float(s.ysiz)/float(difference(s.yrange)))
	
	
	
	
	
	if ((xov gt 0) or (yov gt 0)) then begin
		xshift = xov
		yshift = yov
		conts[xshift,yshift] = cgrid[0:-xshift-1,0:-yshift-1]
;		conts[xshift,yshift] = cgrid
	
	endif else begin
		conts[0,0] = grid[abs(xov):min([n_elements(cgrid[*,0])-1,s.xsiz+abs(xov)-1]),abs(yov):min([n_elements(cgrid[0,*])-1,s.ysiz+abs(yov)-1])]
		;gonts[xov-1:*,*] = range[1]
		;gonts[*,yov-1:*] = range[1]
	
	endelse
	
	;conts = holdim
	;conts[0,0] = cgrid
	;conts = shiftf(conts, round((min(cont.x) - xrange[0])*float(xsiz)/float(difference(xrange))), round((min(cont.y) - yrange[0])*float(ysiz)/float(difference(yrange))))
	
;	p2 = contour(conts, c_value=c_value, /current, /overplot, thick=c_thick, c_thick=c_thick, xmajor=0, xminor=0, ymajor=0, myminor=0)
	
	ind = indgen(n_elements(cont))
	contour, alog10(cont[ind].vor_mu), cont[ind].x, cont[ind].y, /irregular, nlevels=1000, c_colors=fsc_color('white'), path_xy=pathxy, path_info=pathinfo, /path_data_coords
stop


	erase
	cleanplot, /silent
	maxdens = minmax(alog10(cont.vor_mu))
;	psopen, 'density' + roundx(maxdens[0],2)+'_'+roundx(maxdens[1],2), /encapsulated, xs=6, ys=6, /inches, /color
	!p.charsize = 1.3
	!p.charthick=3
	sauron_colormap
	device, decomposed=0
	!p.multi=[0,1,1]
	plot_velfield, cont.x/float(s.rnorm), cont.y/float(s.rnorm), alog10(cont.vor_mu), range=maxdens, xtitle=textoidl('X'), ytitle=textoidl('Y'), /nodots, title='Density: '+roundx(maxdens[0],1)+'/'+roundx(maxdens[1],1)+' ['+textoidl('Log10(# kpc^{-2})]')
	arg = alog10(cont.vor_mu)
	nlevels = 9.
	stop
	contour, arg, cont.x/float(s.rnorm), cont.y/float(s.rnorm), /irregular, c_colors=fsc_color('blk7'), c_thick=0.1, /overplot, levels=(findgen(nlevels)+1)/nlevels*(max(arg)-min(arg)) + min(arg)
	;sharpcorners, thick=3
;	psclose

	stop


	;plot, cont[ind].x, cont[ind].y, psym=3, xrange=[-5,5], yrange=[-5,5], /iso
	nisophotes = n_elements(s.isophote_level)
	str = '((abs(median(dcon.er,/even) - '+strtrim(s.isophote_level[0],2)+'*reff) lt 0.5)'
	if (nisophotes gt 1) then foreach reff_mult, s.isophote_level[1:*] do str = str + ' or (abs(median(dcon.er,/even) - '+strtrim(reff_mult,2)+'*reff) lt 0.5)'
	str = str[0] +')'
	str = '('+str+' and (min(dcon.x)*max(dcon.x) lt 0) and (min(dcon.y)*max(dcon.y) lt 0) )'
	
	for i=0,(n_elements(pathinfo)-1) do begin														;cycle through each contour
		if ((pathinfo[i].type eq 1) and (pathinfo[i].n gt 25)) then begin							;	continue if this is a closed contour with a sufficient number of verticies
			indcon = [indgen(pathinfo[i].n), 0] + pathinfo[i].offset								;		retrieve the indices of this contour's coordinates
			dcon = arr_struct({x:pathxy[0,indcon], y:pathxy[1,indcon]})								;		save to a structure
			q = min([max(abs(pathxy[1,indcon]))/max(abs(pathxy[0,indcon])),1])						;		estimate the axis ratio
			dummy  = ellippar(dcon.x, dcon.y, q, 0d, struc=dcon)									;		generate elliptical coordinates
			dummy = execute('continu = '+str)
			;print, median(dcon.er,/even), continu
			;if (abs(max(dcon.er) - reff) lt 0.5) then begin
			if continu then begin
				dcon = struct_addtags(dcon, replicate({q:q, maxa:max(dcon.a), maxb:max(dcon.b), $		;		add tags to the structure
					maxer:max(dcon.er), meder:median(dcon.er,/even)},n_elements(dcon)))
				con = struct_addtags(con, create_struct('contour'+strtrim(i,2), dcon))					;		save to a master structure
			endif
		endif
	endfor
	;for i=0,(n_tags(con)-1) do oplot, con.(i).x/float(s.rnorm), con.(i).y/float(s.rnorm)

;annuli = fltarr(8,1d4) ;create array to hold xcon_0.5re,ycon_0.5re, xcon_1.5re,ycon_1.5re, etc.
	
	if 0 then begin																					;plot all the contours
		for i=0,(n_tags(con)-1) do begin
			xcon = (((con.(i).x)/float(s.rnorm) - 0*s.xrange[0])*float(s.xsiz)/float(difference(s.xrange))) + s.xsiz/2.
			ycon = (((con.(i).y)/float(s.rnorm) - 0*s.yrange[0])*float(s.ysiz)/float(difference(s.yrange))) + s.ysiz/2.
			p2 = plot(smooth(xcon,s.smcon), smooth(ycon,s.smcon), position=pos, /current, /overplot, thick=s.c_thick)
		endfor
	endif else begin																				;plot just the closest contour to Reff
		er_arr = fltarr(n_tags(con))
		for i=0,(n_tags(con)-1) do er_arr[i] = con.(i)[0].meder
;		k=0 ;start counter for recording kth isophote contour
		foreach mult_fac, s.isophote_level do begin
			mindist = min(abs(er_arr - mult_fac[0]*reff),imin)
			if (mindist lt 0.1) then begin
				;oplot, con.(imin).x, con.(imin).y, color=fsc_color('red')
				;xcon = smooth((((con.(imin).x)/float(s.rnorm) - s.xrange[0]*0)*float(s.xsiz)/float(difference(s.xrange))) + s.xsiz/2., s.smcon, /edge_truncate)
				;ycon = smooth((((con.(imin).y)/float(s.rnorm) - s.yrange[0]*0)*float(s.ysiz)/float(difference(s.yrange))) + s.ysiz/2., s.smcon, /edge_truncate)
				xcon = (((con.(imin).x)/float(s.rnorm) - s.xrange[0]*0)*float(s.xsiz)/float(difference(s.xrange))) + s.xsiz/2
        xcon = smooth([xcon,xcon[0]], s.smcon, /edge_wrap)
        ycon = (((con.(imin).y)/float(s.rnorm) - s.yrange[0]*0)*float(s.ysiz)/float(difference(s.yrange))) + s.ysiz/2.
        ycon = smooth([ycon,ycon[0]], s.smcon, /edge_wrap)
        p2 = plot(xcon, ycon, position=pos, /current, /overplot, thick=s.c_thick, axis_style=0)
				;xcon,ycon = contours for each isophote (1 & 4 R_e)
				;be sure to require "buffer isophotes" in xmap_justin, i.e. 0.5,1.5,3.5,4.5 R_e.
;        annuli[k,*] = xcon ;feed xcon to annulus coordinates array
;        annuli[k+1,*] = ycon ;" ycon ""

				;d = con.(imin)
				;if (mult_fac eq s.isophote_level[-1]) then stop			
				;p2 = plot([xcon, xcon[0]], [ycon, ycon[0]], position=pos, /current, /overplot, thick=s.c_thick, axis_style=0)
			
			endif
		endforeach
	endelse
	
;  cloudyplot, obj, annuli ; initiate cloudplot procedure
	
	
	
	ticknames = s.xticknames
	tickvalues = ticknames*float(s.xsiz)/float(difference(s.xrange)) + xmid
	;topaxis = axis('X', location=[0,s.ysiz], tickdir=1, position=pos, target=p2, textpos=0, ticklen=s.xticklen, subticklen=0.5, $
topaxis = axis('X', location=[0,s.ysiz], tickdir=1, target=p2, textpos=0, ticklen=s.xticklen, subticklen=0.5, $
		tickfont_size=s.tickfont_size, tickvalues=tickvalues, tickname=strarr(n_elements(tickvalues))+' ', minor=s.xminor)
;	if (n_elements(nox) eq 0) then topaxis = axis('X', location=[0,0], tickdir=0, position=pos, target=p2, textpos=0, ticklen=s.xticklen, $
  if (n_elements(nox) eq 0) then topaxis = axis('X', location=[0,0], tickdir=0, target=p2, textpos=0, ticklen=s.xticklen, $

									subticklen=0.5, tickfont_size=s.tickfont_size, tickvalues=tickvalues, tickname=ticknames, minor=s.xminor, title=s.xtitle) $
							else topaxis = axis('X', location=[0,0], tickdir=0,  target=p2, textpos=0, ticklen=s.xticklen, subticklen=0.5, $
								tickfont_size=s.tickfont_size, tickvalues=tickvalues, tickname=strarr(n_elements(tickvalues))+' ', minor=s.xminor)
		
	
	
	ticknames = s.yticknames
	tickvalues = ticknames*float(s.ysiz)/float(difference(s.yrange)) + ymid
	if (n_elements(noy) eq 0) then yaxis = axis('Y', location=[0,0], tickdir=0,  target=p2, ticklen=s.yticklen, subticklen=0.5, $
										tickvalues=tickvalues, tickfont_size=s.tickfont_size, tickname=ticknames, major=1, title=s.ytitle, minor=s.yminor) $
								else yaxis = axis('Y', location=[0,0], tickdir=0,  target=p2, ticklen=s.yticklen, subticklen=0.5, $
									tickvalues=tickvalues, tickfont_size=3, tickname=strarr(n_elements(tickvalues))+' ', major=1, minor=s.yminor)
	yaxis = axis('Y', location=[s.xsiz,0], tickdir=1,  target=p2, ticklen=s.yticklen, subticklen=0.5, tickfont_size=s.tickfont_size, $
			tickvalues=tickvalues, tickname=strarr(n_elements(tickvalues))+' ', major=1, minor=s.yminor)
	




	
endif

;sharpcorner, p2, thick=1, position=pos, dimension=s.dimensions
;print, con.(imin)[0].maxer



end
