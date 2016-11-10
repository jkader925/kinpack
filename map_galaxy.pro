pro map_galaxy, galname, s, p2=p2, nox=nox, noy=noy, point=point, current=current, $
	pos=pos, outline=outline, dataloc=dataloc, isophote=isophote, nodiscrete=nodiscrete, $
	nostellar=nostellar, nopne=nopne, nogcs=nogcs, nosauron=nosauron, nolong=nolong, $
	nosmeag=nosmeag, doreff=doreff, norm=norm, sneb=sneb, dneb=dneb, fill=fill, $
	veldisp=veldisp, h3=h3, h4=h4, dss=dss, invert=invert, log=log, varplot=varplot, $
	new=new, maxPoints=maxPoints, minerr=minerr, red=red, blue=blue, green=green, z=z, $
	noextrap=noextrap, holdim=holdim, gsm=gsm, bb=bb, maxDist=maxdist, vor=vor, $ 
	rms=rms, vmax=vmax, medianing=medianing,vs_fitting=vs_fitting

	map = obj_new('map_'+galname, s, doreff=doreff, norm=norm, gsm=gsm, rms=rms, veldisp=veldisp, h3=h3, h4=h4, maxDist=maxDist, vor=vor, dss=dss)			;instantiate the object
	map->get_gal_info,vs_fitting						;get galaxy specific information
	
;		if 0*(n_elements(dss) gt 0) then map->get_dss, vs_fitting, invert=invert, log=log $
;	else 
	map->get_data, nopne=nopne, vs_fitting=vs_fitting, nogcs=nogcs, $						;get the galaxy specific data
		nodiscrete=nodiscrete, nostellar=nostellar, $
		nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, $
		red=red, blue=blue;, green=green

	if (n_elements(point) gt 0) then map->point_symmetrize				;point-symmetrize if the point keyword is set
	if 1+0*(n_elements(dss) eq 0) then map->get_variogram, varplot=varplot $;, maxPoints=maxPoints 
		else map->map_obj, sneb=sneb, dneb=dneb					;data smoothing
;			map->map_obj, sneb=sneb, dneb=dneb         ;data smoothing

;	if (n_elements(dss) gt 0) then map->map_obj, sneb=sneb, dneb=dneb;, invert=invert
	map->plotting_init, s, pos=pos, current=current						;initialize the plotting variables

	map->make_grid, new=new, z=z, varz=varz, minerr=minerr, noextrap=noextrap, holdim=holdim, fill=fill				;triangulate onto a grid
	
	map->change_units													;convert to the plotting units

	map->plot, p2=p2, nox=nox, noy=noy, s=s, fill=fill, galname						;plot the grid

;	if (n_elements(outline) gt 0) then $								;outline the data if the outline keyword is set
;		map->outline_data, p2=p2, bb=bb

    if (n_elements(vmax) gt 0) then $                   ;extract data in rectangle along major axis, plot as function of position and get vmax -jk
    map->extract_vmax, nopne=nopne, nogcs=nogcs, $            
    nodiscrete=nodiscrete, nostellar=nostellar, $
    nosauron=nosauron, nolong=nolong, nosmeag=nosmeag, $
    red=red, blue=blue
  

;	if (n_elements(dataloc) gt 0) then $								;plot the data locations if the dataloc keyword is set
;		map->plot_data_loc, p2=p2
;
;
	if (n_elements(isophote) gt 0) then $								;plot the isophotes if the isophote keyword is set
		map->isophote

;	map->sharpcorner, p2=p2
p2.refresh ;seems to initiate "/overplotted" or possibly "/current" plot elements to graphics window
;	obj_destroy, map


end