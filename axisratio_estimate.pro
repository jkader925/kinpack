PRO Axisratio_estimate, data, reff, sige, q, diagnostics=diagnostics, cosmo

; Isolate 1 R_e contour has two methods:
;   Method 1: Calculate semimajor axis for each contour, select contour with
;     smallest offset between its semimajor axis length and the x = +1 R_e.
;
;   Method 2: Calculate semimajor axis length for each contour. Isolate contours
;     with semimajor axes that differ from 1 R_e in length by the amount 'del' in
;     Kpc. The 1 R_e contour is the one with the largest number of points in this
;     subset, (the most well-behaved presumably.)
;
;   
; Calculating the axis ratio:
;     Take points from within +/- xbox and +/- ybox from the 1 R_e contour. Store the
;     absolute value of the x/y-offsets from the center. The median
;     x- or y- distances of each of these (two) subsets of points from the center is
;     a robust estimate of the semi-major and semi-minor axis lengths of the contour.   
;     Therefore, q = median(semi-minor axis length / semi-major axis length)

;---------------------------------
; axis ratio estimator program variables:

  ;turn on/off method 1
  meth1 = 0
  ;turn on/off method 2
  meth2 = 1
  
  
  ;For method 2 of finding the 1 R_e contour, use only contours with
  ;the following offset [Kpc] between their semimajor axis length and x=+1 R_e:
  del = 0.5

  ;When calculating the axis ratio of the 1 R_e contour, only points within the following
  ;ranges will be used:
  xbox = 1.0
  ybox = 0.5
;---------------------------------


N=15d  ; number of bins in x,y directions for 2d density array
;if n_elements(cosmo) eq 1 then N=5d
siz=N/2.  ; data used will be -N/2 < data.x,data.y < N/2
range = [-1,1]*siz
rnorm = 1.
xsiz=siz
xrange=[-1,1]*siz
yrange=[-1,1]*siz
dat=data[where(data.x GT xrange[0] and data.x LT xrange[1] and data.y GT yrange[0] and data.y LT yrange[1])]
ysiz = fix(difference(yrange)/float(difference(xrange))*xsiz)

xmid = round((xsiz-1)/2.)
ymid = round((ysiz-1)/2.)
holdim = fltarr(xsiz, ysiz) + range[1]
dim = round(difference(minmax(dat.x)/float(rnorm)/float(difference(xrange))*xsiz))
dim = [dim, round(dim * (max(dat.y)-min(dat.y))/float(max(dat.x)-min(dat.x)))]


dat.x = dat.x + siz
dat.y = dat.y + siz


smoothed = fltarr(N,N)  ; this will be the smoothed array. has dimensions = dimensions of griddata grid

scale = max(dat.x)/N  ;scaling between x,y data range and smoothing array dimensions





smoothed = replicate({x:0d, y:0d, rho:0d}, N^2d)

for i = 0,n_elements(smoothed)-1 do smoothed[i].x = i mod N
for i = 0,n_elements(smoothed)-1 do smoothed[i].y = (i - smoothed[i].x)/N

i=0
j=0
;k=0
l=0 
D = fltarr(n_elements(smoothed),(n_elements(dat.x)))
 a1 = systime()
print,'starting the Gaussian smoothing...'
for i = 0, n_elements(smoothed)-1 do begin  ;for each grid bin center
;    for j = 0, n_elements(dat.x)-1 do begin  ;cycle through the data points
;        D[i,j] = (smoothed[i].x-dat[j].x)^2d + (smoothed[i].y-dat[j].y)^2d ;distance-squared from current (i,j)th element to the (k,l)th column
        D[i,*] = (smoothed[i].x-dat.x)^2d + (smoothed[i].y-dat.y)^2d
;    endfor
endfor
a2 = systime()
i=0
for i=0,n_elements(smoothed)-1 do begin
   smoothed[i].rho = total(smoothed[i].rho + dat.mass*exp(-0.5d*D[i,*]/((0.5d)*(0.05*1.5)^2d))) ; D has a D[i] (a vector) for each smoothed[i] grid point!
endfor
print,'smoothing complete!'




      ;------calc time elapsed
      c1=strmid(a1,11,8)
      c1=strjoin(strsplit(c1,':',/extract),' ')
      h1 = strmid(c1,0,2)    ;hours
      m1 = strmid(c1,3,2)    ;minutes
      s1 = strmid(c1,6,2)    ;seconds
      c2=strmid(a2,11,8)
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

      print,strcompress(string(dh)+' hour(s)'+string(dm)+' minute(s)'+string(ds)+' second(s)'+' elapsed for smoothing algorithm')
      ;------

;end






inds=fltarr(N,N)
for i = 0,N-1 do begin
  for j = 0,N-1 do begin
    inds[i,j] = [smoothed[where(smoothed.x eq i and smoothed.y eq j)].rho]
  endfor
endfor
    
datx = data.x + N/2.
daty = data.y + N/2.

inds=alog(inds)

if diagnostics GT 0 then p=plot(datx,daty,linestyle=6,symbol='dot',sym_size=0.5,xrange=[0,N],yrange=[0,N],xtitle='Y [kpc]', ytitle='Z [kpc]',title='Hoffman 1:3, orbit m32')
if diagnostics GT 0 then c=contour(inds,n_levels=6,rgb_table=33,transparency=45,overplot=p,/fill)
isocontour,inds,outverts,outconn,n_levels=60,level_values=vals,/double  ; isocontour:outverts has (x,y)-coords of contours. These have minmax = dimensions of SMOOTHED ARRAY!
unique = uniq(outverts[2,*])

;'isocontour' returns a list called 'outverts', which has x-coordinate (col.1), y-coordinate (col.2), and a float-type ID (col.3) of each contour. For each contour the x,y-coords
;are ordered sequentially starting from quadrant III and proceeding in CCW sense.

;---------------------------------------------------------------------------------------------------------------------------------------------
struc = replicate({val:0d,num:0d},1) ;a new list to keep track of the number of (x,y)-pairs belonging to each contour.
c=0

for h=0,n_elements(outverts[2,*])-1 do begin
  c=c+1
  if h ne 0 and outverts[2,h] ne outverts[2,h-1] then begin ; IF (x,y)-pairs belonging to the same contour ID, (and not the 0th (x,y)-pair)...
    num=c-1
    val=float(outverts[2,h-1])
    append=replicate({val:0d,num:0d},1)
    append.val = val ;ID
    append.num = num ; number of (x,y)-pairs belonging to each ID.
    struc = struct_append(struc,append)
    c=0
  endif
endfor
;---------------------------------------------------------------------------------------------------------------------------------------------
;Find the effective radius contour. I'm just using the contour with semi-major axis length closest to 1 R_e.

  majrad = fltarr(2,n_elements(struc)) ; column 1: max x-value for each unique contour, column 2: contour ID, List length: number of unique contours
  majrad[1,*] = struc.val ;fill column 2 with IDs
  for i = 1,n_elements(struc)-1 do begin
    majrad[0,i] = abs(max(abs(outverts[0,where(float(outverts[2,*]) eq struc[i].val)])) - (reff+N/2)) 
    endfor
    ;Go through all contours. For each contour, take only the absolute value of x-coordinates.
    ;Find the point of the contour that has largest abs(x). Calculate the x-offset of this point from x = (center)+(1*R_e).


  majrad = majrad[*,1:n_elements(majrad[0,*])-1] ;throw out first contour
  majrad = majrad[*,sort(majrad[1,*])] ; sort by float-type ID

  minind = majrad[1,where(majrad[0,*] eq min(majrad[0,*]))] ; obtain ID of contour with smallest difference between its semimajor axis and the value 1 R_e
  minind = minind[0] ; I hate IDL sometimes... (converted 1-element array into scalar)

;alternate method of finding 1 R_e contour... 

IF meth2 EQ 1 THEN BEGIN
  ;using only contours with semi-major axes 1 R_e +/- 0.9 Kpc 
  mincontvals = majrad[1,where(majrad[0,*] LE del)] ; IDs of contours which have offsets of less than 0.9 Kpc between semimajor axis and 1 R_e.
  c_vals = float(outverts[2,*]) ; the float-type IDs
  c_vals = c_vals[*]; type conversion from 1xM array to 1-d list
  holdval = fltarr(2,n_elements(mincontvals))
  for u = 0,n_elements(mincontvals)-1 do begin
    holdval[*,u] = [mincontvals[u],n_elements(c_vals[where(c_vals eq mincontvals[u])])] ; each row in holdval contains 'mincontval' contour IDs (col.1) and that contour's # of constituent points (col.2) 
  endfor
  holdval=holdval[*,sort(holdval[1,*])]
  minind2 = holdval[0,where(holdval[1,*] eq max(holdval[1,*]))] ;ID of close-to-1R_e contour, with the most amount of points
  minind2 = minind2[0] ; in case there are multiple output for the same contour. Why is this happening? ---> Maybe some contours have same # points
ENDIF

;ind = where(unique eq minind) ;float-type ID of 1 R_e contour from method 1
;ind2 = where(unique eq minind2) ;float-type ID of 1 R_e contour from method 2

;---------------------------------------------------------------------------------------------------------------------------------------------
;Calculate the axis ratio -- proceed using 1 R_e contour from method 1:
IF meth1 EQ 1 THEN BEGIN
reiso = outverts[*,where(float(outverts[2,*]) eq minind)] ;extract (x,y)-pairs for appropriate contour, by cross-matching IDs
reisoplot2 = reiso
ENDIF

IF meth2 EQ 1 THEN BEGIN
;Calculate the axis ratio -- proceed using 1 R_e contour from method 2:
reiso = outverts[*,where(float(outverts[2,*]) eq minind2)] ;extract (x,y)-pairs for appropriate contour, by cross-matching IDs
reisoplot2 = reiso
ENDIF


;axis len.  = center - (x/y-coord of contour points between y/x=(center)+[-1,1] Kpc)
minorradius = (N/2d) - (reisoplot2[1,where(reisoplot2[0,*] GE (N/2d)-xbox and reisoplot2[0,*] LE (N/2d)+xbox)])
majorradius = (N/2d) - (reisoplot2[0,where(reisoplot2[1,*] GE (N/2d)-ybox and reisoplot2[0,*] LE (N/2d)+ybox)])
minorradius = abs(minorradius[uniq(minorradius)]) ;collect positive, unique, distances for points
majorradius = abs(majorradius[uniq(majorradius)])
minorradius = minorradius[where(minorradius GT 0)] ;collect only positive, non-zero y/x-offsets
majorradius = majorradius[where(majorradius GT 1)]  ;why was this GT 1??

q = max(minorradius)/max(majorradius) ; choose min, since the estimator tends to overpredict the axis ratio!

if q GE 1.0 then q = 0.99
;---------------------------------------------------------------------------------------------------------------------------------------------
;Plotting


if diagnostics GT 0 then p=plot(reisoplot2[0,*],reisoplot2[1,*],color='green',thick=2,/overplot)
if diagnostics GT 0 then t=text(1,27,strcompress('q ='+string(q)),color='red',target=p,/DATA)
print,strcompress('photometric q ='+string(q))
;----------------------------------------------------------------------------------------------------------------------------------------------
;plot isophote ellipse

xsiz = N + 1;357.;500;
ysiz = N + 1;285.7
;s.xrange = [-7,7]
yrange = [0,N]
xrange = [0,N]
;-----------------------------------------------------------------------------------------------
nbound = 100
xmid = ((xsiz-1)/2.)                       ;x-midpoint of the map
ymid = ((ysiz-1)/2.)                        ;y-midpoint of the map






maxa = 1*reff  ; DONT KNOW BUT PRETTY SURE THIS IS ACTUALLY 2 R_e; shorten slightly to not stretch map
xpa = 0

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



;  p2 = plot(xhigh1, yhigh1, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0,color='lime green',thick=2)
  
;  p2 = plot(xhigh2, yhigh2, position=pos, /current, /overplot, xmajor=0, xminor=0, ymajor=0)
END