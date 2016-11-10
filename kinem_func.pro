;+
;PURPOSE
;	Calcualte a velocity given position angles and kinemetric coefficients
;SYNTAX
;	kinem_func(position angle, [rot. amplitude, kinematic PA, kin. ellipticity, systemic velocity])
;Written by Jacob A. Arnold, 8-22-09, UCSC
;-

function kinem_func, x, coef

PA = (x - coef[1])

if (n_elements(coef) eq 5) then PA = PA*2d


result = coef[0] / sqrt(1d + (tan(PA)/coef[2])^(2d))								;ellipsoidal figure
result = result * sign(cos(PA))						;unless specified, this is an odd function
result = result + coef[3]															;add on the systemic value

;!p.multi=[0,1,2]
;plot, x, result, psym=3
;plot, x, sign(cos(pa)), yrange=[-2,2]

;if (coef[0] eq -50) then stop


return, result





;return, coef[0] / sqrt(1d + (tan(PA)/coef[2])^(2d)) * cos(PA)/abs(cos(PA)) + coef[3]
;
;return, coef[0] / sqrt(1d + (tan(X - coef[1])/coef[2])^(2d)) * cos(X - coef[1])/abs(cos(X - coef[1])) + coef[3]
			
end
