restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/cs.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/h.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/n2.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/va.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/vax.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/vaz.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/wc.sav'

cs0=cs
h0=h
n20=n2
va0=va
vax0=vax
vaz0=vaz
wc0=wc

cs=dblarr(1531,2556)
h=cs
n2=cs
va=cs
vax=cs
vaz=cs
wc=cs


for i = 0, 1530 do begin
    cs[i,*]=interpol(reform(cs0[i,*]),2556,/spline)
    h[i,*]=interpol(reform(h0[i,*]),2556,/spline)
    n2[i,*]=interpol(reform(n20[i,*]),2556,/spline)
    ;va[i,*]=interpol(reform(va0[i,*]),2556,/spline)
    vax[i,*]=interpol(reform(vax0[i,*]),2556,/spline)
    vaz[i,*]=interpol(reform(vaz0[i,*]),2556,/spline)
    wc[i,*]=interpol(reform(wc0[i,*]),2556,/spline)
endfor

cs0=cs
h0=h
n20=n2
va0=va
vax0=vax
vaz0=vaz
wc0=wc

cs=dblarr(7651,2556)
h=cs
n2=cs
va=cs
vax=cs
vaz=cs
wc=cs

for j = 0, 2555 do begin
    cs[*,j]=interpol(reform(cs0[*,j]),7651,/spline)
    h[*,j]=interpol(reform(h0[*,j]),7651,/spline)
    n2[*,j]=interpol(reform(n20[*,j]),7651,/spline)   
    ;va[*,j]=interpol(reform(va0[*,j]),7651,/spline)  
    vax[*,j]=interpol(reform(vax0[*,j]),7651,/spline) 
    vaz[*,j]=interpol(reform(vaz0[*,j]),7651,/spline)  
    wc[*,j]=interpol(reform(wc0[*,j]),7651,/spline)  
endfor
va=sqrt(vax^2+vaz^2)

end