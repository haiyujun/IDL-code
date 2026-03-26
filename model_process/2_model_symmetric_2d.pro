restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_old/cs.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_old/H_1.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_old/N2_1.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_old/va.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_old/vax.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_old/vaz.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_old/ωc_1.sav'

cs0=cs
h0=h
n20=n2
va0=va
vax0=vax
vaz0=vaz
wc0=wc

cs=dblarr(1531,512)
h=cs
n2=cs
va=cs
vax=cs
vaz=cs
wc=cs

cs[510:765,*]=reverse(cs0)
h[510:765,*]=reverse(h0)
n2[510:765,*]=reverse(n20)
va[510:765,*]=reverse(va0)
vax[510:765,*]=-1*reverse(vax0);va是矢量，vax中轴两边反号！
vaz[510:765,*]=reverse(vaz0)
wc[510:765,*]=reverse(wc0)

cs[765:1020,*]=cs0
h[765:1020,*]=h0
n2[765:1020,*]=n20
va[765:1020,*]=va0
vax[765:1020,*]=vax0
vaz[765:1020,*]=vaz0
wc[765:1020,*]=wc0

for i = 0, 509 do begin
    cs[i,*]=cs0[255,*]
    h[i,*]=h0[255,*]
    n2[i,*]=n20[255,*]
    va[i,*]=va0[255,*]
    vax[i,*]=-1*vax0[255,*];va是矢量，vax经过中轴后反号！
    vaz[i,*]=vaz0[255,*]
    wc[i,*]=wc0[255,*]
endfor

for i = 1021, 1530 do begin
    cs[i,*]=cs0[255,*]
    h[i,*]=h0[255,*]
    n2[i,*]=n20[255,*]
    va[i,*]=va0[255,*]
    vax[i,*]=vax0[255,*]
    vaz[i,*]=vaz0[255,*]
    wc[i,*]=wc0[255,*]
endfor

;ct = reverse(COLORTABLE(72))
;im=image(vaz^2+vax^2-va^2,RGB_TABLE=ct,position=[0.1,0.2,0.9,1]);,max_value=6,min_value=2)
;ip=plot(reform(r),reform(z),/nodata,/current,xrange=[min(r),max(r)],yrange=[min(z),max(z)],xstyle=1,ystyle=1,position=[0.3,0.2,0.7,1],xtitle='r/Mm',ytitle='z/Mm')
;c=colorbar(target=im)
end