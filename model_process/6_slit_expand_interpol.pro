restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_cs.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_h.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_n2.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_wc.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_va.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_vaxx.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_vayy.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_3D/d3_vaz.sav'

;step0:角度与bounce点选取
theta=0*!pi/6.0;单位：弧度？
x0=0*1d6;单位：m

;step1:切片兼扩展
cs_1=dblarr(1531,512)
h_1=cs_1
n2_1=cs_1
wc_1=cs_1
va_1=cs_1
var_1=cs_1;从今以后二维vax改称var
vaz_1=cs_1




for i = 0, 765 do begin;3d模型(jx,jy)的中心像素是255，此处2d(i)的是765，此处的水平分辨率为96000m，注意此处并没有截取到3d模型的边界

    jx=255+round(x0+i*96000.0*cos(theta))/96000.0
    jy=255+round(i*96000.0*sin(theta))/96000.0

    if (jx-255.0)^2 + (jy-255.0)^2 LT 254.0^2 then begin
        cs_1[765+i,*]=d3_cs[jx,jy,*]
        h_1[765+i,*]=d3_h[jx,jy,*]
        n2_1[765+i,*]=d3_n2[jx,jy,*]
        wc_1[765+i,*]=d3_wc[jx,jy,*]
        va_1[765+i,*]=d3_va[jx,jy,*]
        vaz_1[765+i,*]=d3_vaz[jx,jy,*]

        var_1[765+i,*]=d3_vayy[jx,jy,*]*sin(theta)+d3_vaxx[jx,jy,*]*cos(theta)
        ;print,765+i,cs_1[765+i,480],jx,jy
    endif else begin
        cs_1[765+i,*]=cs_1[765+i-1,*]
        h_1[765+i,*]=h_1[765+i-1,*]
        n2_1[765+i,*]=n2_1[765+i-1,*]
        wc_1[765+i,*]=wc_1[765+i-1,*]
        va_1[765+i,*]=va_1[765+i-1,*]
        vaz_1[765+i,*]=vaz_1[765+i-1,*]
        var_1[765+i,*]=var_1[765+i-1,*]
        ;print,'     ',765+i-1
    endelse
endfor


for i = 0, 765 do begin;3d模型(jx,jy)的中心像素是255，此处2d(i)的是765，此处的水平分辨率为96000m

    jx=255+round(x0-i*96000.0*cos(theta))/96000.0
    jy=255+round(-1*i*96000.0*sin(theta))/96000.0

    if (jx-255.0)^2 + (jy-255.0)^2 LT 254.0^2 then begin
        cs_1[765-i,*]=d3_cs[jx,jy,*]
        h_1[765-i,*]=d3_h[jx,jy,*]
        n2_1[765-i,*]=d3_n2[jx,jy,*]
        wc_1[765-i,*]=d3_wc[jx,jy,*]
        va_1[765-i,*]=d3_va[jx,jy,*]
        vaz_1[765-i,*]=d3_vaz[jx,jy,*]

        var_1[765-i,*]=d3_vayy[jx,jy,*]*sin(theta)+d3_vaxx[jx,jy,*]*cos(theta)
    endif else begin
        cs_1[765-i,*]=cs_1[765-i+1,*]
        h_1[765-i,*]=h_1[765-i+1,*]
        n2_1[765-i,*]=n2_1[765-i+1,*]
        wc_1[765-i,*]=wc_1[765-i+1,*]
        va_1[765-i,*]=va_1[765-i+1,*]
        vaz_1[765-i,*]=vaz_1[765-i+1,*]
        var_1[765-i,*]=var_1[765-i+1,*]
    endelse
endfor

;step2:interpol
cs=dblarr(1531,2556)
h=cs
n2=cs
va=cs
var=cs
vaz=cs
wc=cs

for i = 0, 1530 do begin
    cs[i,*]=interpol(reform(cs_1[i,*]),2556,/spline)
    h[i,*]=interpol(reform(h_1[i,*]),2556,/spline)
    n2[i,*]=interpol(reform(n2_1[i,*]),2556,/spline)
    ;va[i,*]=interpol(reform(va_1[i,*]),2556,/spline)
    var[i,*]=interpol(reform(var_1[i,*]),2556,/spline)
    vaz[i,*]=interpol(reform(vaz_1[i,*]),2556,/spline)
    wc[i,*]=interpol(reform(wc_1[i,*]),2556,/spline)
endfor

cs0=cs
h0=h
n20=n2
va0=va
var0=var
vaz0=vaz
wc0=wc

cs=dblarr(7651,2556)
h=cs
n2=cs
va=cs
var=cs
vaz=cs
wc=cs

for j = 0, 2555 do begin
    cs[*,j]=interpol(reform(cs0[*,j]),7651,/spline)
    h[*,j]=interpol(reform(h0[*,j]),7651,/spline)
    n2[*,j]=interpol(reform(n20[*,j]),7651,/spline)   
    ;va[*,j]=interpol(reform(va0[*,j]),7651,/spline)  
    var[*,j]=interpol(reform(var0[*,j]),7651,/spline) 
    vaz[*,j]=interpol(reform(vaz0[*,j]),7651,/spline)  
    wc[*,j]=interpol(reform(wc0[*,j]),7651,/spline)  
endfor
va=sqrt(var^2+vaz^2)

save,cs,filename='cs.sav'
save,h,filename='h.sav'
save,n2,filename='n2.sav'
save,va,filename='va.sav'
save,var,filename='var.sav'
save,vaz,filename='vaz.sav'
save,wc,filename='wc.sav'
end