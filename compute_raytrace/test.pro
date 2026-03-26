z=mrdfits('/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/z.fits');单位：Mm;z=0.032*j-15.6088
;！！！！！！全部转化为国际单位制！！！！！！
;x=x*1d6
z=z*1d6
x=dindgen(1531,increment=96000,start=-255*9.6*1d4*3) 
;以r=0为轴将模型扩展后的1531✖️512轴对称模型
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/cs.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/h.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/n2.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/va.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/vax.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/vaz.sav'
restore,'/Users/lihaiyu/Desktop/IDL code/forward_modeling/param_SI/param_2Dexpand/wc.sav'

;计算各参数随x，z变化的导数
dcs_dx=dblarr(1531,512)
dh_dx=dblarr(1531,512)
dn2_dx=dblarr(1531,512)
dva_dx=dblarr(1531,512)
dvax_dx=dblarr(1531,512)
dvaz_dx=dblarr(1531,512)
dwc_dx=dblarr(1531,512)

dcs_dz=dblarr(1531,512)
dh_dz=dblarr(1531,512)
dn2_dz=dblarr(1531,512)
dva_dz=dblarr(1531,512)
dvax_dz=dblarr(1531,512)
dvaz_dz=dblarr(1531,512)
dwc_dz=dblarr(1531,512)

for i = 0, 1530 do begin
    dcs_dz[i,*]=deriv(reform(z),cs[i,*])
    dh_dz[i,*]=deriv(reform(z),h[i,*])
    dn2_dz[i,*]=deriv(reform(z),n2[i,*])
    dva_dz[i,*]=deriv(reform(z),va[i,*])
    dvax_dz[i,*]=deriv(reform(z),vax[i,*])
    dvaz_dz[i,*]=deriv(reform(z),vaz[i,*])
    dwc_dz[i,*]=deriv(reform(z),wc[i,*])
endfor

for i = 0, 511 do begin
    dcs_dx[*,i]=deriv(reform(x),cs[*,i])
    dh_dx[*,i]=deriv(reform(x),h[*,i])
    dn2_dx[*,i]=deriv(reform(x),n2[*,i])
    dva_dx[*,i]=deriv(reform(x),va[*,i])
    dvax_dx[*,i]=deriv(reform(x),vax[*,i])
    dvaz_dx[*,i]=deriv(reform(x),vaz[*,i])
    dwc_dx[*,i]=deriv(reform(x),wc[*,i])
endfor


;step0：模拟参数设置
omg=0.001*4.5*2*!pi;ω
;F=omg^4-omg^2*(cs^2+va^2)*(p^2+q^2)+cs^2*(p^2+q^2)*(vax*p+vaz*q)^2-wc^2*(omg^2-cs^2*q^2)+cs^2*N2*p^2=0
step=15000;总步数
l=5000.0;步长
y=dblarr(step+1,4);第一维是步数，第二维分别代表——0:x;1:z;2:p or kx;3:q or kz。
S=dblarr(step+1);相位=∫kxdx+kzdz，除以omg即为时间
;step1：初始值设置

y[0,0]=-40*1d6;x0(m)
y[0,1]=-7.2819*1d6;z0(m)
;y[0,2]:p0(m^-1)，用初始时刻的色散关系求
y[0,3]=0*1d-9;q0=0

jx=round((y[0,0]+73440000.0)/96000.0)
jz=round((y[0,1]+15608800.0)/32000.0)
cs_o=cs[jx,jz]
H_o=H[jx,jz]
N2_o=N2[jx,jz]
va_o=va[jx,jz]
vax_o=vax[jx,jz]
vaz_o=vaz[jx,jz]
wc_o=wc[jx,jz]

coeffs=[omg^4 - omg^2*(cs_o^2+va_o^2)*y[0,3]^2 + cs_o^2*vaz_o^2*y[0,3]^4 - wc_o^2*(omg^2-cs_o^2*y[0,3]^2), $;0阶系数
2*cs_o^2*vax_o*vaz_o*y[0,3]^3, $;1阶系数
-omg^2*(cs_o^2+va_o^2) + cs_o^2*va_o^2*y[0,3]^2 + cs_o^2*N2_o, $;2阶系数
2*cs_o^2*vax_o*vaz_o*y[0,3], $;3阶系数
cs_o^2*vax_o^2];4阶系数

root=FZ_ROOTS(coeffs)
y[0,2]=root[1];!!!!!!!!!!!!!!!!!!!!!!!
;y[0,2]=-1.3232367*1d-6

;step2：计算当前x,z坐标所在的格子,提取其各项参数，然后rk4
k1=dblarr(step+1,4)
k2=dblarr(step+1,4)
k3=dblarr(step+1,4)
k4=dblarr(step+1,4);rk4参数

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
disp=dblarr(step+1);色散关系D=0
jxx=dblarr(step+1)
jzz=dblarr(step+1)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for i = 0, step-1 do begin
    
    jx=round((y[i,0]+73440000.0)/96000.0)
    jz=round((y[i,1]+15608800.0)/32000.0)

    jxx[i]=jx
    jzz[i]=jz

    cs_o=cs[jx,jz]
    H_o=H[jx,jz]
    N2_o=N2[jx,jz]
    va_o=va[jx,jz]
    vax_o=vax[jx,jz]
    vaz_o=vaz[jx,jz]
    wc_o=wc[jx,jz]

    dcs_dx_o=dcs_dx[jx,jz]
    dh_dx_o=dh_dx[jx,jz]
    dN2_dx_o=dn2_dx[jx,jz]
    dva_dx_o=dva_dx[jx,jz]
    dvax_dx_o=dvax_dx[jx,jz]
    dvaz_dx_o=dvaz_dx[jx,jz]
    dwc_dx_o=dwc_dx[jx,jz]

    dcs_dz_o=dcs_dz[jx,jz]
    dh_dz_o=dh_dz[jx,jz]
    dN2_dz_o=dn2_dz[jx,jz]
    dva_dz_o=dva_dz[jx,jz]
    dvax_dz_o=dvax_dz[jx,jz]
    dvaz_dz_o=dvaz_dz[jx,jz]
    dwc_dz_o=dwc_dz[jx,jz]

    ;rk4
    ;dydl=[$
;    -2*(omg^2)*(cs_o^2+va_o^2)*y[i,2] + 2*(cs_o^2)*y[i,2]*(vax_o*y[i,2]+vaz_o*y[i,3])^2 + 2*(cs_o^2)*(y[i,2]^2+y[i,3]^2)*(vax_o*y[i,2]+vaz_o*y[i,3])*vax_o + 2*(cs_o^2)*N2_o*y[i,2],$;dx/dl=∂F/∂p
;    -2*(omg^2)*(cs_o^2+va_o^2)*y[i,3] + 2*(cs_o^2)*y[i,3]*(vax_o*y[i,2]+vaz_o*y[i,3])^2 + 2*(cs_o^2)*(y[i,2]^2+y[i,3]^2)*(vax_o*y[i,2]+vaz_o*y[i,3])*vaz_o + 2*(cs_o^2)*(wc_o^2)*y[i,3],$;dz/dl=∂F/∂q
;    -2*(omg^2)*(cs_o*dcs_dx_o+va_o*dva_dx_o)*(y[i,2]^2+y[i,3]^2) + 2*cs_o*dcs_dx_o*(vax_o*y[i,2]+vaz_o*y[i,3])^2*(y[i,2]^2+y[i,3]^2) + 2*cs_o^2*(vax_o*y[i,2]+vaz_o*y[i,3])*(y[i,2]*dvax_dx_o+y[i,3]*dvaz_dx_o)*(y[i,2]^2+y[i,3]^2) - 2*wc_o*dwc_dx_o*(omg^2-cs_o^2*y[i,3]^2) + 2*wc_o^2*cs_o*dcs_dx_o*y[i,3]^2 + 2*cs_o*dcs_dx_o*N2_o*y[i,2]^2 + cs_o^2*dN2_dx_o*y[i,2]^2,$
;    -2*(omg^2)*(cs_o*dcs_dz_o+va_o*dva_dz_o)*(y[i,2]^2+y[i,3]^2) + 2*cs_o*dcs_dz_o*(vax_o*y[i,2]+vaz_o*y[i,3])^2*(y[i,2]^2+y[i,3]^2) + 2*cs_o^2*(vax_o*y[i,2]+vaz_o*y[i,3])*(y[i,2]*dvax_dz_o+y[i,3]*dvaz_dz_o)*(y[i,2]^2+y[i,3]^2) - 2*wc_o*dwc_dz_o*(omg^2-cs_o^2*y[i,3]^2) + 2*wc_o^2*cs_o*dcs_dz_o*y[i,3]^2 + 2*cs_o*dcs_dz_o*N2_o*y[i,2]^2 + cs_o^2*dN2_dz_o*y[i,2]^2$
;    ]
;   
    
    k1[i,0]=-2*(omg^2)*(cs_o^2+va_o^2)*y[i,2] + 2*(cs_o^2)*y[i,2]*(vax_o*y[i,2]+vaz_o*y[i,3])^2 + 2*(cs_o^2)*(y[i,2]^2+y[i,3]^2)*(vax_o*y[i,2]+vaz_o*y[i,3])*vax_o + 2*(cs_o^2)*N2_o*y[i,2];dx/dl=∂F/∂p
    k1[i,1]=-2*(omg^2)*(cs_o^2+va_o^2)*y[i,3] + 2*(cs_o^2)*y[i,3]*(vax_o*y[i,2]+vaz_o*y[i,3])^2 + 2*(cs_o^2)*(y[i,2]^2+y[i,3]^2)*(vax_o*y[i,2]+vaz_o*y[i,3])*vaz_o + 2*(cs_o^2)*(wc_o^2)*y[i,3];dz/dl=∂F/∂q
    k1[i,2]=(-2*(omg^2)*(cs_o*dcs_dx_o+va_o*dva_dx_o)*(y[i,2]^2+y[i,3]^2) + 2*cs_o*dcs_dx_o*(vax_o*y[i,2]+vaz_o*y[i,3])^2*(y[i,2]^2+y[i,3]^2) + 2*cs_o^2*(vax_o*y[i,2]+vaz_o*y[i,3])*(y[i,2]*dvax_dx_o+y[i,3]*dvaz_dx_o)*(y[i,2]^2+y[i,3]^2) - 2*wc_o*dwc_dx_o*(omg^2-cs_o^2*y[i,3]^2) + 2*wc_o^2*cs_o*dcs_dx_o*y[i,3]^2 + 2*cs_o*dcs_dx_o*N2_o*y[i,2]^2 + cs_o^2*dN2_dx_o*y[i,2]^2)*(-1);-dp/dl
    k1[i,3]=(-2*(omg^2)*(cs_o*dcs_dz_o+va_o*dva_dz_o)*(y[i,2]^2+y[i,3]^2) + 2*cs_o*dcs_dz_o*(vax_o*y[i,2]+vaz_o*y[i,3])^2*(y[i,2]^2+y[i,3]^2) + 2*cs_o^2*(vax_o*y[i,2]+vaz_o*y[i,3])*(y[i,2]*dvax_dz_o+y[i,3]*dvaz_dz_o)*(y[i,2]^2+y[i,3]^2) - 2*wc_o*dwc_dz_o*(omg^2-cs_o^2*y[i,3]^2) + 2*wc_o^2*cs_o*dcs_dz_o*y[i,3]^2 + 2*cs_o*dcs_dz_o*N2_o*y[i,2]^2 + cs_o^2*dN2_dz_o*y[i,2]^2)*(-1);-dq/dl

    k2[i,0]=-2*(omg^2)*(cs_o^2+va_o^2)*(y[i,2]+l*k1[i,2]/2) + 2*(cs_o^2)*(y[i,2]+l*k1[i,2]/2)*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))^2 + 2*(cs_o^2)*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2)*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))*vax_o + 2*(cs_o^2)*N2_o*(y[i,2]+l*k1[i,2]/2)
    k2[i,1]=-2*(omg^2)*(cs_o^2+va_o^2)*(y[i,3]+l*k1[i,3]/2) + 2*(cs_o^2)*(y[i,3]+l*k1[i,3]/2)*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))^2 + 2*(cs_o^2)*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2)*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))*vaz_o + 2*(cs_o^2)*(wc_o^2)*(y[i,3]+l*k1[i,3]/2);dz/dl=∂F/∂q
    k2[i,2]=(-2*(omg^2)*(cs_o*dcs_dx_o+va_o*dva_dx_o)*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2) + 2*cs_o*dcs_dx_o*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))^2*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2) + 2*cs_o^2*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))*((y[i,2]+l*k1[i,2]/2)*dvax_dx_o+(y[i,3]+l*k1[i,3]/2)*dvaz_dx_o)*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2) - 2*wc_o*dwc_dx_o*(omg^2-cs_o^2*(y[i,3]+l*k1[i,3]/2)^2) + 2*wc_o^2*cs_o*dcs_dx_o*(y[i,3]+l*k1[i,3]/2)^2 + 2*cs_o*dcs_dx_o*N2_o*(y[i,2]+l*k1[i,2]/2)^2 + cs_o^2*dN2_dx_o*(y[i,2]+l*k1[i,2]/2)^2)*(-1);-dp/dl
    k2[i,3]=(-2*(omg^2)*(cs_o*dcs_dz_o+va_o*dva_dz_o)*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2) + 2*cs_o*dcs_dz_o*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))^2*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2) + 2*cs_o^2*(vax_o*(y[i,2]+l*k1[i,2]/2)+vaz_o*(y[i,3]+l*k1[i,3]/2))*((y[i,2]+l*k1[i,2]/2)*dvax_dz_o+(y[i,3]+l*k1[i,3]/2)*dvaz_dz_o)*((y[i,2]+l*k1[i,2]/2)^2+(y[i,3]+l*k1[i,3]/2)^2) - 2*wc_o*dwc_dz_o*(omg^2-cs_o^2*(y[i,3]+l*k1[i,3]/2)^2) + 2*wc_o^2*cs_o*dcs_dz_o*(y[i,3]+l*k1[i,3]/2)^2 + 2*cs_o*dcs_dz_o*N2_o*(y[i,2]+l*k1[i,2]/2)^2 + cs_o^2*dN2_dz_o*(y[i,2]+l*k1[i,2]/2)^2)*(-1);-dq/dl 
 
    k3[i,0]=-2*(omg^2)*(cs_o^2+va_o^2)*(y[i,2]+l*k2[i,2]/2) + 2*(cs_o^2)*(y[i,2]+l*k2[i,2]/2)*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))^2 + 2*(cs_o^2)*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2)*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))*vax_o + 2*(cs_o^2)*N2_o*(y[i,2]+l*k2[i,2]/2)
    k3[i,1]=-2*(omg^2)*(cs_o^2+va_o^2)*(y[i,3]+l*k2[i,3]/2) + 2*(cs_o^2)*(y[i,3]+l*k2[i,3]/2)*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))^2 + 2*(cs_o^2)*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2)*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))*vaz_o + 2*(cs_o^2)*(wc_o^2)*(y[i,3]+l*k2[i,3]/2);dz/dl=∂F/∂q
    k3[i,2]=(-2*(omg^2)*(cs_o*dcs_dx_o+va_o*dva_dx_o)*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2) + 2*cs_o*dcs_dx_o*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))^2*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2) + 2*cs_o^2*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))*((y[i,2]+l*k2[i,2]/2)*dvax_dx_o+(y[i,3]+l*k2[i,3]/2)*dvaz_dx_o)*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2) - 2*wc_o*dwc_dx_o*(omg^2-cs_o^2*(y[i,3]+l*k2[i,3]/2)^2) + 2*wc_o^2*cs_o*dcs_dx_o*(y[i,3]+l*k2[i,3]/2)^2 + 2*cs_o*dcs_dx_o*N2_o*(y[i,2]+l*k2[i,2]/2)^2 + cs_o^2*dN2_dx_o*(y[i,2]+l*k2[i,2]/2)^2)*(-1);-dp/dl
    k3[i,3]=(-2*(omg^2)*(cs_o*dcs_dz_o+va_o*dva_dz_o)*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2) + 2*cs_o*dcs_dz_o*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))^2*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2) + 2*cs_o^2*(vax_o*(y[i,2]+l*k2[i,2]/2)+vaz_o*(y[i,3]+l*k2[i,3]/2))*((y[i,2]+l*k2[i,2]/2)*dvax_dz_o+(y[i,3]+l*k2[i,3]/2)*dvaz_dz_o)*((y[i,2]+l*k2[i,2]/2)^2+(y[i,3]+l*k2[i,3]/2)^2) - 2*wc_o*dwc_dz_o*(omg^2-cs_o^2*(y[i,3]+l*k2[i,3]/2)^2) + 2*wc_o^2*cs_o*dcs_dz_o*(y[i,3]+l*k2[i,3]/2)^2 + 2*cs_o*dcs_dz_o*N2_o*(y[i,2]+l*k2[i,2]/2)^2 + cs_o^2*dN2_dz_o*(y[i,2]+l*k2[i,2]/2)^2)*(-1);-dq/dl 
    
    k4[i,0]=-2*(omg^2)*(cs_o^2+va_o^2)*(y[i,2]+l*k3[i,2]) + 2*(cs_o^2)*(y[i,2]+l*k3[i,2])*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))^2 + 2*(cs_o^2)*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2)*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))*vax_o + 2*(cs_o^2)*N2_o*(y[i,2]+l*k3[i,2])
    k4[i,1]=-2*(omg^2)*(cs_o^2+va_o^2)*(y[i,3]+l*k3[i,3]) + 2*(cs_o^2)*(y[i,3]+l*k3[i,3])*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))^2 + 2*(cs_o^2)*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2)*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))*vaz_o + 2*(cs_o^2)*(wc_o^2)*(y[i,3]+l*k3[i,3]);dz/dl=∂F/∂q
    k4[i,2]=(-2*(omg^2)*(cs_o*dcs_dx_o+va_o*dva_dx_o)*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2) + 2*cs_o*dcs_dx_o*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))^2*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2) + 2*cs_o^2*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))*((y[i,2]+l*k3[i,2])*dvax_dx_o+(y[i,3]+l*k3[i,3])*dvaz_dx_o)*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2) - 2*wc_o*dwc_dx_o*(omg^2-cs_o^2*(y[i,3]+l*k3[i,3])^2) + 2*wc_o^2*cs_o*dcs_dx_o*(y[i,3]+l*k3[i,3])^2 + 2*cs_o*dcs_dx_o*N2_o*(y[i,2]+l*k3[i,2])^2 + cs_o^2*dN2_dx_o*(y[i,2]+l*k3[i,2])^2)*(-1);dp/dl
    k4[i,3]=(-2*(omg^2)*(cs_o*dcs_dz_o+va_o*dva_dz_o)*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2) + 2*cs_o*dcs_dz_o*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))^2*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2) + 2*cs_o^2*(vax_o*(y[i,2]+l*k3[i,2])+vaz_o*(y[i,3]+l*k3[i,3]))*((y[i,2]+l*k3[i,2])*dvax_dz_o+(y[i,3]+l*k3[i,3])*dvaz_dz_o)*((y[i,2]+l*k3[i,2])^2+(y[i,3]+l*k3[i,3])^2) - 2*wc_o*dwc_dz_o*(omg^2-cs_o^2*(y[i,3]+l*k3[i,3])^2) + 2*wc_o^2*cs_o*dcs_dz_o*(y[i,3]+l*k3[i,3])^2 + 2*cs_o*dcs_dz_o*N2_o*(y[i,2]+l*k3[i,2])^2 + cs_o^2*dN2_dz_o*(y[i,2]+l*k3[i,2])^2)*(-1);dq/dl 
    
    y[i+1,*]=y[i,*] + l*(k1[i,*]+2*k2[i,*]+2*k3[i,*]+k4[i,*])/6
 
    S[i+1]=S[i] + abs(y[i,2]*(y[i+1,0]-y[i,0])+y[i,3]*(y[i+1,1]-y[i,1]))
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    disp[i]=omg^4 - (omg^2)*(cs_o^2+va_o^2)*(y[i,2]^2+y[i,3]^2) + cs_o^2*(y[i,2]^2+y[i,3]^2)*((vax_o*y[i,2]+vaz_o*y[i,3])^2) - (wc_o^2)*(omg^2-cs_o^2*y[i,3]^2) + cs_o^2*N2_o*y[i,2]^2
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    jx1=round((y[i+1,0]+73440000.0)/96000.0)
    jz1=round((y[i+1,1]+15608800.0)/32000.0)
    k111=-2*(omg^2)*(cs[jx1,jz1]^2+va[jx1,jz1]^2)*y[i,3] + 2*(cs[jx1,jz1]^2)*y[i,3]*(vax[jx1,jz1]*y[i,2]+vaz[jx1,jz1]*y[i,3])^2 + 2*(cs[jx1,jz1]^2)*(y[i,2]^2+y[i,3]^2)*(vax[jx1,jz1]*y[i,2]+vaz[jx1,jz1]*y[i,3])*vaz[jx1,jz1] + 2*(cs[jx1,jz1]^2)*(wc[jx1,jz1]^2)*y[i,3];dz/dl=∂F/∂q
    if (k1[i,1] GT 0) AND (k111 LT 0) AND (y[i,3] LT 0) then begin
        l=5000/(sqrt(k1[i,0]^2+k1[i,1]^2)) 
        print,'gonnalong',i
    endif else begin
        l=5000
    endelse;1️⃣动态调整步长
    
    ;print,i,jzz[i],' next l=',l
endfor

;print,'part'
tmp=0
x_upMm=dblarr(2000)
z_upMm=x_upMm
i_upMm=x_upMm
for i = 1, step-1 do begin
    if (y[i,1] GT y[i-1,1]) and (y[i,1] GT y[i+1,1]) then begin
        x_upMm[tmp]=y[i,0]/1000000.0
        z_upMm[tmp]=y[i,1]/1000000.0
        i_upMm[tmp]=i
        print,'x=',x_upMm[tmp],' z=',z_upMm[tmp],' kx=',y[i,2],' kz=',y[i,3],' i=',i,' jz=',jzz[i]
        tmp=tmp+1
    endif
endfor
;print,'down-turningpoint'
;for i = 1, step-1 do if (y[i,1] LT y[i-1,1]) and (y[i,1] LT y[i+1,1]) then print,'x=',y[i,0]/1000000.0,' z=',y[i,1]/1000000.0,' kx=',y[i,2],' kz=',y[i,3],' i=',i

f=dblarr(511,512)
for i = 0, 510 do begin
    for j = 0, 511 do begin
        if abs(wc[i,j]/(0.001*4.0*2*!pi)-1) LE 0.05 then f[i,j]=1
    endfor
endfor

;ct = reverse(COLORTABLE(72))
;im=image(wc,RGB_TABLE=ct,position=[0.1,0.2,0.9,1]);,max_value=6,min_value=2)
;ip=plot(reform(x)/1000000.0,reform(z)/1000000.0,/nodata,/current,xrange=[min(x),max(x)]/1000000.0,yrange=[min(z),max(z)]/1000000.0,xstyle=1,ystyle=1,position=[0.3,0.2,0.7,1],xtitle='r/Mm',ytitle='z/Mm')
;c=colorbar(target=im)

plot,y[*,0]/1d6,y[*,1]/1d6
print,'dx=',x_upMm[1]-x_upMm[0] 
print,'dt=',(s[i_upMm[1]]-s[i_upMm[0]])/omg  


end