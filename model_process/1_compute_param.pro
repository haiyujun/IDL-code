;二维情况的模拟.x:r黑子径向;z:深度
;注意：！！！！！！全部转化为国际单位制！！！！！！
r=mrdfits('/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/r.fits');单位：Mm
z=mrdfits('/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/z.fits');单位：Mm

fits2map, '/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/br.fits',map
br=map.data*0.0001;1G=0.0001T
fits2map, '/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/bz.fits',map
bz=map.data*0.0001;1G=0.0001T

fits2map, '/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/T.fits',map
T=map.data;K
fits2map, '/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/p.fits',map
p=map.data*0.1;1bar=0.1pa
fits2map, '/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/gamma1.fits',map
gamma1=map.data;无量纲
fits2map, '/Users/lihaiyu/Desktop/IDL code/forward_modeling/model/rho.fits',map
rho=map.data*1000;1g/cm^3=1000kg/m^3

;计算cs,va
cs=sqrt(gamma1*p/rho);计算声速,m/s
cs2=gamma1*p/rho
vax=br/sqrt(0.000001256637*rho)
vaz=bz/sqrt(0.000001256637*rho)
va=sqrt(vax^2+vaz^2);总Alfven速度,m/s
va2=br^2/(0.000001256637*rho)+bz^2/(0.000001256637*rho)

;计算标高，然后计算ωc、N
dp_dz=dblarr(256,512)
drho_dz=dblarr(256,512)
for i = 0, 255 do begin
    dp_dz[i,*]=deriv(reform(z),p[i,*])
    drho_dz[i,*]=deriv(reform(z),rho[i,*])
endfor
H1=-p/dp_dz
H2=-rho/drho_dz;0723:暂时先用密度标高;Mm
;0731:密度标高算出来的N2经常为负，该用压强标高
H2[*,511]=H2[*,510];最顶上一排出现了导数为负的情况，将次最上层延伸至该层
H1[*,510]=H1[*,509]
H1[*,511]=H1[*,509]
H=H1*1000000.0;单位:m

g=6673.0*19890.0/(696.0*696.0) ;6.67430*10^(-11)*1.989*10^30/(696000000)^2
;H3=13806.49*T/(1.75*1.66*g);1.380649*10^(-23)*T/(1.75*1.66*10^(-27)*g)

wc=0.5*cs/H
N=sqrt(g/H-g*g/(cs*cs))
N2=g/H-g*g/(cs*cs)

f=dblarr(256,512)
for i = 0, 255 do begin
    for j = 0, 511 do begin
        if abs(wc[i,j]/(0.001*4.0*2*!pi)-1) LE 0.1 then f[i,j]=1
    endfor
endfor

b=sqrt(br^2+bz^2)
;plot image
ct = reverse(COLORTABLE(72))
im=image(vax^2+vaz^2-va^2,RGB_TABLE=ct,position=[0.1,0.2,0.9,1]);,max_value=6,min_value=2)
ip=plot(reform(r),reform(z),/nodata,/current,xrange=[min(r),max(r)],yrange=[min(z),max(z)],xstyle=1,ystyle=1,position=[0.3,0.2,0.7,1],xtitle='r/Mm',ytitle='z/Mm')
c=colorbar(target=im);,title='c!Ds!N(m/s)')

;cqv=dblarr(256,512)
;for i = 0, 255 do begin
;    for j = 0, 511 do begin
;        if abs(va[i,j]/cs[i,j]-1) LE 0.1 then begin
;            cqv[i,j]=1
;        endif
;    endfor
;endfor
; 
;ct = reverse(COLORTABLE(0))
;im=image(f,RGB_TABLE=ct,position=[0.1,0.2,0.9,1])
;ip=plot(reform(r),reform(z),/nodata,/current,xrange=[min(r),max(r)],yrange=[min(z),max(z)],xstyle=1,ystyle=1,position=[0.3,0.2,0.7,1],xtitle='r/Mm',ytitle='z/Mm')

end