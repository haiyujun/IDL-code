a = mrdfits('/Users/lihaiyu/Downloads/Oct10_2_V.fits')
;contour,a(*,*,0),/cell_fill
ft=complexarr(640)
for i = 0, 99 do begin
    for j = 0, 99 do begin
        ft=ft+fft(a(i,j,*))
    endfor
endfor
ft=ft/10000
plot,findgen(320)*1000/(640*45),abs(ft),xtitle='frequncy(mHz)',ytitle='amplitude(m/s)',max_value=40
end