data = mrdfits('/Users/lihaiyu/Downloads/Oct10_2_V.fits')
loadct,3
b=data[100,*,*]
ftb=reform(abs(fft(b)))
tvim,ftb[0:255,0:319],range = [0,8,.1],xrange = [0,255]/(12.15*0.03*512),yrange = [0,319]*2*!pi/(640*45)
end