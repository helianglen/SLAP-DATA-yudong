
# See how calibrated Tb over land varies with assumed Tb over lake 

# Tb
skyTv= 10   # K
skyTh= 10   # K

#Counts 
skyCv= 2.34E6  
skyCh= 2.55E6  

#over soil counts
soilCv = 3.22E6
soilCh = 3.43E6

#foambox counts 
foamCv = 3.54E6
foamCh = 3.68E6

#lake counts
lakeCv = 2.78E6
lakeCh = 2.86E6

#lake Tb range, assuming no polariation 
lakeTv=50:150
lakeTh=50:150

#
slope_v = (lakeCv - skyCv) / (lakeTv - skyTv) 
slope_h = (lakeCh - skyCh) / (lakeTh - skyTh) 

soilTv = skyTv + (soilCv - skyCv) / slope_v 
soilTh = skyTh + (soilCh - skyCh) / slope_h 

foamTv = skyTv + (foamCv - skyCv) / slope_v 
foamTh = skyTh + (foamCh - skyCh) / slope_h 

x11(bg='white', width=10, height=15)
par(mfrow=c(3, 1))
#par(mfg=c(1, 1, 1, 1))

plot(lakeTv, soilTv, type='l', col='blue', xlab='Lake Tb (K)', ylab='Tb (K)', 
     main='No counts-swap', lwd=2)
lines(lakeTh, soilTh, col="green", lty=1, lwd=2)
lines(lakeTv, foamTv, col="blue", lty=2, lwd=2)
lines(lakeTv, foamTh, col="green", lty=2, lwd=2)
legend('topleft', c("soil TbV", "soil TbH", "foam TbV", "foam TbH"), 
       col=c('blue', 'green', 'blue', 'green'), lty=c(1, 1, 2, 2))

#savePlot('soilTb_vs_lakeTb.png') 

# swap all counts 
tmp=skyCv
skyCv=skyCh
skyCh=tmp 

tmp=soilCv
soilCv=soilCh
soilCh=tmp 

tmp=lakeCv
lakeCv=lakeCh
lakeCh=tmp 

tmp=foamCv
foamCv=foamCh
foamCh=tmp 

slope_v = (lakeCv - skyCv) / (lakeTv - skyTv)
slope_h = (lakeCh - skyCh) / (lakeTh - skyTh)

soilTv = skyTv + (soilCv - skyCv) / slope_v
soilTh = skyTh + (soilCh - skyCh) / slope_h

foamTv = skyTv + (foamCv - skyCv) / slope_v 
foamTh = skyTh + (foamCh - skyCh) / slope_h 

#x11(bg='white', width=8, height=6)

plot(lakeTv, soilTv, type='l', col='blue', xlab='Lake Tb (K)', ylab='Tb (K)', 
     main='All counts swapped', lwd=2)
lines(lakeTh, soilTh, col="green", lty=1, lwd=2)
lines(lakeTv, foamTv, col="blue", lty=2, lwd=2)
lines(lakeTv, foamTh, col="green", lty=2, lwd=2)
legend('topleft', c("soil TbV", "soil TbH", "foam TbV", "foam TbH"), 
       col=c('blue', 'green', 'blue', 'green'), lty=c(1, 1, 2, 2))

#savePlot('soilTb_vs_lakeTb-all-swapped.png')


# swap sky only counts
#tmp=skyCv
#skyCv=skyCh
#skyCh=tmp

# other swapped back to original 
tmp=soilCv
soilCv=soilCh
soilCh=tmp

tmp=lakeCv
lakeCv=lakeCh
lakeCh=tmp

tmp=foamCv
foamCv=foamCh
foamCh=tmp

slope_v = (lakeCv - skyCv) / (lakeTv - skyTv)
slope_h = (lakeCh - skyCh) / (lakeTh - skyTh)

soilTv = skyTv + (soilCv - skyCv) / slope_v
soilTh = skyTh + (soilCh - skyCh) / slope_h

foamTv = skyTv + (foamCv - skyCv) / slope_v
foamTh = skyTh + (foamCh - skyCh) / slope_h


plot(lakeTv, soilTv, type='l', col='blue', xlab='Lake Tb (K)', ylab='Tb (K)',
     main='Sky counts swapped', lwd=2)
lines(lakeTh, soilTh, col="green", lty=1, lwd=2)
lines(lakeTv, foamTv, col="blue", lty=2, lwd=2)
lines(lakeTv, foamTh, col="green", lty=2, lwd=2)
legend('topleft', c("soil TbV", "soil TbH", "foam TbV", "foam TbH"),
       col=c('blue', 'green', 'blue', 'green'), lty=c(1, 1, 2, 2))

savePlot('soilTb_vs_lakeTb-w-swap.png')


