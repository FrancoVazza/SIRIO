
;....SIRIO
;....by Franco Vazza, Annalisa Bonafede (orig. 2012, new version 2017)
;....takes as input a simulated radio image [erg/(s Hz)]
;....computes the observable flux muJ/pixel at distance dL
;....takes as input a dataset for a given interferometer configuration
;....simulates filling in the uv plane due to earth rotation (angles are
;.......treated in an over-simplifed way!
;....computes UV coverage and dirty beam
;....clean the dirty image via iterations
;....IMPORTANT: this basic version assumes that pixel=beamsise in the input sky model!
;
;
;...calling sequence 
;...sirio,file_sky_model='input.fits',number_point_sources=20,live_plot='y', sigma_rms=1e-6     etc...

pro sirio,file_sky_model=file_sky_model,folder=folder,file_antenna=file_antenna,number_point_sources=number_point_sources,live_plot=live_plot,max_iter=max_iter,sigma_rms=sigma_rms,hour=hour

print,'>>STARTING SIRIO<<'

  if not IS_DEF(file_sky_model) then file_sky_model='relic1.fits'  ;...input sky model 
   if not IS_DEF(folder) then folder='/Users/francovazza/Downloads/SIRIO/'   ;...main folder containing SIRIO
   if not IS_DEF(file_antenna) then file_antenna='VLA.reg'   ;...file with atenna positions
   if not IS_DEF(number_point_sources) then number_point_sources=2   ;...additional pointlike sources to be generated
   if not IS_DEF(live_plot) then live_plot='y'  ;...plots on 'x' device while computing 
   if not IS_DEF(max_iter) then max_iter=1e4   ;..max iterations in cleaning
   if not IS_DEF(sigma_rms) then sigma_rms=8e-5 ;...Jy/beam
   if not IS_DEF(hour) then hour=10  ;...hours of integration
   

print,'Reading sky model' 
   imag=read_image_fits(folder,file_sky_model)
   
   ;imag=read_image_bin(n0,folder,namef)   alternative reader for binary 2D datasets
   ns=size(imag)
   n0=ns(1)
   print,'the size of the input image is ',ns(1:2) 

print,'the number of additional point sources is ',number_point_sources
  if number_point_sources gt 0 then begin
  imas=generate_radio_gal(number_point_sources,n0)
  norm_gal=1. ;..normalisation factor for radio galaxy emission [erg/s/Hz], can be changed for different purposes
  imag+=norm_gal*imas  
  writefits,folder+'/output/sky_model_with_sources.fits',imag
  endif


 if live_plot eq 'y' then begin
  print,'live plotting of all main steps in SIRIO'
  print,'please notice for > 256x256 images the live plotting introduces a slow-down of everything'
  set_plot,'x' ;...show the image 
  !p.multi=[0,2,4]
  !p.font=1
  window,1,xsize=600,ysize=1200
  contour,sqrt(imag),nlevels=128,/fill,title='input sky model',charsize=3
  endif

 print,'writing the sky model as .fits file in /output'
 writefits,folder+'/output/input_map.fits',abs(imag)

 print,'FFT transforming the dataset' 
 imaff=shift(fft(imag,1),n0*0.5-1,n0*0.5-1)  ;...simple FFT transform of the image, centre in 0.5*n-1,0.5*n-1 (the conversion in idl is to centre it at 0,0)
 if live_plot eq 'y' then contour,sqrt(abs(real_part(imaff))),/fill,nlevels=128,title='input sky model - FFT transformed',charsize=3
 
 a=read_antenna(folder,file_antenna,n_ant,x,y)

 if live_plot eq 'y' then plot,x,y,title='antenna positions -JVLA',psym=4,xtitle='km',ytitle='km',charsize=3

  a=generate_uv(n0,hour,uv,x,y,uv_t,rmax)

  if live_plot eq 'y' then contour,uv_t,title='UV coverage (visibilites)',charsize=3
  writefits,folder+'/output/uv_time_t.fits',uv_t

 beam=generate_beam(uv_t,n0)   ;...generate beam of the radio observation 

 if live_plot eq 'y' then contour,sqrt(beam),nlevels=10,title='beam',charsize=3
 writefits,folder+'/output/beam.fits',beam
                                          ;...we need the usual shift for the FFT in IDL
 imaff=imaff*uv_t ;...here the image in the Fourier space is convolved by the beam
 ;...THIS IS THE MAP OF VISIBILITES (IMA-FFT CONVOLVED BY THE BEAM)
 writefits,folder+'/output/visib_map.fits',abs(real_part(imaff))

if live_plot eq 'y' then contour,sqrt(abs(imaff)),/fill,nlevels=128,title='visibilities (ima-FFT x beam)',charsize=3

  ; GRIDDING IN THE UV PLANE 
  ; imaff=gridding(imaff,rmax)
   ;....under development, not clear what interpolation technique the real radio observation use

 a=show_visib(folder,imaff,n0)

 noise=add_noise(sigma_rms,n0)
  writefits,folder+'/output/noise.fits',noise

 
  imac=create_dirty(folder,imaff,noise)
  if live_plot eq 'y' then contour,sqrt(abs(imac)),nlevels=128,/fill,title='dirty image',charsize=3
  writefits,folder+'/output/dirty_image.fits',imac

 
  rima=imac
  a=cleaning_iterations(ima_clean,rima,ima_components,beam,n0,sigma_rms,max_iter,live_plot)
  

;....here we produce fits files for the restored image, the residual map and the components
 
if live_plot eq 'y' then begin

set_plot,'x' 
print,'showing the final result'
!p.multi=[0,3,0]
!p.font=1
window,12,xsize=1200,ysize=400
tvscl,[ima_clean(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1),rima(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1)]
;contour,((ima_clean(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1))),nlevels=128,/fill,title='clean image',charsize=3
;contour,((rima(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1))),nlevels=128,/fill,title='residuals',charsize=3
;contour,((ima_components(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1))),nlevels=128,/fill,title='components',charsize=3
 endif
 
 writefits,folder+'/output/ima_clean_final.fits',ima_clean(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1)
 writefits,folder+'/output/ima_residuals.fits',(rima(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1))
 writefits,folder+'/output/ima_components.fits',(ima_components(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1))
 
print,'run done'
print,'>>all data written in /output<<'

end


function uv_rotate,uv_t,hour
;...dummy function which makes the assumed array rotating, but applies
;...only if the array is at the north pole and the source is close to the zenith

s=size(uv_t)
n=s(1)
n_ang=360       ;...number of angular elements to sample the rotation
dt=2*!pi/float(n_ang) ;...corresponding angular resolution
       
obs_angle=2.*!pi*hour/float(24.) ;...total earth rotation during observation
tmax=uint(obs_angle/float(dt))
uv=uv_t
for t=1L,tmax do begin   
ang=360*t*dt/float(2*!pi)  ;...total rot.angle since begin. of observation
uv_t=uv_t+rot(uv,ang)
endfor
return,uv_t
end


function uv_rotate_real,uv_t,hour,del,lat

;....this routine would provide the correct way of filling the UV plane
;....making use of declination of the source, latidue and hour angle
;...still not working...

hour=1 ;....integration time

u=x
v=y
obs_angle=2*!pi*12/float(24)  ;...staring observation angle

u=uvcx+((u-uvcx))*(sin(obs_angle)+cos(obs_angle))
v=uvcy+((v-uvcy))*(-1*sin(del)*cos(obs_angle)+sin(obs_angle)*sin(del))

window,6
plot,u,v,psym=4,xrange=[0,2*n0-1],yrange=[0,2*n0-1]
oplot,x,y,psym=2

for t=0,1e2 do begin
obs_angle=obs_angle+2*!pi/float(9e2)
print,obs_angle

ut=uvcx+(u-uvcx)*(sin(obs_angle)+cos(obs_angle))
vt=uvcy+(v-uvcy)*(-1*sin(del)*cos(obs_angle)+sin(obs_angle)*sin(del))

u=ut
v=vt
oplot,u,v,psym=3
endfor


end


function gridding,imaff,rmax

;...this is the gridding procedure to fill in the UV plane
;...it is not working properly, so it is skipped for the moment
;...unluckily, in IDL there are plenty of routines to choose, not sure
;...of what is actually done in the real radio reduction
;...it assumes the centre of the UV plane is at [n0*0.5-1,n0*0.5-1]
s=size(imaff)
n0=s(1)
iuv=where(imaff ne 0,nn)
puv=complexarr(nn,3) ;...array of final covered "pixels" in the uv plane
                 ;...done in this way to use easily the IDL function grid_tps 
for i=0,nn-1 do begin
x=array_indices(imaff,iuv(i))
puv(i,0)=x(0)
puv(i,1)=x(1) 
puv(i,2)=imaff(x(0),x(1))
endfor

;imaffr=grid_tps(puv(*,0),puv(*,1),puv(*,2),ngrid=[n0,n0],start=[0,0],delta=[1,1])

;imaffr=grid_tps(puv(*,0),puv(*,1),float(puv(*,2)),ngrid=[n0,n0],start=[0,0],delta=[1,1])
;imaffc=grid_tps(puv(*,0),puv(*,1),imaginary(puv(*,2)),ngrid=[n0,n0],start=[0,0],delta=[1,1])

imaffr=griddata(puv(*,0),puv(*,1),float(puv(*,2)),dimension=n0,start=[0,0],delta=[1,1],function_type=0,smoothing=2)
imaffc=griddata(puv(*,0),puv(*,1),imaginary(puv(*,2)),dimension=n0,start=[0,0],delta=[1,1],function_type=0,smoothing=2)

imaff=complex(imaffr,imaffc)
;...we do not allow grid_tps to extrapolate outside the disk
for l=0,n0-1 do begin
for i=0,n0-1 do begin
 if sqrt((i-n0*0.5+1)^2.+(l-n0*0.5+1)^2.) gt rmax then imaff(i,l)=complex(0,0)
endfor
endfor
contour,sqrt(abs(imaff)),/fill,nlevels=128,title='visibilities (ima-FFT x beam), gridded',charsize=3

return,imaff
end


function read_image_bin,n0,folder,namefile
imag=fltarr(n0,n0)  ;...input image
file_input=folder+'/input/'+namefile
openr,4,file_input  ;...the input image is read
readu,4,imag
close,4
return,imag
end


function read_image_fits,folder,namefile
 imag=readfits(folder+'/input/'+namefile,/noscale,h)
 return,imag
end



function read_antenna,folder,name_ant,n_ant,x,y


  ;..............read in the configuration of a given array of n_ant antennas
  ;.....warning: not a real configuration file, the position of antennas is just approximate!

  file_input=folder+'/input/'+name_ant
 
  b=1
  c=1
 
  close,4
  openr,4,file_input
  readf,4,n_ant
  
  x=intarr(n_ant)
  y=intarr(n_ant)
  for a=0,n_ant-1 do begin
    readf,4,b,c
    x(a)=b
    y(a)=c
  endfor
  close,4

  ;....positions are referred to the centre of the array 
  x=x-total(x)/float(n_ant)
  y=y-total(y)/float(n_ant)



end


function generate_uv,n0,hour,uv,x,y,uv_t,rmax

uv=fltarr(n0,n0)  ;....the UV plane is initialized here to fill with antennas projected positions
;....we need som juggling to convert the input array into the UV plane
uv(*,*)=0.

rmax=max(sqrt(x(*)^2.+y(*)^2.)) ;..maximum radius of the array, it must be normalized to the real km scale

;....now the centre of the array is in [0.5*n0-1,0.5*n0-1]
uvcx=n0*0.5-1. ;..centre of the UV domain
uvcy=n0*0.5-1.
x=uvcx+n0*0.5*(x)/float(rmax)-1
y=uvcy+n0*0.5*(y)/float(rmax)-1
ns=size(x)
n_ant=ns(1)

;....in this way, the starting location of the array is expressed in the plane of
;....the image [0,n0]x[0,n0], with each pixel having about the right angular resolution
;....of the VLA-D  25 " = 25 kpc
;...this has to be changed for images produced at different resolution / redshift

lat=34 ;is the observer latitude in degrees
del=75 ;is the source declination in degrees. since we deal with simulation, we can assume
; the axes of the image are "confortably" aligned with the image plane, so the latitude and the earth
; rotation are providing the only angles
; JUST A SIMPLIFYING STARTING ASSUMPTION, WE NEED TO REFINE ON THIS!

lat=2*!pi*lat/float(360.) ;..turns to radians
del=2*!pi*del/float(360.)              ;..."


y(*)=uvcy+sin(lat)*(y(*)-uvcy) ;...assume some projection at the latitude L
;...not what happens in real observation, just a start

for a=0,n_ant-1 do begin ;...we fill the UV plane with the initial configuration here
  b=uint(x(a))
  c=uint(y(a))
  uv(b,c)=1
endfor

xi=x    ;....the initial position of antennas in the UV plane are stored here as xi and yi
yi=y
x=0
y=0

;writefits,folder+'/out/uv_plane_t0.fits',uv

;....TIME-DEPENDENT UV COVERAGE IS SIMULATED HERE

uv_t=uv ;...the time dependent UV coverage is initialized here
;..we fill the UV coverage at the end of the observation time
uv_t=uv_rotate(uv_t,hour)
i2=where(uv_t gt 0) ;...1 only within the UV coverage of the observation, 0 elsewhere
uv_t(i2)=1.

;...we compute the maximum extent of the uv disk here.
;...necessary below to avoid extrapolation outside of this disk in the "gridding"
rmax=0
r2=fltarr(2)
for l=0,n0-1 do begin
  for i=0,n0-1 do begin
    if uv_t(i,l) gt 0 then begin
      rr=sqrt((i-n0*0.5+1)^2.+(l-n0*0.5+1)^2.)
      r2(0)=rr
      r2(1)=rmax
      rmax=max(r2)
    endif
  endfor
endfor

end


function generate_beam,uv_t,n0
  ;....BEAM CREATION
  beam=fft(shift(uv_t,-n0*0.5+1,-n0*0.5+1),-1)  ;...the beam in real space is created here
;  beam=fft(uv_t,-1)  ;...the beam in real space is created here

  beam=shift(beam,n0*0.5,n0*0.5)
  for l=0,n0-1 do begin
    for i=0,n0-1 do begin
      rr=sqrt((i-n0*0.5+1)^2.+(l-n0*0.5+1)^2.)
      beam(i,l)=beam(i,l)/float(1.+5e-2*rr) ;..trick to reduce edge effects in the outer beam pattern
    endfor
  endfor
  return,beam
  end

function show_visib,folder,imaff,n0

  ;....plot visibilities as in radio data reduction
;  set_plot,'ps'
;  !p.multi=0
;  device,filename=folder+'output/visib_radio.ps',/color
  o=indgen(2)
  loadct,13
  plot,o,o,/nodata,xrange=[1,70],yrange=[1,1e6*max(imaff)],xtitle='baselines', ytitle='flux [!4l!6Jy/baseline]',charsize=3,title='visibilities',/ylog
  for k2=0,n0-1 do begin
    for k1=0,n0-1 do begin
      rk=(1+sqrt((k2-(n0*0.5-1))^2.+(k1-(n0*0.5-1))^2.))
      mo=(float(imaff(k1,k2)))^2.+(imaginary(imaff(k1,k2)))^2.
     
      plots,rk,1e6*sqrt(real_part(mo)),psym=3,noclip=0
    
    endfor
  endfor
 ; device,/close
  
  end
  
  
  function add_noise,sigma_rms,n0
  
    noise=sigma_rms*randomn(12312313,n0,n0)  ;...generate 2D map of noise (assumed normal distribution)
    
    return,noise
    end
    
  function create_dirty,folder,imaff,noise
  print,'creating dirty image '
 
  imagk=fft(imaff,-1)                   ;...and FFT transformed in real space image

  imac=sqrt((real_part(imagk))^2.+(imaginary(imagk))^2.)+noise

  return,imac
  end


function cleaning_iterations,ima_clean,rima,ima_components,beam,n0,sigma_rms,max_iter,live_plot

  
  ;....to mimic real data reduction, we need to paste our image into a 4x4 larger domain,
  ;....since we need to subtract the beam even at the edges
  
  ima=fltarr(2*n0,2*n0)
  ima(*,*)=0
  ima(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1)=rima ;...2 times larger image
  imadum=(beam)                            ;...we create a temporary image to paste beam and original
  ;...image into the 4x4 larger volume
  beaml=fltarr(2*n0,2*n0)
  beaml(*,*)=0
  beaml(0.5*n0:1.5*n0-1,0.5*n0:1.5*n0-1)=imadum  ;...2 times larger beam image
  imadum=0

  ;......CLEANING ALGORITHM
  print,"NOW CLEANING DIRTY IMAGE WITH", max_iter,' iterations'
  rima=ima  ;...maps of residuals (initialized here)

  max_beam=max(beaml)
  pi_beam=where(beaml eq max_beam)
  pi_beam2=array_indices(beaml,pi_beam)  ;...we ensure that locate centre of the beam image is the centre of the domain
  beaml=abs(1.*beaml/float(max_beam)) ;...beam is 1 at his maximum


if live_plot eq 'y' then begin
  set_plot,'x'
  window,6,xsize=2*n0,ysize=n0
endif
  ima_clean=fltarr(2*n0,2*n0) ;...map of components
  ima_clean(*,*)=0.

  it_max=max_iter

  x1=0.5*n0
  x2=1.5*n0-1
  g=0.03 ;...loop gain

  for it=0L,it_max do  begin
  if it/float(1000) eq uint(it/float(1000)) then print,'done ',it,' iterations'
    mi=max((rima(x1:x2,x1:x2)))   ;...maximum in the maps of resiudal
    pix=where(rima eq mi)         ;...locate maximum

    pix2=array_indices(rima,pix)

    ibeam=g*mi*shift(beaml,pix2(0)-pi_beam2(0),pix2(1)-pi_beam2(1))  ;...the beam pattern is multiplied by the maximum and by g
    ima_clean(pix2(0),pix2(1))=ima_clean(pix2(0),pix2(1))+g*mi ;...the map of components is updated by this contribution
    rima=rima-ibeam ;...the map of residuals is reduced by this amount

    if it/float(100) eq uint(it/float(100)) then begin
      ;...here we compute run-time distribution of total real flux, flux in the dirty
      ;...image and flux in the component map

  if live_plot eq 'y' then    tvscl,[rima(x1:x2,x1:x2),ima_clean(x1:x2,x1:x2)]

     ; print,it
    endif


    if rima(pix2(0),pix2(1)) lt 1.1*sigma_rms then break ;..we exit from the loop when the maximum in the residual map>3 sigma noise
  endfor

  imag:
  print,'end of cleaning after',it,'   iterations'

  ;.....WE MUST CONVOLVE HERE THE MAP OF COMPONENTS WITH A GAUSSIAN OF FWHM=TO THAT OF THE PRIMARY BEAM

  beam_fwhm=2 ;...needs to be refined
  ima_components=ima_clean
  ima_clean=smooth(ima_clean,beam_fwhm)+rima

end
  
  
  
  
function generate_radio_gal,number_point_sources,n0

imas=fltarr(n0,n0)
ngal=number_point_sources

xc=n0*randomu(12313,ngal)  ;..random x and y positions
yc=n0*randomu(1255313,ngal)
ig=randomn(1141,ngal)   ;...random distribution of fluxes

normg=0.1/float(max(ig))   ;..we normalize the sources so that their max is 0.1 Jy/beam @1.4 GHz
ig*=normg   

for i=0,ngal-1 do begin
imas(xc(i),yc(i))+=ig(i)
endfor

imas=smooth(imas,1) ;...the sources are smoothed to make them roundish

return,imas

end
