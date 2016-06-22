pro rh1d_fig_iform,l0,onscreen=onscreen,ptitle=ptitle,debug=debug,vrange=vrange,$
                 trange=trange,zrange=zrange,irange=irange,scalecf=scalecf

if n_params() lt 1 then begin
print,'pro rh1d_fig_iform,l0,onscreen=onscreen,ptitle=ptitle,debug=debug,vrange=vrange,'
print,'                 trange=trange,zrange=zrange,irange=irange,scalecf=scalecf'
return
endif


@input.common
@geometry.common
@files.common
@atmos.common
@spectrum.common
@opacity.common

read_rh2,'./',sp,atom,atmos,geometry


; lambda
lambda=sp.lambda
nl=n_elements(lambda)

;read intensity
intensity=sp.i[*,4]

;ndep
nz=geometry.ndep

; read chi,snu,jnu
chi=fltarr(nl,nz)
snu=fltarr(nl,nz)
jnu=fltarr(nl,nz)
for i=0,nl-1 do begin
   opacandsource,waveno=i,ray=4,tmp1,tmp2,jnu=tmp3
   chi[i,*]=tmp1
   snu[i,*]=tmp2
   jnu[i,*]=tmp3
endfor

z=geometry.height
tg=atmos.t

vz=geometry.vz

filename='iform.eps'
if keyword_set(ptitle) then $
   title='(ix,iy)=('+string3(ix)+','+string3(iy)+')' $
else $
   title=''

if not keyword_set(vrange) then vrange=[-50,50]
if not keyword_set(trange) then trange=[2,10]
if not keyword_set(zrange) then zrange=[0,max(z/1e3)]

fig_iform,z,chi,snu,vz,tg,intensity,lambda,l0,vrange=vrange,filename=filename,$
          trange=trange,zrange=zrange,onscreen=onscreen,title=title,debug=debug,$
          irange=irange,scalecf=scalecf

if keyword_set(debug) then stop

end

pro fig_iform,height2,chi2,snu2,vz2,temp2,ie2,ll,l0,$
              zrange=zrange,vrange=vrange,debug=debug,$
              trange=trange,irange=irange,$
              title=title,novel=novel,filename=filename,$
              onscreen=onscreen,scalecf=scalecf
;
;jorritl: May 16 2012 
;makes a 4-panel figure,explaining the vertical intensity as in
; Carlsson & Stein,1997,ApJ.,481,500
;
; height            [m]
; chi(lambda,height)[m^-1]
; snu(lambda,height)[J s^-1 m^-2 Hz^-1 sr^-1]
; ie(lambda)        [J s^-1 m^-2 Hz^-1 sr^-1]
; vz(lambda,height) [m s^-1]
; ll                [nm]
; l0                [nm]
;
; zrange gives the plotting range of the z axis, whole domain is
; plotted if absent. (in km)
; vrange is the same for the frequency scale (in km/s)
; trange is the same for S and T in [K]
; irange is the same for profile, in [K] converted from [J s^-1 m^-2 Hz^-1 sr^-1]
; title: adds a title in the upper left corner of the image
; novel suppresses plotting of the vertical velocity
;-


if (n_params(0) lt 8) then begin
   print,' fig_iform,height,chi2,snu2,vz2,ie2,ll,l0,$'
   print,'   zrange=zrange,vrange=vrange,debug=debug,$'
   print,'   trange=trange,irange=irange,$'
   print,'   title=title,novel=novel,,filename=filename,$'
   print,'   onscreen=onscreen,scalecf=scalecf'
   return
endif

; constant definitions
CLIGHT     = 2.99792458E+08
m_to_km=1e-3
km_to_mm=1e-3
k2kk=1e-3
labelsize=1.2

; set filename
if not keyword_set(filename) then filename='fig_iform.eps'

; all arrays to expected number of dimensions
height=height2
chi=reform(chi2)
snu=reform(snu2)
vz =reform(vz2 )
ie=reform(ie2)
temp=reform(temp2)

; plot velocity curve or not?
plotvz=1
if keyword_set(novel) then plotvz=0

z=height*m_to_km
if not keyword_set(zrange) then zrange=[min(z),max(z)]

; compute tau scale
sz=size(chi)
nl=sz[1]
nv=nl
nz=sz[2]
tau=fltarr(nl,nz)
for i=0,nl-1 do begin
   vertint,chi[i,*],height,dum
   tau[i,*]=dum
endfor

; velocity scale,negative is blueshift [km/s]
dv = clight*(ll-l0)/l0*m_to_km
if not keyword_set(vrange) then vrange=[min(vz),max(vz)]

;velocity, positive is an upflow [km/s]
vz*=m_to_km

 ; temperature in kK
temp*=k2kk
if not keyword_set(trange) then trange=[min(temp),max(temp)]


; chi/tau [m-1]
xt=chi/tau
xt[*,0]=0 ;avoid spurious value at top of atmopshere.

; tau*exp(-tau)
tet=tau*exp(-tau) 

;total contribution function
cfc=xt*tet*snu
if keyword_set(scalecf) then for i=0,nv-1 do cfc[i,*]/=max(cfc[i,*])   ;-trapez(z,cfc[i,*])

; source function as radiation temperature [kK]
for i=0,nl-1 do snu[i,*]=tradiation(snu[i,*],ll[i])*k2kk

; tau=1 height [km]
t1h=fltarr(nv)
for i=0,nv-1 do t1h[i]=interpol(height,reform(tau[i,*]),1.)*m_to_km

;source function at t1h
st1=fltarr(nv)
for i=0,nv-1 do st1[i]=interpol(snu[i,*],reform(tau[i,*]),1.)

; split in blue and red part
iib=max(where(dv lt 3*vrange[0]))
iir=min(where(dv gt 3*vrange[1]))
iic=where(t1h eq max(t1h[iib:iir]))

; fine height grid to account for sudden jumps in t1h
np=200
t1hb=findgen(np)/float(np-1)*(t1h[iic]-t1h[iib])[0] + t1h[iib]
t1hr=findgen(np)/float(np-1)*(t1h[iic]-t1h[iir])[0] + t1h[iir]

;interpolate corresponding dv on fine height grid
dvb=interpol(dv[iib:iic],t1h[iib:iic],t1hb)
dvr=interpol(dv[iic:iir],t1h[iic:iir],t1hr)

; source function at tau=1 height
xint=interpol(findgen(nv),dv,dvb)
yint=interpol(findgen(nz),z,t1hb)
st1b=fltarr(np)
for i=0,np-1 do st1b[i]=bilinear(snu,xint[i],yint[i])

xint=interpol(findgen(nv),dv,dvr)
yint=interpol(findgen(nz),z,t1hr)
st1r=fltarr(np)
for i=0,np-1 do st1r[i]=bilinear(snu,xint[i],yint[i])

; intensity in radiation temperature
for i=0,nv-1 do ie[i]=tradiation(ie[i],ll[i])*k2kk
if not keyword_set(irange) then irange=[min(ie),max(ie)]

; set z axis to Mm
z*=km_to_mm
zrange*=km_to_mm
t1hb*=km_to_mm
t1hr*=km_to_mm
t1h*=km_to_mm

; reverse v axis direction
vrange=reverse(vrange)
dv=-dv

; data ready, now plot

; ---------------------------------------------------------------
; main plot parameters
     nxx=2                      ; number of plots in x-direction
     nyy=2                      ; number of plots in y-direction
     xleft=2.5                  ; left margin in cm
     xright=2.5                 ; right margin in cm
     ybottom=2.5                ; bottom margin in cm
     ytop=2.0                   ; bottom margin in cm
     xpicsz=8.                  ; x-size of single image in cm
     ypicsz=8                   ; y-size of single image in cm
     xdist=0.4                  ; distance between plots in x-direction in cm
     ydist=0.4                  ; distance between plots in y-direction in cm
     xtot=xleft+nxx*xpicsz+(nxx-1)*xdist+xright
     ytot=ybottom+nyy*ypicsz+(nyy-1)*ydist+ytop
     aspect=ytot/xtot
; ---------------------------------------------------------------
; calculate normalized plot window positions
     calc_position,v_pos,v_cpos,xs,ys,$
                   nx=nxx,ny=nyy,$
                   xplotsz=xpicsz,xleft=xleft,xright=xright,xdist=xdist,$
                   yplotsz=ypicsz,ytop=ytop,ybottom=ybottom,ydist=ydist

; ---------------------------------------------------------------
; graphical niceties
     !p.ticklen=-0.02
     !p.thick=2
     empty=replicate(' ',30)
; color table index.                                                            
     rainbow=13
     BW=0
     loadct,bw,/silent
; ---------------------------------------------------------------

if keyword_set(onscreen) then begin
   set_plot,'x'
   !p.font=1
   window,/free,xs=1024,ys=1024
   bgblack=1
   !p.charsize=3
   labelsize=4
endif else begin
   bgblack=0
  !p.font=1
  

; open file
     outpath='./'
     set_plot,'ps'
     device, filename=outpath+filename,$
             bits_per_pixel=8,font_size=10, set_font='Helvetica',/encapsulated,$
             /tt_font,/color,xsize=12,ysize=12*aspect

  endelse

; chi/tau
     im=bytscl(xt)
     x=dv
     y=z
     mplot_image,im,x,y,$
                position=v_pos[0,*],$
                xtitle=' ',ytitle='z [Mm]',xtickname=empty,$
                xthick=2,ythick=2,/noerase,$
                 xrange=vrange,yrange=zrange,$
                 xstyle=1,ystyle=1,bgblack=bgblack
     loadct,13,/silent
     oplot,dv,t1h,color=250,linestyle=2  
     loadct,0,/silent
     if plotvz then oplot,vz,y,color=250
     oplot,[0,0],[min(y),max(y)],color=255,thick=1 

     xyouts,0.95*vrange[0],0.9*zrange[1], '!9c!dn!n!x/!9t!dn!n!x',$
            color=255,charsize=labelsize

; Snu
     im=bytscl(snu > trange[0] < trange[1])
     mplot_image,im,x,y,$
                position=v_pos[1,*],$
                xtitle=' ',ytitle=' ',xtickname=empty,ytickname=empty,$
                xthick=2,ythick=2,/noerase,xstyle=9,$
                 xrange=vrange,yrange=zrange,$
                 ystyle=1,bgblack=bgblack
     loadct,13,/silent
     oplot,dv,t1h,color=250,linestyle=2   ; tau=1 height
     loadct,0,/silent
     if plotvz then oplot,vz,y,color=250 ; vertical velocity

     xyouts,0.95*vrange[0],0.9*zrange[1], 'S!d!9n!n!x',color=255,charsize=labelsize
  
      plot,reverse(trange),zrange,position=v_pos[1,*],/noerase,$
          xstyle=13,ystyle=13,/nodata,xrange=reverse(trange)
     oplot,temp,y,line=1,color=255,thick=1
loadct,13,/s
oplot,st1b,t1hb,color=50
oplot,st1r,t1hr,color=250
loadct,0,/s

;     oplot,snu[nv/2,*],y,line=2,color=255,thick=1
     axis,0,zrange[1],xaxis=1,xstyle=1,xtitle='T [kK]',xticklen=-0.02
  

; tau * exp(-tau)
     im=bytscl(tet)
     mplot_image,im,x,y,$
                position=v_pos[2,*],$
                xtitle='!9Dn!x [km/s]',ytitle='z [Mm]',$
                xthick=2,ythick=2,/noerase ,$
                xrange=vrange,yrange=zrange,xstyle=1,ystyle=1,bgblack=bgblack
     loadct,13,/silent
     oplot,dv,t1h,color=250,linestyle=2     ; tau=1 height
     loadct,0,/silent
     if plotvz then oplot,vz,y,color=250

  
     xyouts,0.95*vrange[0],0.9*zrange[1], '!9t!dn!x!nexp(-!9t!dn!x!n)',$
            color=255,charsize=labelsize
    
; contribution function

;     im=bytscl(cfc)
im=cfc
     mplot_image,im,x,y,$
                position=v_pos[3,*],$
                xtitle='!9Dn!x [km/s]',ytitle=' ',ytickname=empty,$
                xthick=2,ythick=2,/noerase,ystyle=9,$
                 xrange=vrange,yrange=zrange,bgblack=bgblack
     loadct,13,/silent
     oplot,dv,t1h,color=250  ,linestyle=2     ; tau=1 height
     loadct,0,/silent
     if plotvz then oplot,vz,y,color=250
  
     xyouts,0.95*vrange[0],0.9*zrange[1], 'C!dI!n',color=255,charsize=labelsize
     

     plot,dv,ie,position=v_pos[3,*],/noerase,$
       xstyle=13,ystyle=13,/nodata,$
       xrange=vrange,yrange=irange
     oplot,dv,ie,color=255
     axis,vrange[1],0,yaxis=1,ystyle=1,yticklen=-0.02,$
       ytitle='I!d!9n!n!x [kK]'


     ; PLOT TITLE OF REQUESTED
     if n_elements(title) gt 0 then begin
        xyouts,0.05,0.95,/norm,title
     endif


     if  keyword_set(onscreen) then begin  
        !p.font=-1
        !p.charsize=1
     endif else begin
        !p.font=-1
        stopplot
     endelse

     if keyword_set(debug) then stop
  
end

FUNCTION tradiation, intensity, lambda0

;+
; NAME:
;	TRADIATION
;
; PURPOSE:
;	This function returns the radiation temperature of corresponding
;       to the intensities INTENSITY at central wavelength LAMBDA0.
;
; CATEGORY:
;	Radiative transfer
;
; CALLING SEQUENCE:
;	Result = TRADIATION(intensity, lambda0)
;
; INPUTS:
;	Intensity:  intensity [J m^-2 s^-1 Hz^-1 sr^-1]
;       Lambda0:    central wavelength [nm]
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	Returns radiation temperature in K.
;
; OPTIONAL OUTPUTS:
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
; 	Written by:    Han Uitenbroek
;
;   --- Last modified: Thu Apr  9 12:33:46 1998 --
;-

  IF (n_params(0) LT 2) THEN BEGIN
    print, "Usage: Trad = TRADIATION(intensity, lambda0)"
    return, 0.0
  ENDIF

  CLIGHT     = 2.99792458E+08
  HPLANCK    = 6.626176E-34
  KBOLTZMANN = 1.380662E-23
  NM_TO_M    = 1.0E-09

  IF (n_elements(lambda0) GT 1) THEN BEGIN
    lambda = avg(lambda0)
    print, FORMAT='("Using average wavelength: ", F8.2, " [nm])', lambda
    lambda = lambda * NM_TO_M
  ENDIF ELSE $ 
   lambda = lambda0 * NM_TO_M

  trad = (HPLANCK*CLIGHT) / ((lambda*KBOLTZMANN) * $
             alog(1.0 + (2.0*HPLANCK*CLIGHT)/(lambda^3 * intensity)))

  return, trad
END


pro calc_position,pos,cpos,xs,ys,nx=nx,ny=ny,$
  xplotsz=xplotsz,xleft=xleft,xright=xright,xdist=xdist,$
  yplotsz=yplotsz,ytop=ytop,ybottom=ybottom,ydist=ydist,cbv=cbv
 
; calculates plot positions in normalized coordinates for multiple
; plots on a page.  outputs pos, an [nx * ny,4] float array, with the
; normalized plot position. The first index goes from left to right
; and from top to bottom.  if cpos keyword is present, it also outputs
; cpos, an[nx * ny,4] float array, containing the normalized
; coordinates of a colorbar position to the left of the image, 4 mm
; from the image to the left. make sure that xleft is larger or equal
; to xdist and that there is enough space (typical 4 cm to put in the
; bar)
; set cbv to get colorbar above panels
 
  on_error,2
  
  if (n_elements(nx) eq 0) then nx=1.
  if (n_elements(ny) eq 0) then ny=1.
  if (n_elements(xplotsz) eq 0) then xplotsz=5.
  if (n_elements(yplotsz) eq 0) then yplotsz=5.
  if (n_elements(xleft) eq 0) then xleft=2.
  if (n_elements(xright) eq 0) then xright=1.
  if (n_elements(ybottom) eq 0) then ybottom=2.
  if (n_elements(ytop) eq 0) then ytop=1.
  if (n_elements(xdist) eq 0) then xdist=0.
  if (n_elements(ydist) eq 0) then ydist=0.

  xs=xleft+xright+nx*xplotsz+(nx-1)*xdist
  ys=ytop+ybottom+ny*yplotsz+(ny-1)*ydist
  pos=fltarr(nx*ny,4)

  for x=0,nx-1 do begin
      for y=0,ny-1 do begin
          x1=(xleft+x*xplotsz+x*xdist)/xs
          x2=(xleft+(x+1)*xplotsz+x*xdist)/xs
          y1=(ybottom+(ny-y-1)*yplotsz+(ny-y-1)*ydist)/ys
          y2=(ybottom+(ny-y)*yplotsz+(ny-y-1)*ydist)/ys          
          pos[x+nx*y,*]=[x1,y1,x2,y2]  
      endfor 
  endfor

  if not keyword_set(cbv) then begin
      cpos=fltarr(nx*ny,4)
      for x=0,nx-1 do begin
          for y=0,ny-1 do begin
              x1=(xleft+(x+1)*xplotsz+x*xdist+0.3)/xs
              x2=(xleft+(x+1)*xplotsz+x*xdist+0.7)/xs
              y1=(ybottom+(ny-y-1)*yplotsz+(ny-y-1)*ydist)/ys
              y2=(ybottom+(ny-y)*yplotsz+(ny-y-1)*ydist)/ys          
              cpos[x+nx*y,*]=[x1,y1,x2,y2]  
          endfor 
      endfor 
  endif else begin
      cpos=fltarr(nx*ny,4)
        for x=0,nx-1 do begin
          for y=0,ny-1 do begin
              x1=(xleft+(x)*xplotsz+x*xdist)/xs
              x2=(xleft+(x+1)*xplotsz+x*xdist)/xs
              y1=(ybottom+(ny-y)*yplotsz+(ny-y-1)*ydist+0.1)/ys
              y2=(ybottom+(ny-y)*yplotsz+(ny-y-1)*ydist+0.3)/ys          
              cpos[x+nx*y,*]=[x1,y1,x2,y2]  
          endfor 
      endfor 
  endelse

end


pro vertint,x,z,c,c0=c0

if n_params() lt 3 then begin
message,/info,'pro vertint,x,z,c,c0=c0'
return
endif

sz=size(x)
if keyword_set(c0) then cn=c0 else cn=1e-30

if sz[0] eq 1 then begin
   nz=sz[1]
   c=fltarr(nz)
   c[0]=cn
   for k=1,nz-1 do c[k]=c[k-1] + 0.5*(x[k]+x[k-1])* abs(z[k]-z[k-1])
endif

if sz[0] eq 2 then begin
   nx=sz[1]
   nz=sz[2]
   c=fltarr(nx,nz)
   c[*,0]=cn
   for k=1,nz-1 do c[*,k]=c[*,k-1] + 0.5*(x[*,k]+x[*,k-1])* abs(z[k]-z[k-1])
 endif

if sz[0] eq 3 then begin
   nx=sz[1]
   ny=sz[2]
   nz=sz[3]
   c=fltarr(nx,ny,nz)
   c[*,*,0]=cn
   for k=1,nz-1 do c[*,*,k]=c[*,*,k-1] + 0.5*(x[*,*,k]+x[*,*,k-1])* abs(z[k]-z[k-1])
endif

end
