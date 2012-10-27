c vim: set filetype=fortran et ts=2 sw=2 sts=2 :
c       Heavily Modified a from Phil Fischer's original code by 
c       Now has the right corrections for sub-pixel integration
c           Erin Sheldon
c
c       whyflag values:
c       0: abs(e1-e1old) .lt. tol1 .and. abs(e2-e2old).lt.tol1
c            .and. abs(m(1,1)/m11old-1.).le.tol2
c       2**0: sum  <= 0
c       2**1: abs(xcen-xcenorig) > shiftmax or abs(ycen-ycenorig) > shiftmax
c       2**2: another sum <= 0, but what's the diff?
c       2**3: m(1,1) <= 0 and m(2,2) <= 0
c       2**4: detm.le.1.e-7
c       2**5: detn.le.0
c       2**6: w(1,1).lt.0..or.w(2,2).lt.0.
c       2**7: imom.eq.maxit
c       2**8: detw <= 0
c       2**9: nsub is not a postive integer


      subroutine ad_mom(image,nx,ny,sky,sigsky,ax,ay,nel,shiftmax,
     &nsub,
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,s2n,numiter,whyflag)

c     calculates adaptive moments
c     applies pixelization corrections

      implicit none

      
      integer nel,nx,ny
      real*8 image(nx,ny)
      real*8 ax(nel),ay(nel)
      real*8 ixx(nel),iyy(nel),ixy(nel)
      real*8 sky(nel),uncer(nel),s2n(nel),sigsky(nel),rho4(nel)
      real*8 shiftmax
      real*8 wcenx(nel),wceny(nel)
      integer numiter(nel),whyflag(nel)

      real*8 x,y,xl,yl,xx,xx2,yy,yy2,td,e1,e2
      real*8 tol1,tol2
      integer maxit

      integer*4 nsub
      real*8 stepsize, offset

      integer imom
      integer ix1,ix2,iy1,iy2,i,j,ii,jj
      real*8 w(2,2),m(2,2),n(2,2)
      real*8 xcen,ycen
      real*8 sumx,sumy,grad,expon,weight,detm,detw,detn,sums4
      real*8 spi,sumxx,sumyy,sumxy,sum,wsum2tot,w1,w2,w12

      integer kk
      real*8 e1old,e2old,m11old
      real*8 xcenorig,ycenorig

      real*8 ymod

      real*8 wsum, w2sum, wwsumx, wwsumy
      real*8 wwsumxx, wwsumxy, wwsumyy, wwexpon2sum
      real*8 weight2
  
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
      parameter(tol1=0.001,tol2=0.01,maxit=100)
      

      if (nsub <= 0) then
        do kk=1,nel
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                uncer(kk),s2n(kk))
          whyflag(kk)=2**9
        enddo
        return
      endif

c     for sub-pixel corrections
      stepsize = 1./nsub
      offset = (nsub-1)*stepsize/2.

c     The main do loop now

      spi=sqrt(atan(1.)*4.)
      do kk=1,nel

        w(1,1)=max(1.,ixx(kk))
        w(2,2)=max(1.,iyy(kk))
        w(1,2)=ixy(kk)

        e1old=10.
        e2old=10.
        m11old=1.e6
        xcen=ax(kk)
        ycen=ay(kk)

        xcenorig=xcen
        ycenorig=ycen
        
        imom=0
        do while(imom.lt.maxit)
          imom=imom+1

c         4-sigma region around object, but within image
          grad=4.*sqrt(max(w(1,1),w(2,2)))
          ix1=nint(max(xcen-grad-0.5,1.))
          iy1=nint(max(ycen-grad-0.5,1.))
          ix2=nint(min(xcen+grad+0.5,float(nx)))
          iy2=nint(min(ycen+grad+0.5,float(ny)))

          sumxx=0.
          sumyy=0.
          sumxy=0.
          sumx=0
          sumy=0.
          sum=0.
          sums4=0.

          detw=w(1,1)*w(2,2)-w(1,2)*w(1,2)

          if(detw.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 9090
          endif

          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw

          ! first get the weighted centroid
          do i=ix1,ix2
            x=i-xcen
            xl=x-offset
            do j=iy1,iy2
              y=j-ycen
              yl=y-offset

              ! work over a 4x4 sub pixel grid
              ! to compute pixel corrections
              wsum=0
              w2sum=0
              wwsumx=0
              wwsumy=0

              xx=xl
              do ii=1,nsub
                xx2=xx*xx
                yy=yl
                do jj=1,nsub
                  yy2=yy*yy

                  expon= xx2*w2 + yy2*w1 - 2.*xx*yy*w12

                  ! the fortran I'm using is actually ok with very
                  ! large exponents
                  weight=exp(-0.5*expon)
                  weight2 = weight*weight

                  wsum=wsum+weight
                  w2sum = w2sum+weight2
                  wwsumx=wwsumx + weight2*(xx+xcen)
                  wwsumy=wwsumy + weight2*(yy+ycen)


                  yy = yy + stepsize
                enddo ! loop y sub pixels
                xx = xx + stepsize
              enddo ! loop x sub pixels

              if (wsum .gt. 0) then
                ymod = image(i,j)-sky(kk)
                sumx = sumx + ymod*wwsumx/wsum
                sumy = sumy + ymod*wwsumy/wsum
                sum = sum + ymod*w2sum/wsum
              endif

            enddo
          enddo

          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**0
            goto 9090
          endif

          xcen=sumx/sum
          ycen=sumy/sum

          if(abs(xcen-xcenorig).gt.shiftmax.or.
     &    abs(ycen-ycenorig).gt.shiftmax)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**1
            goto 9090
          endif
          
          ! now with the new centroid, measure the weighted moments
          ! with sub-pixel corrections
          sum=0.
          wsum2tot=0.
          do i=ix1,ix2
            x=i-xcen
            xl=x-offset
            do j=iy1,iy2
              y=j-ycen
              yl=y-offset

              ! derive the correction factors from the subpixel
              ! grid
              wsum=0
              w2sum=0
              wwsumxx=0
              wwsumxy=0
              wwsumyy=0
              wwexpon2sum=0

              xx=xl
              do ii=1,nsub
                xx2=xx*xx
                yy=yl
                do jj=1,nsub
                  yy2=yy*yy

                  expon= xx2*w2 + yy2*w1 - 2.*xx*yy*w12

                  ! the fortran I'm using is actually ok with very
                  ! large exponents
                  weight=exp(-0.5*expon)
                  weight2 = weight*weight

                  wsum = wsum+weight
                  w2sum = w2sum+weight2
                  wwsumxx = wwsumxx + weight2*xx2
                  wwsumxy = wwsumxy + weight2*xx*yy
                  wwsumyy = wwsumyy + weight2*yy2
                  wwexpon2sum = wwexpon2sum + weight2*expon*expon

                  yy = yy + stepsize
                enddo ! loop y sub pixels
                xx = xx + stepsize
              enddo ! loop x sub pixels

              if (wsum .gt. 0) then
                ymod = image(i,j)-sky(kk)
                sumxx = sumxx + ymod*wwsumxx/wsum
                sumxy = sumxy + ymod*wwsumxy/wsum
                sumyy = sumyy + ymod*wwsumyy/wsum
                sums4 = sums4 + ymod*wwexpon2sum/wsum
                sum = sum + ymod*w2sum/wsum
                wsum2tot = wsum2tot + (w2sum/wsum)**2
              endif

            enddo
          enddo

          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**2
            goto 9090
          endif
          
          m(1,1)=sumxx/sum
          m(2,2)=sumyy/sum
          m(1,2)=sumxy/sum
          if(m(1,1).le.0..and.m(2,2).le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**3
            goto 9090
          endif

          td=w(1,1)+w(2,2)
          e1=(w(1,1)-w(2,2))/td
          e2=2.*w(1,2)/td

          !print *,imom,sum,td,e1,e2

          if(abs(e1-e1old).lt.tol1.and.abs(e2-e2old).lt.tol1.and.
     &    abs(m(1,1)/m11old-1.).le.tol2)then

            ! convergence criteria met
            ixx(kk)=w(1,1)
            iyy(kk)=w(2,2)
            ixy(kk)=w(1,2)
            rho4(kk)=sums4/sum
            detw=((w(1,1)*w(2,2)-w(1,2)*w(1,2)))**0.25
            whyflag(1)=0
            if(4.*sum-sums4.gt.0.)then
              uncer(kk)=4.*spi*sigsky(kk)*detw/(4.*sum-sums4)
            else
              uncer(kk)=9999.
            endif
            if (wsum2tot .gt. 0.) then
              s2n(kk)=sum/sqrt(wsum2tot)/sigsky(kk)
            else
              s2n(kk)=-9999.
            endif
            wcenx(kk)=xcen
            wceny(kk)=ycen
            numiter(kk)=imom
            goto 9090

          else

            detm=(m(1,1)*m(2,2)-m(1,2)*m(1,2))
            if(detm.le.1.e-7)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                    uncer(kk),s2n(kk))
              numiter(kk)=imom
              whyflag(kk)=2**4
              goto 9090
            endif
            
            ! no convergence, set a new weight function from the
            ! difference of the measured covar and the weight
            ! covar, inverted
            detm=1./detm
            detw=1./detw
            n(1,1)=m(2,2)*detm-w(2,2)*detw
            n(2,2)=m(1,1)*detm-w(1,1)*detw
            n(1,2)=-m(1,2)*detm+w(1,2)*detw
            detn=n(1,1)*n(2,2)-n(1,2)*n(1,2)
            
            if(detn.le.0.)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                    uncer(kk),s2n(kk))
              numiter(kk)=imom
              whyflag(kk)=2**5
              goto 9090
            endif
            
            detn=1./detn
            w(1,1)=n(2,2)*detn
            w(2,2)=n(1,1)*detn
            w(1,2)=-n(1,2)*detn
            e1old=e1
            e2old=e2
            m11old=m(1,1)

          endif

          if(w(1,1).lt.0..or.w(2,2).lt.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**6
            goto 9090
          endif

        enddo
        
        if(imom.eq.maxit)then
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                uncer(kk),s2n(kk))
          numiter(kk)=imom
          whyflag(kk)=2**7
          goto 9090
        endif
        
 9090 continue
      enddo

      return

      end




      subroutine setbad(ixx,iyy,ixy,rho4,uncer,s2n)
      
      real*8 ixx,iyy,ixy,rho4,uncer,s2n

      ixx=-9999.
      iyy=-9999.
      ixy=-9999.
      rho4=-9999.
      uncer=9999.
      s2n=-9999.

      return

      end




c add moments from a psf to the moments we fit, in order
c to account for PSF effects
      subroutine ad_mom_1psf(image,nx,ny,sky,sigsky,ax,ay,nel,shiftmax,
     &ixx_psf1,ixy_psf1,iyy_psf1,
     &nsub,
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,s2n,numiter,whyflag)

c     calculates adaptive moments
c     applies pixelization corrections

      implicit none

      
      integer nel,nx,ny
      real*8 image(nx,ny)
      real*8 ax(nel),ay(nel)
      real*8 ixx(nel),iyy(nel),ixy(nel)
      real*8 ixx_psf1(nel),iyy_psf1(nel),ixy_psf1(nel)
      real*8 sky(nel),uncer(nel),s2n(nel),sigsky(nel),rho4(nel)
      real*8 shiftmax
      real*8 wcenx(nel),wceny(nel)
      integer numiter(nel),whyflag(nel)

      real*8 x,y,xl,yl,xx,xx2,yy,yy2,td,e1,e2
      real*8 tol1,tol2
      integer maxit

      integer*4 nsub
      real*8 stepsize, offset

      integer imom
      integer ix1,ix2,iy1,iy2,i,j,ii,jj
      real*8 w(2,2),m(2,2),n(2,2)
      real*8 xcen,ycen
      real*8 sumx,sumy,grad,expon,weight,detm,detw,detn,sums4
      real*8 spi,sumxx,sumyy,sumxy,sum,wsumtot,w1,w2,w12

      integer kk
      real*8 e1old,e2old,m11old
      real*8 xcenorig,ycenorig

      real*8 ymod

      real*8 wsum, w2sum, wwsumx, wwsumy
      real*8 wwsumxx, wwsumxy, wwsumyy, wwexpon2sum
      real*8 weight2
  
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
      parameter(tol1=0.001,tol2=0.01,maxit=100)
      

      if (nsub <= 0) then
        do kk=1,nel
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                uncer(kk),s2n(kk))
          whyflag(kk)=2**9
        enddo
        return
      endif

c     for sub-pixel corrections
      stepsize = 1./nsub
      offset = (nsub-1)*stepsize/2.

c     The main do loop now

      spi=sqrt(atan(1.)*4.)
      do kk=1,nel

        w(1,1)=max(1.,ixx(kk)+ixx_psf1(kk))
        w(2,2)=max(1.,iyy(kk)+iyy_psf1(kk))
        w(1,2)=ixy(kk)+ixy_psf1(kk)

        e1old=10.
        e2old=10.
        m11old=1.e6
        xcen=ax(kk)
        ycen=ay(kk)

        xcenorig=xcen
        ycenorig=ycen
        
        imom=0
        do while(imom.lt.maxit)
          imom=imom+1

c         4-sigma region around object, but within image
          grad=4.*sqrt(max(w(1,1),w(2,2)))
          ix1=nint(max(xcen-grad-0.5,1.))
          iy1=nint(max(ycen-grad-0.5,1.))
          ix2=nint(min(xcen+grad+0.5,float(nx)))
          iy2=nint(min(ycen+grad+0.5,float(ny)))

          sumxx=0.
          sumyy=0.
          sumxy=0.
          sumx=0
          sumy=0.
          sum=0.
          sums4=0.

          detw=w(1,1)*w(2,2)-w(1,2)*w(1,2)

          if(detw.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 9090
          endif

          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw

          ! first get the weighted centroid
          do i=ix1,ix2
            x=i-xcen
            xl=x-offset
            do j=iy1,iy2
              y=j-ycen
              yl=y-offset

              ! work over a 4x4 sub pixel grid
              ! to compute pixel corrections
              wsum=0
              w2sum=0
              wwsumx=0
              wwsumy=0

              xx=xl
              do ii=1,nsub
                xx2=xx*xx
                yy=yl
                do jj=1,nsub
                  yy2=yy*yy

                  expon= xx2*w2 + yy2*w1 - 2.*xx*yy*w12

                  ! the fortran I'm using is actually ok with very
                  ! large exponents
                  weight=exp(-0.5*expon)
                  weight2 = weight*weight

                  wsum=wsum+weight
                  w2sum = w2sum+weight2
                  wwsumx=wwsumx + weight2*(xx+xcen)
                  wwsumy=wwsumy + weight2*(yy+ycen)


                  yy = yy + stepsize
                enddo ! loop y sub pixels
                xx = xx + stepsize
              enddo ! loop x sub pixels

              if (wsum .gt. 0) then
                ymod = image(i,j)-sky(kk)
                sumx = sumx + ymod*wwsumx/wsum
                sumy = sumy + ymod*wwsumy/wsum
                sum = sum + ymod*w2sum/wsum
              endif

            enddo
          enddo

          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**0
            goto 9090
          endif

          xcen=sumx/sum
          ycen=sumy/sum

          if(abs(xcen-xcenorig).gt.shiftmax.or.
     &    abs(ycen-ycenorig).gt.shiftmax)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**1
            goto 9090
          endif
          
          ! now with the new centroid, measure the weighted moments
          ! with sub-pixel corrections
          sum=0.
          wsumtot=0.
          do i=ix1,ix2
            x=i-xcen
            xl=x-offset
            do j=iy1,iy2
              y=j-ycen
              yl=y-offset

              ! derive the correction factors from the subpixel
              ! grid
              wsum=0
              w2sum=0
              wwsumxx=0
              wwsumxy=0
              wwsumyy=0
              wwexpon2sum=0

              xx=xl
              do ii=1,nsub
                xx2=xx*xx
                yy=yl
                do jj=1,nsub
                  yy2=yy*yy

                  expon= xx2*w2 + yy2*w1 - 2.*xx*yy*w12

                  ! the fortran I'm using is actually ok with very
                  ! large exponents
                  weight=exp(-0.5*expon)
                  weight2 = weight*weight

                  wsum = wsum+weight
                  w2sum = w2sum+weight2
                  wwsumxx = wwsumxx + weight2*xx2
                  wwsumxy = wwsumxy + weight2*xx*yy
                  wwsumyy = wwsumyy + weight2*yy2
                  wwexpon2sum = wwexpon2sum + weight2*expon*expon

                  yy = yy + stepsize
                enddo ! loop y sub pixels
                xx = xx + stepsize
              enddo ! loop x sub pixels

              if (wsum .gt. 0) then
                ymod = image(i,j)-sky(kk)
                sumxx = sumxx + ymod*wwsumxx/wsum
                sumxy = sumxy + ymod*wwsumxy/wsum
                sumyy = sumyy + ymod*wwsumyy/wsum
                sums4 = sums4 + ymod*wwexpon2sum/wsum
                sum = sum + ymod*w2sum/wsum
                wsumtot = wsumtot + w2sum/wsum
              endif

            enddo
          enddo

          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**2
            goto 9090
          endif
          
          m(1,1)=sumxx/sum
          m(2,2)=sumyy/sum
          m(1,2)=sumxy/sum
          if(m(1,1).le.0..and.m(2,2).le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**3
            goto 9090
          endif

          td=w(1,1)+w(2,2)
          e1=(w(1,1)-w(2,2))/td
          e2=2.*w(1,2)/td

          !print *,imom,sum,td,e1,e2

          if(abs(e1-e1old).lt.tol1.and.abs(e2-e2old).lt.tol1.and.
     &    abs(m(1,1)/m11old-1.).le.tol2)then

            ! convergence criteria met
            ixx(kk)=w(1,1)-ixx_psf1(kk)
            iyy(kk)=w(2,2)-iyy_psf1(kk)
            ixy(kk)=w(1,2)-ixy_psf1(kk)
            rho4(kk)=sums4/sum
            detw=((w(1,1)*w(2,2)-w(1,2)*w(1,2)))**0.25
            whyflag(1)=0
            if(4.*sum-sums4.gt.0.)then
              uncer(kk)=4.*spi*sigsky(kk)*detw/(4.*sum-sums4)
            else
              uncer(kk)=9999.
            endif
            if (wsumtot .gt. 0.) then
              s2n(kk)=sum/sqrt(wsumtot)/sigsky(kk)
            else
              s2n(kk)=-9999.
            endif
            wcenx(kk)=xcen
            wceny(kk)=ycen
            numiter(kk)=imom
            goto 9090

          else

            detm=(m(1,1)*m(2,2)-m(1,2)*m(1,2))
            if(detm.le.1.e-7)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                    uncer(kk),s2n(kk))
              numiter(kk)=imom
              whyflag(kk)=2**4
              goto 9090
            endif
            
            ! no convergence, set a new weight function from the
            ! difference of the measured covar and the weight
            ! covar, inverted
            detm=1./detm
            detw=1./detw
            n(1,1)=m(2,2)*detm-w(2,2)*detw
            n(2,2)=m(1,1)*detm-w(1,1)*detw
            n(1,2)=-m(1,2)*detm+w(1,2)*detw
            detn=n(1,1)*n(2,2)-n(1,2)*n(1,2)
            
            if(detn.le.0.)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                    uncer(kk),s2n(kk))
              numiter(kk)=imom
              whyflag(kk)=2**5
              goto 9090
            endif
            
            detn=1./detn
            w(1,1)=n(2,2)*detn
            w(2,2)=n(1,1)*detn
            w(1,2)=-n(1,2)*detn
            e1old=e1
            e2old=e2
            m11old=m(1,1)

          endif

          if(w(1,1).lt.0..or.w(2,2).lt.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                  uncer(kk),s2n(kk))
            numiter(kk)=imom
            whyflag(kk)=2**6
            goto 9090
          endif

        enddo
        
        if(imom.eq.maxit)then
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),
     &                uncer(kk),s2n(kk))
          numiter(kk)=imom
          whyflag(kk)=2**7
          goto 9090
        endif
        
 9090 continue
      enddo

      return

      end



