c vim: set filetype=fortran et ts=2 sw=2 sts=2 :
c       Heavily Modified a little from Phil Fischer's original code by 
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
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,numiter,whyflag)

c     calculates adaptive moments
c     applies pixelization corrections

      implicit none

      
      integer nel,nx,ny
      real image(nx,ny)
      real ax(nel),ay(nel)
      real ixx(nel),iyy(nel),ixy(nel)
      real sky(nel),uncer(nel),sigsky(nel),rho4(nel)
      real shiftmax
      real wcenx(nel),wceny(nel)
      integer numiter(nel),whyflag(nel)

      real x,y,xl,yl,xx,xx2,yy,yy2,td,e1,e2
      real tol1,tol2
      integer maxit

      integer*4 nsub
      real*4 stepsize, offset

      integer imom
      integer ix1,ix2,iy1,iy2,i,j,ii,jj
      real w(2,2),m(2,2),n(2,2)
      real xcen,ycen
      real sumx,sumy,grad,expon,weight,detm,detw,detn,sums4
      real spi,sumxx,sumyy,sumxy,sum,w1,w2,w12

      integer kk
      real e1old,e2old,m11old
      real xcenorig,ycenorig

      real ymod

      real wsum, w2sum, wwsumx, wwsumy
      real wwsumxx, wwsumxy, wwsumyy, wwexpon2sum
      real weight2
  
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
      parameter(tol1=0.001,tol2=0.01,maxit=100)
      

      if (nsub <= 0) then
        do kk=1,nel
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
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
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 8080
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
                  ! large exponents, but just to be safe...
                  if (expon .le. 100) then
                    weight=exp(-0.5*expon)
                    weight2 = weight*weight

                    wsum=wsum+weight
                    w2sum = w2sum+weight2
                    wwsumx=wwsumx + weight2*(xx+xcen)
                    wwsumy=wwsumy + weight2*(yy+ycen)

                  endif

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
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**0
            goto 8080
          endif

          xcen=sumx/sum
          ycen=sumy/sum

          if(abs(xcen-xcenorig).gt.shiftmax.or.
     &    abs(ycen-ycenorig).gt.shiftmax)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**1
            goto 8080
          endif
          
          ! now with the new centroid, measure the weighted moments
          ! with sub-pixel corrections
          sum=0.
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
                  ! large exponents, but just to be safe...
                  if (expon .le. 100) then
                    weight=exp(-0.5*expon)
                    weight2 = weight*weight

                    wsum = wsum+weight
                    w2sum = w2sum+weight2
                    wwsumxx = wwsumxx + weight2*xx2
                    wwsumxy = wwsumxy + weight2*xx*yy
                    wwsumyy = wwsumyy + weight2*yy2
                    wwexpon2sum = wwexpon2sum + weight2*expon*expon

                  endif

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
              endif

            enddo
          enddo

          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**2
            goto 8080
          endif
          
          m(1,1)=sumxx/sum
          m(2,2)=sumyy/sum
          m(1,2)=sumxy/sum
          if(m(1,1).le.0..and.m(2,2).le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**3
            goto 8080
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
            wcenx(kk)=xcen
            wceny(kk)=ycen
            numiter(kk)=imom
            goto 8080

          else

            detm=(m(1,1)*m(2,2)-m(1,2)*m(1,2))
            if(detm.le.1.e-7)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=2**4
              goto 8080
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
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=2**5
              goto 8080
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
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**6
            goto 8080
          endif

        enddo
        
        if(imom.eq.maxit)then
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
          numiter(kk)=imom
          whyflag(kk)=2**7
          goto 8080
        endif
        
 8080 continue
      enddo

      return

      end




      subroutine setbad(ixx,iyy,ixy,rho4,uncer)
      
      real ixx,iyy,ixy,rho4,uncer

      ixx=-9999.0
      iyy=-9999.0
      ixy=-9999.0
      rho4=-9999.0
      uncer=9999.0

      return

      end





      subroutine ad_mom_nosub(image,nx,ny,sky,sigsky,ax,ay,nel,shiftmax,
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,numiter,whyflag)


      implicit none

      real tol1,tol2
      real td,e1,e2
      integer maxit
      
      integer nel,nx,ny
      real image(nx,ny)
      real ax(nel),ay(nel)
      real ixx(nel),iyy(nel),ixy(nel)
      real sky(nel),uncer(nel),sigsky(nel),shiftmax,rho4(nel)
      real wcenx(nel),wceny(nel)
      integer numiter(nel),whyflag(nel)

      integer imom
      integer ix1,ix2,iy1,iy2,i,j
      real w(2,2),m(2,2),n(2,2)
      real xcen,ycen
      real sumx,sumy,grad,expon,weight,detm,detw,detn,sums4
      real spi,sumxx,sumyy,sumxy,sum,w1,w2,w12,x,x2,y,y2,xy

      integer kk

      real e1old,e2old,m11old
      real xcenorig,ycenorig
      real ymod

  
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
      parameter(tol1=0.001,tol2=0.01,maxit=100)
      

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

          ! 4-sigma region
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
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 6060
          endif


          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw

          do i=ix1,ix2
            x=i-xcen
            do j=iy1,iy2
              y=j-ycen

              expon=x*x*w2 + y*y*w1 -2.*x*y*w12

              if(expon.le.100.)then

                weight=exp(-0.5*expon)
                ymod=(image(i,j)-sky(kk))*weight
                sumx=sumx+ymod*float(i)
                sumy=sumy+ymod*float(j)
                sum=sum+ymod

              endif
            enddo

          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**0
            goto 6060
          endif
          xcen=sumx/sum
          ycen=sumy/sum
          if(abs(xcen-xcenorig).gt.shiftmax.or.
     &    abs(ycen-ycenorig).gt.shiftmax)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**1
            goto 6060
          endif
          
          sum=0.
          do i=ix1,ix2
            x=i-xcen
            x2=x*x
            do j=iy1,iy2
              y=j-ycen
              y2=y*y

              xy=x*y
              expon=x2*w2+y2*w1-2.*xy*w12

              if(expon.le.100.)then

                weight=exp(-0.5*expon)
                ymod=(image(i,j)-sky(kk))*weight
                sumxx=sumxx+x2*ymod
                sumyy=sumyy+y2*ymod
                sumxy=sumxy+xy*ymod
                sums4=sums4+expon*expon*ymod
                sum=sum+ymod

              endif
            enddo

          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**2
            goto 6060
          endif
          
          m(1,1)=sumxx/sum
          m(2,2)=sumyy/sum
          m(1,2)=sumxy/sum
          if(m(1,1).le.0..and.m(2,2).le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**3
            goto 6060
          endif
          td=w(1,1)+w(2,2)
          e1=(w(1,1)-w(2,2))/td
          e2=2.*w(1,2)/td
          if(abs(e1-e1old).lt.tol1.and.abs(e2-e2old).lt.tol1.and.
     &    abs(m(1,1)/m11old-1.).le.tol2)then

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
            wcenx(kk)=xcen
            wceny(kk)=ycen
            numiter(kk)=imom
            goto 6060

          else
            detm=(m(1,1)*m(2,2)-m(1,2)*m(1,2))
            if(detm.le.1.e-7)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=2**4
              goto 6060
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
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=2**5
              goto 6060
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
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**6
            goto 6060
          endif
        enddo
        
        if(imom.eq.maxit)then
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
          numiter(kk)=imom
          whyflag(kk)=2**7
          goto 6060
        endif
        
 6060 continue
      enddo

      return

      end



