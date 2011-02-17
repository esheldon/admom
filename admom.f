c vim: set filetype=fortran et ts=2 sw=2 sts=2 :
c       Modified a little from Phil Fischer's original code by 
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

      subroutine ad_mom(data,nx,ny,sky,sigsky,ax,ay,nel,shiftmax,
     &allowinterp,
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,numiter,interpolated,whyflag)

c     calculates adaptive moments
c     interpolates weights for poorly sampled objects

      implicit none

      real tol1,tol2,xinterp2,xl,xh,yl,yh,xxx,yyy,xy,td,e1,e2,xinterp
      integer maxit

      logical allowinterp
      
      integer nel,nx,ny
      real ax(nel),ay(nel)
      real ixx(nel),iyy(nel),ixy(nel)
      real sky(nel),uncer(nel),sigsky(nel),shiftmax,rho4(nel)
      integer imom
      integer ix1,ix2,iy1,iy2,i,j,ii,jj
      real w(2,2),m(2,2),n(2,2)
      real xcen,ycen
      real sumx,sumy,grad,expon,weight,detm,detw,detn
      real xcenorig,ycenorig,sums4
      real data(nx,ny),spi,sumxx,sumyy,sumxy,sum,w1,w2,w12,xx,xx2,yy,yy2
      integer numiter(nel),whyflag(nel),kk
      real wcenx(nel),wceny(nel),e1old,e2old,m11old
      real ymod,tymod
      integer*2 interpolated(nel)
      logical interpflag
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
c     This xinterp converts to sigma as sqrt(xinterp).  3 was too small
c     it turns out.  I got 5.3 from monte-carlos, which is
c     sigma=2.3 pixels
      parameter(tol1=0.001,tol2=0.01,maxit=100,xinterp=5.3,
     &xinterp2=xinterp*xinterp)

      
c     The main do loop now


      spi=sqrt(atan(1.)*4.)
      do kk=1,nel
        w(1,1)=max(1.,ixx(kk))
        w(2,2)=max(1.,iyy(kk))
        w(1,2)=ixy(kk)
c     the starting moments of the weight
        e1old=10.
        e2old=10.
        m11old=1.e6
        xcen=ax(kk)
        ycen=ay(kk)
        xcenorig=xcen
        ycenorig=ycen
c     the starting centroids
        
        imom=0
        do while(imom.lt.maxit)
          interpolated(kk)=0
          interpflag=.false.
          imom=imom+1
c     the number of iterations
          grad=4.*sqrt(max(w(1,1),w(2,2)))
          ix1=nint(max(xcen-grad-0.5,1.))
          iy1=nint(max(ycen-grad-0.5,1.))
          ix2=nint(min(xcen+grad+0.5,float(nx)))
          iy2=nint(min(ycen+grad+0.5,float(ny)))
c     defined the neighborhood of the object
          sumxx=0.
          sumyy=0.
          sumxy=0.
          sumx=0
          sumy=0.
          sum=0.
          sums4=0.
c          rmax=0.
c     loop over the neighborhood pixels
          detw=w(1,1)*w(2,2)-w(1,2)*w(1,2)
          if(detw.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 9090
          endif

          if (allowinterp) then
            if(w(1,1).lt.xinterp.or.w(2,2).lt.xinterp.or.
     &          detw.lt.xinterp2) then
              interpflag=.true.
              interpolated(kk)=1
            endif
          endif
          !interpflag=.false.
          !interpolated(kk)=0

          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.375
            xh=xx+0.375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.375
                yh=yy+0.375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then

                  tymod = data(i,j)-sky(kk)

                  ! work over a 4x4 sub pixel grid
                  xxx=xl
                  do ii=1,4
                    xx2=xxx*xxx
                    yyy=yl
                    do jj=1,4
                      yy2=yyy*yyy

                      expon= xx2*w2 + yy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)

                      ymod=tymod*weight/16.

                      !sumx=sumx+ymod*float(i)
                      !sumy=sumy+ymod*float(j)
                      sumx=sumx+ymod*(xxx+xcen)
                      sumy=sumy+ymod*(yyy+ycen)
                      sum=sum+ymod

                      yyy = yyy + 0.25
                    enddo ! loop y sub pixels
                    xxx = xxx + 0.25
                  enddo ! loop x sub pixels

                endif
              else
                expon=xx2*w2+yy2*w1-2.*xx*yy*w12
                if(expon.le.9..and.expon.gt.0.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumx=sumx+ymod*float(i)
                  sumy=sumy+ymod*float(j)
                  sum=sum+ymod
                endif
              endif
            enddo
          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**0
            goto 9090
          endif
          xcen=sumx/sum
          ycen=sumy/sum
          if(abs(xcen-xcenorig).gt.shiftmax.or.
     &    abs(ycen-ycenorig).gt.shiftmax)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**1
            goto 9090
          endif
          
          sum=0.
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.375
            xh=xx+0.375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.375
                yh=yy+0.375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then

                  tymod = data(i,j)-sky(kk)

                  ! work over a 4x4 sub pixel grid
                  xxx=xl
                  do ii=1,4
                    xx2=xxx*xxx
                    yyy=yl
                    do jj=1,4
                      yy2=yyy*yyy

                      expon= xx2*w2 + yy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)

                      ymod=tymod*weight/16.

                      sumxx=sumxx+xx2*ymod
                      sumyy=sumyy+yy2*ymod
                      sumxy=sumxy+xxx*yyy*ymod
                      sums4=sums4+expon*expon*ymod
                      sum=sum+ymod

                      yyy = yyy + 0.25
                    enddo ! loop y sub pixels
                    xxx = xxx + 0.25
                  enddo ! loop x sub pixels

                endif
              else
                xy=xx*yy
                expon=xx2*w2+yy2*w1-2.*xy*w12
                if(expon.le.9..and.expon.gt.0.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumxx=sumxx+xx2*ymod
                  sumyy=sumyy+yy2*ymod
                  sumxy=sumxy+xy*ymod
                  sums4=sums4+expon*expon*ymod
                  sum=sum+ymod
                endif
              endif
            enddo
          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**2
            goto 9090
          endif
          
          m(1,1)=sumxx/sum
          m(2,2)=sumyy/sum
          m(1,2)=sumxy/sum
          if(m(1,1).le.0..and.m(2,2).le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**3
            goto 9090
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
c           Don't reset the x,y inputs
c            ax(kk)=xcen
c            ay(kk)=ycen
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
            goto 9090
          else
            detm=(m(1,1)*m(2,2)-m(1,2)*m(1,2))
            if(detm.le.1.e-7)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=2**4
              goto 9090
            endif
            
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
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**6
            goto 9090
          endif
        enddo
        
        if(imom.eq.maxit)then
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
          numiter(kk)=imom
          whyflag(kk)=2**7
          goto 9090
        endif
        

 9090 enddo

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




      subroutine ad_mom_sub4(data,nx,ny,sky,sigsky,ax,ay,nel,shiftmax,
     &allowinterp,
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,numiter,interpolated,whyflag)

c     calculates adaptive moments
c     interpolates weights for poorly sampled objects

      implicit none

      real tol1,tol2,xinterp2,xl,xh,yl,yh
      real xxx,xxx2,yyy,yyy2,xy,td,e1,e2,xinterp
      integer maxit

      logical allowinterp
      
      integer nel,nx,ny
      real ax(nel),ay(nel)
      real ixx(nel),iyy(nel),ixy(nel)
      real sky(nel),uncer(nel),sigsky(nel),shiftmax,rho4(nel)
      integer imom
      integer ix1,ix2,iy1,iy2,i,j,ii,jj
      real w(2,2),m(2,2),n(2,2)
      real xcen,ycen
      real sumx,sumy,grad,expon,weight,detm,detw,detn
      real xcenorig,ycenorig,sums4
      real data(nx,ny),spi,sumxx,sumyy,sumxy,sum,w1,w2,w12,xx,xx2,yy,yy2
      integer numiter(nel),whyflag(nel),kk
      real wcenx(nel),wceny(nel),e1old,e2old,m11old
      real ymod
      integer*2 interpolated(nel)
      logical interpflag

      real wsum, w2sum, wwsumx, wwsumy
      real wwsumxx, wwsumxy, wwsumyy, wwexpon2sum
      real weight2
  
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
c     This xinterp converts to sigma as sqrt(xinterp).  3 was too small
c     it turns out.  I got 5.3 from monte-carlos, which is
c     sigma=2.3 pixels
      parameter(tol1=0.001,tol2=0.01,maxit=100,xinterp=5.3,
     &xinterp2=xinterp*xinterp)

      
c     The main do loop now


      spi=sqrt(atan(1.)*4.)
      do kk=1,nel
        w(1,1)=max(1.,ixx(kk))
        w(2,2)=max(1.,iyy(kk))
        w(1,2)=ixy(kk)
c     the starting moments of the weight
        e1old=10.
        e2old=10.
        m11old=1.e6
        xcen=ax(kk)
        ycen=ay(kk)
        xcenorig=xcen
        ycenorig=ycen
c     the starting centroids
        
        imom=0
        do while(imom.lt.maxit)
          interpolated(kk)=0
          interpflag=.false.
          imom=imom+1
c     the number of iterations
          grad=4.*sqrt(max(w(1,1),w(2,2)))
          ix1=nint(max(xcen-grad-0.5,1.))
          iy1=nint(max(ycen-grad-0.5,1.))
          ix2=nint(min(xcen+grad+0.5,float(nx)))
          iy2=nint(min(ycen+grad+0.5,float(ny)))
c     defined the neighborhood of the object
          sumxx=0.
          sumyy=0.
          sumxy=0.
          sumx=0
          sumy=0.
          sum=0.
          sums4=0.
c          rmax=0.
c     loop over the neighborhood pixels
          detw=w(1,1)*w(2,2)-w(1,2)*w(1,2)
          if(detw.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 8080
          endif

          if (allowinterp) then
            if(w(1,1).lt.xinterp.or.w(2,2).lt.xinterp.or.
     &          detw.lt.xinterp2) then
              interpflag=.true.
              interpolated(kk)=1
            endif
          endif
          !interpflag=.false.
          !interpolated(kk)=0

          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.375
            xh=xx+0.375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.375
                yh=yy+0.375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then


                  ! work over a 4x4 sub pixel grid
                  wsum=0
                  w2sum=0
                  wwsumx=0
                  wwsumy=0

                  xxx=xl
                  do ii=1,4
                    xxx2=xxx*xxx
                    yyy=yl
                    do jj=1,4
                      yyy2=yyy*yyy


                      expon= xxx2*w2 + yyy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)
                      weight2 = weight*weight

                      wsum=wsum+weight
                      w2sum = w2sum+weight2
                      wwsumx=wwsumx + weight2*(xxx+xcen)
                      wwsumy=wwsumy + weight2*(yyy+ycen)

                      yyy = yyy + 0.25
                    enddo ! loop y sub pixels
                    xxx = xxx + 0.25
                  enddo ! loop x sub pixels

                  ymod = data(i,j)-sky(kk)
                  sumx = sumx + ymod*wwsumx/wsum
                  sumy = sumy + ymod*wwsumy/wsum
                  sum = sum + ymod*w2sum/wsum

                endif
              else
                expon=xx2*w2+yy2*w1-2.*xx*yy*w12
                if(expon.le.9..and.expon.gt.0.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumx=sumx+ymod*float(i)
                  sumy=sumy+ymod*float(j)
                  sum=sum+ymod
                endif
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
          
          sum=0.
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.375
            xh=xx+0.375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.375
                yh=yy+0.375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then


                  ! derive the correction factors from the subpixel
                  ! grid
                  wsum=0
                  w2sum=0
                  wwsumxx=0
                  wwsumxy=0
                  wwsumyy=0
                  wwexpon2sum=0

                  xxx=xl
                  do ii=1,4
                    xxx2=xxx*xxx
                    yyy=yl
                    do jj=1,4
                      yyy2=yyy*yyy

                      expon= xxx2*w2 + yyy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)
                      weight2 = weight*weight

                      wsum = wsum+weight
                      w2sum = w2sum+weight2
                      wwsumxx = wwsumxx + weight2*xxx2
                      wwsumxy = wwsumxy + weight2*xxx*yyy
                      wwsumyy = wwsumyy + weight2*yyy2
                      wwexpon2sum = wwexpon2sum + weight2*expon*expon

                      yyy = yyy + 0.25
                    enddo ! loop y sub pixels
                    xxx = xxx + 0.25
                  enddo ! loop x sub pixels

                  ymod = data(i,j)-sky(kk)
                  sumxx = sumxx + ymod*wwsumxx/wsum
                  sumxy = sumxy + ymod*wwsumxy/wsum
                  sumyy = sumyy + ymod*wwsumyy/wsum
                  sums4 = sums4 + ymod*wwexpon2sum/wsum
                  sum = sum + ymod*w2sum/wsum
                  !sum=sum+ymod*wsum/16.


                  ! print comparison
                  !xy=xx*yy
                  !expon=xx2*w2+yy2*w1-2.*xy*w12
                  !weight=exp(-expon)
                  !print *,"(",xx,",",yy,") x2 fac:",(wwsumxx/wsum)

                endif
              else
                xy=xx*yy
                expon=xx2*w2+yy2*w1-2.*xy*w12
                if(expon.le.9..and.expon.gt.0.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumxx=sumxx+xx2*ymod
                  sumyy=sumyy+yy2*ymod
                  sumxy=sumxy+xy*ymod
                  sums4=sums4+expon*expon*ymod
                  sum=sum+ymod
                endif
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
          if(abs(e1-e1old).lt.tol1.and.abs(e2-e2old).lt.tol1.and.
     &    abs(m(1,1)/m11old-1.).le.tol2)then
            ixx(kk)=w(1,1)
            iyy(kk)=w(2,2)
            ixy(kk)=w(1,2)
            rho4(kk)=sums4/sum
c           Don't reset the x,y inputs
c            ax(kk)=xcen
c            ay(kk)=ycen
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
        

 8080 enddo

      return

      end


      subroutine ad_mom_sub8(data,nx,ny,sky,sigsky,ax,ay,nel,shiftmax,
     &allowinterp,
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,numiter,interpolated,whyflag)

c     calculates adaptive moments
c     interpolates weights for poorly sampled objects

      implicit none

      real tol1,tol2,xinterp2,xl,xh,yl,yh
      real xxx,xxx2,yyy,yyy2,xy,td,e1,e2,xinterp
      integer maxit

      logical allowinterp
      
      integer nel,nx,ny
      real ax(nel),ay(nel)
      real ixx(nel),iyy(nel),ixy(nel)
      real sky(nel),uncer(nel),sigsky(nel),shiftmax,rho4(nel)
      integer imom
      integer ix1,ix2,iy1,iy2,i,j,ii,jj
      real w(2,2),m(2,2),n(2,2)
      real xcen,ycen
      real sumx,sumy,grad,expon,weight,detm,detw,detn
      real xcenorig,ycenorig,sums4
      real data(nx,ny),spi,sumxx,sumyy,sumxy,sum,w1,w2,w12,xx,xx2,yy,yy2
      integer numiter(nel),whyflag(nel),kk
      real wcenx(nel),wceny(nel),e1old,e2old,m11old
      real ymod
      integer*2 interpolated(nel)
      logical interpflag

      real wsum, w2sum, wwsumx, wwsumy
      real wwsumxx, wwsumxy, wwsumyy, wwexpon2sum
      real weight2
  
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
c     This xinterp converts to sigma as sqrt(xinterp).  3 was too small
c     it turns out.  I got 5.3 from monte-carlos, which is
c     sigma=2.3 pixels
      parameter(tol1=0.001,tol2=0.01,maxit=100,xinterp=5.3,
     &xinterp2=xinterp*xinterp)

      
c     The main do loop now


      spi=sqrt(atan(1.)*4.)
      do kk=1,nel
        w(1,1)=max(1.,ixx(kk))
        w(2,2)=max(1.,iyy(kk))
        w(1,2)=ixy(kk)
c     the starting moments of the weight
        e1old=10.
        e2old=10.
        m11old=1.e6
        xcen=ax(kk)
        ycen=ay(kk)
        xcenorig=xcen
        ycenorig=ycen
c     the starting centroids
        
        imom=0
        do while(imom.lt.maxit)
          interpolated(kk)=0
          interpflag=.false.
          imom=imom+1
c     the number of iterations
          grad=4.*sqrt(max(w(1,1),w(2,2)))
          ix1=nint(max(xcen-grad-0.5,1.))
          iy1=nint(max(ycen-grad-0.5,1.))
          ix2=nint(min(xcen+grad+0.5,float(nx)))
          iy2=nint(min(ycen+grad+0.5,float(ny)))
c     defined the neighborhood of the object
          sumxx=0.
          sumyy=0.
          sumxy=0.
          sumx=0
          sumy=0.
          sum=0.
          sums4=0.
c          rmax=0.
c     loop over the neighborhood pixels
          detw=w(1,1)*w(2,2)-w(1,2)*w(1,2)
          if(detw.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 7070
          endif

          if (allowinterp) then
            if(w(1,1).lt.xinterp.or.w(2,2).lt.xinterp.or.
     &          detw.lt.xinterp2) then
              interpflag=.true.
              interpolated(kk)=1
            endif
          endif
          !interpflag=.false.
          !interpolated(kk)=0

          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.4375
            xh=xx+0.4375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.4375
                yh=yy+0.4375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then


                  ! work over a 4x4 sub pixel grid
                  wsum=0
                  w2sum=0
                  wwsumx=0
                  wwsumy=0

                  xxx=xl
                  do ii=1,4
                    xxx2=xxx*xxx
                    yyy=yl
                    do jj=1,4
                      yyy2=yyy*yyy


                      expon= xxx2*w2 + yyy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)
                      weight2 = weight*weight

                      wsum=wsum+weight
                      w2sum = w2sum+weight2
                      wwsumx=wwsumx + weight2*(xxx+xcen)
                      wwsumy=wwsumy + weight2*(yyy+ycen)

                      yyy = yyy + 0.125
                    enddo ! loop y sub pixels
                    xxx = xxx + 0.125
                  enddo ! loop x sub pixels

                  ymod = data(i,j)-sky(kk)
                  sumx = sumx + ymod*wwsumx/wsum
                  sumy = sumy + ymod*wwsumy/wsum
                  sum = sum + ymod*w2sum/wsum

                endif
              else
                expon=xx2*w2+yy2*w1-2.*xx*yy*w12
                if(expon.le.9..and.expon.gt.0.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumx=sumx+ymod*float(i)
                  sumy=sumy+ymod*float(j)
                  sum=sum+ymod
                endif
              endif
            enddo
          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**0
            goto 7070
          endif
          xcen=sumx/sum
          ycen=sumy/sum
          if(abs(xcen-xcenorig).gt.shiftmax.or.
     &    abs(ycen-ycenorig).gt.shiftmax)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**1
            goto 7070
          endif
          
          sum=0.
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-0.4375
            xh=xx+0.4375
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-0.4375
                yh=yy+0.4375
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                if(expon.le.9.)then


                  ! derive the correction factors from the subpixel
                  ! grid
                  wsum=0
                  w2sum=0
                  wwsumxx=0
                  wwsumxy=0
                  wwsumyy=0
                  wwexpon2sum=0

                  xxx=xl
                  do ii=1,4
                    xxx2=xxx*xxx
                    yyy=yl
                    do jj=1,4
                      yyy2=yyy*yyy

                      expon= xxx2*w2 + yyy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)
                      weight2 = weight*weight

                      wsum = wsum+weight
                      w2sum = w2sum+weight2
                      wwsumxx = wwsumxx + weight2*xxx2
                      wwsumxy = wwsumxy + weight2*xxx*yyy
                      wwsumyy = wwsumyy + weight2*yyy2
                      wwexpon2sum = wwexpon2sum + weight2*expon*expon

                      yyy = yyy + 0.125
                    enddo ! loop y sub pixels
                    xxx = xxx + 0.125
                  enddo ! loop x sub pixels

                  ymod = data(i,j)-sky(kk)
                  sumxx = sumxx + ymod*wwsumxx/wsum
                  sumxy = sumxy + ymod*wwsumxy/wsum
                  sumyy = sumyy + ymod*wwsumyy/wsum
                  sums4 = sums4 + ymod*wwexpon2sum/wsum
                  sum = sum + ymod*w2sum/wsum
                  !sum=sum+ymod*wsum/16.


                  ! print comparison
                  !xy=xx*yy
                  !expon=xx2*w2+yy2*w1-2.*xy*w12
                  !weight=exp(-expon)
                  !print *,"(",xx,",",yy,") x2 fac:",(wwsumxx/wsum)

                endif
              else
                xy=xx*yy
                expon=xx2*w2+yy2*w1-2.*xy*w12
                if(expon.le.9..and.expon.gt.0.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumxx=sumxx+xx2*ymod
                  sumyy=sumyy+yy2*ymod
                  sumxy=sumxy+xy*ymod
                  sums4=sums4+expon*expon*ymod
                  sum=sum+ymod
                endif
              endif
            enddo
          enddo
          if(sum.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**2
            goto 7070
          endif
          
          m(1,1)=sumxx/sum
          m(2,2)=sumyy/sum
          m(1,2)=sumxy/sum
          if(m(1,1).le.0..and.m(2,2).le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**3
            goto 7070
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
c           Don't reset the x,y inputs
c            ax(kk)=xcen
c            ay(kk)=ycen
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
            goto 7070
          else
            detm=(m(1,1)*m(2,2)-m(1,2)*m(1,2))
            if(detm.le.1.e-7)then
              call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
              numiter(kk)=imom
              whyflag(kk)=2**4
              goto 7070
            endif
            
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
              goto 7070
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
            goto 7070
          endif
        enddo
        
        if(imom.eq.maxit)then
          call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
          numiter(kk)=imom
          whyflag(kk)=2**7
          goto 7070
        endif
        

 7070 enddo

      return

      end


      subroutine ad_mom_sub(data,nx,ny,sky,sigsky,ax,ay,nel,shiftmax,
     &nsub,
     &ixx,ixy,iyy,rho4,wcenx,wceny,uncer,numiter,interpolated,whyflag)

c     calculates adaptive moments
c     interpolates weights for poorly sampled objects

      implicit none

      real tol1,tol2,xinterp2,xl,xh,yl,yh
      real xxx,xxx2,yyy,yyy2,xy,td,e1,e2,xinterp
      integer maxit

      integer*4 nsub
      real*4 stepsize, offset
      
      integer nel,nx,ny
      real ax(nel),ay(nel)
      real ixx(nel),iyy(nel),ixy(nel)
      real sky(nel),uncer(nel),sigsky(nel),shiftmax,rho4(nel)
      integer imom
      integer ix1,ix2,iy1,iy2,i,j,ii,jj
      real w(2,2),m(2,2),n(2,2)
      real xcen,ycen
      real sumx,sumy,grad,expon,weight,detm,detw,detn
      real xcenorig,ycenorig,sums4
      real data(nx,ny),spi,sumxx,sumyy,sumxy,sum,w1,w2,w12,xx,xx2,yy,yy2
      integer numiter(nel),whyflag(nel),kk
      real wcenx(nel),wceny(nel),e1old,e2old,m11old
      real ymod
      integer*2 interpolated(nel)
      logical interpflag

      real wsum, w2sum, wwsumx, wwsumy
      real wwsumxx, wwsumxy, wwsumyy, wwexpon2sum
      real weight2
  
c     I changed tol1 from 0.01 to 0.001 to agree with the C code
c     This xinterp converts to sigma as sqrt(xinterp).  3 was too small
c     it turns out.  I got 5.3 from monte-carlos, which is
c     sigma=2.3 pixels
      parameter(tol1=0.001,tol2=0.01,maxit=100,xinterp=5.3,
     &xinterp2=xinterp*xinterp)

      

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

      print *,"admom:  nsub:",nsub," stepsize:",stepsize,
     & " offset:",offset

c     The main do loop now

      spi=sqrt(atan(1.)*4.)
      do kk=1,nel
        w(1,1)=max(1.,ixx(kk))
        w(2,2)=max(1.,iyy(kk))
        w(1,2)=ixy(kk)
c     the starting moments of the weight
        e1old=10.
        e2old=10.
        m11old=1.e6
        xcen=ax(kk)
        ycen=ay(kk)
        xcenorig=xcen
        ycenorig=ycen
c     the starting centroids
        
        imom=0
        do while(imom.lt.maxit)
          interpolated(kk)=0
          interpflag=.false.
          imom=imom+1
c     the number of iterations
          grad=4.*sqrt(max(w(1,1),w(2,2)))
          ix1=nint(max(xcen-grad-0.5,1.))
          iy1=nint(max(ycen-grad-0.5,1.))
          ix2=nint(min(xcen+grad+0.5,float(nx)))
          iy2=nint(min(ycen+grad+0.5,float(ny)))
c     defined the neighborhood of the object
          sumxx=0.
          sumyy=0.
          sumxy=0.
          sumx=0
          sumy=0.
          sum=0.
          sums4=0.
c          rmax=0.
c     loop over the neighborhood pixels
          detw=w(1,1)*w(2,2)-w(1,2)*w(1,2)
          if(detw.le.0.)then
            call setbad(ixx(kk),iyy(kk),ixy(kk),rho4(kk),uncer(kk))
            numiter(kk)=imom
            whyflag(kk)=2**8
            goto 6060
          endif

          if (nsub > 1) then
            if(w(1,1).lt.xinterp.or.w(2,2).lt.xinterp.or.
     &          detw.lt.xinterp2) then
              interpflag=.true.
              interpolated(kk)=1
            endif
          endif
          !interpflag=.false.
          !interpolated(kk)=0

          w1=w(1,1)/detw
          w2=w(2,2)/detw
          w12=w(1,2)/detw
          do i=ix1,ix2
            xx=i-xcen
            xx2=xx*xx
            xl=xx-offset
            xh=xx+offset
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
              !if(.false.)then
                yl=yy-offset
                yh=yy+offset
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                !if(expon.le.9.)then
                if(expon.le.100.)then


                  ! work over a 4x4 sub pixel grid
                  ! to compute pixel corrections
                  wsum=0
                  w2sum=0
                  wwsumx=0
                  wwsumy=0

                  xxx=xl
                  do ii=1,nsub
                    xxx2=xxx*xxx
                    yyy=yl
                    do jj=1,nsub
                      yyy2=yyy*yyy

                      expon= xxx2*w2 + yyy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)
                      weight2 = weight*weight

                      wsum=wsum+weight
                      w2sum = w2sum+weight2
                      wwsumx=wwsumx + weight2*(xxx+xcen)
                      wwsumy=wwsumy + weight2*(yyy+ycen)

                      yyy = yyy + stepsize
                    enddo ! loop y sub pixels
                    xxx = xxx + stepsize
                  enddo ! loop x sub pixels

                  ymod = data(i,j)-sky(kk)
                  sumx = sumx + ymod*wwsumx/wsum
                  sumy = sumy + ymod*wwsumy/wsum
                  sum = sum + ymod*w2sum/wsum

                endif
              else
                expon=xx2*w2+yy2*w1-2.*xx*yy*w12
                if(expon.le.100.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumx=sumx+ymod*float(i)
                  sumy=sumy+ymod*float(j)
                  sum=sum+ymod
                endif
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
            xx=i-xcen
            xx2=xx*xx
            xl=xx-offset
            xh=xx+offset
            do j=iy1,iy2
              yy=j-ycen
              yy2=yy*yy
              if(interpflag)then
                yl=yy-offset
                yh=yy+offset
                expon=xl*xl*w2+yl*yl*w1-2.*xl*yl*w12
                expon=max(expon,xh*xh*w2+yh*yh*w1-2.*xh*yh*w12)
                expon=max(expon,xl*xl*w2+yh*yh*w1-2.*xl*yh*w12)
                expon=max(expon,xh*xh*w2+yl*yl*w1-2.*xh*yl*w12)
                !if(expon.le.9.)then
                if(expon.le.100.)then


                  ! derive the correction factors from the subpixel
                  ! grid
                  wsum=0
                  w2sum=0
                  wwsumxx=0
                  wwsumxy=0
                  wwsumyy=0
                  wwexpon2sum=0

                  xxx=xl
                  do ii=1,nsub
                    xxx2=xxx*xxx
                    yyy=yl
                    do jj=1,nsub
                      yyy2=yyy*yyy

                      expon= xxx2*w2 + yyy2*w1 - 2.*xxx*yyy*w12
                      weight=exp(-0.5*expon)
                      weight2 = weight*weight

                      wsum = wsum+weight
                      w2sum = w2sum+weight2
                      wwsumxx = wwsumxx + weight2*xxx2
                      wwsumxy = wwsumxy + weight2*xxx*yyy
                      wwsumyy = wwsumyy + weight2*yyy2
                      wwexpon2sum = wwexpon2sum + weight2*expon*expon

                      yyy = yyy + stepsize
                    enddo ! loop y sub pixels
                    xxx = xxx + stepsize
                  enddo ! loop x sub pixels

                  ymod = data(i,j)-sky(kk)
                  sumxx = sumxx + ymod*wwsumxx/wsum
                  sumxy = sumxy + ymod*wwsumxy/wsum
                  sumyy = sumyy + ymod*wwsumyy/wsum
                  sums4 = sums4 + ymod*wwexpon2sum/wsum
                  sum = sum + ymod*w2sum/wsum

                endif
              else
                xy=xx*yy
                expon=xx2*w2+yy2*w1-2.*xy*w12
                if(expon.le.100.)then
                  weight=exp(-0.5*expon)
                  ymod=(data(i,j)-sky(kk))*weight
                  sumxx=sumxx+xx2*ymod
                  sumyy=sumyy+yy2*ymod
                  sumxy=sumxy+xy*ymod
                  sums4=sums4+expon*expon*ymod
                  sum=sum+ymod
                endif
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



