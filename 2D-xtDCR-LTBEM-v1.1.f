c  ------------------------------------------------------------------------
      program main
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      parameter (pi=3.141592654d0)
      parameter (ns=12)
      parameter (ntime=4)
      parameter (ncpx=5,ncpy=10)
      character*30 outfil,infil,pr
      dimension nsse(10),kode(nto),xsil(10),ysil(10),xsiu(10),ysiu(10),
     &          sqlo1(nto),sqlo2(nto),sqhi1(nto),sqhi2(nto)
      dimension sqm1(nto),sqm2(nto)
      dimension sqn1(nto),sqn2(nto)
      dimension aa(nto,nto),bb(nto),xx(nto),cc(nto),ff(nto),bc(nto)
      dimension ttt(ntime)
      dimension xxint(19), yyint(19)
      dimension ccb(nto),ffb(nto)
      dimension cci(19,19), cci1(19,19), cci2(19,19)
      dimension v(ns)
      common /moduli/ d11,d12,d22
      common /koord/ sqlo1,sqlo2,sqhi1,sqhi2
      common /normal/ sqn1,sqn2
      common /tengah/ sqm1,sqm2
      common /constantD/ CapD
      common /vektorv/ v1, v2
      common /laplace_constant/ s
      common /alpha/ ralpha
      common /lambda/ rlambda
      common /constantKsmall/ SmallK
      common /nequ/ nequ
      data ttt / .5d0, 1.d0, 1.50, 2.00 /
      data xxint /.05d0,.1d0,.15d0,.2d0,.25d0,.3d0,.35d0,.4d0,.45d0,.5d0
     &           ,.55d0,.6d0,.65d0,.7d0,.75d0,.8d0,.85d0,.9d0,.95d0/
      data yyint /.05d0,.1d0,.15d0,.2d0,.25d0,.3d0,.35d0,.4d0,.45d0,.5d0
     &           ,.55d0,.6d0,.65d0,.7d0,.75d0,.8d0,.85d0,.9d0,.95d0/
      call cpu_time ( t1 )

c  CHOOSING THE EQUATION TYPE
 1    print *,"Please choose the number of equation type, and then press
     & enter:"
      print *,"         1. DR (Diffusion-Reaction) equation"
      print *,"         2. DCR (Diffusion-Convection-Reaction) equation"
      read(*,*) nequ
      if (nequ.eq.1) then
      print *,"You choose",nequ,"DR equation"
      elseif (nequ.eq.2) then
      print *,"You choose",nequ,"DCR equation"
      else
      go to 1
      endif
      print *,""

c  INPUT AND OUTPUT FILES
      outfil='2D-xtDCR-LTBEM-v1.1.out'
      open(unit=1,file=outfil,status='old')
      infil='2D-xtDCR-LTBEM-v1.1.dat'
      OPEN(unit=5,file=infil,status='old')

c  DISCRETISATION OF THE BOUNDARY
      nn=0
      read(5,*) nsides
      do i=1,nsides
        read(5,*) nsse(i),xsil(i),ysil(i),xsiu(i),ysiu(i)
 33     format(5x,i2,10x,i3,13x,4f5.2)
        nn=nn+nsse(i)
        dx1=(xsiu(i)-xsil(i))/dfloat(nsse(i))
        dx2=(ysiu(i)-ysil(i))/dfloat(nsse(i))
        do j=1,nsse(i)
          k=k+1
          sqlo1(k)=xsil(i)+dfloat(j-1)*dx1
          sqlo2(k)=ysil(i)+dfloat(j-1)*dx2
          sqhi1(k)=xsil(i)+dfloat(j)*dx1
          sqhi2(k)=ysil(i)+dfloat(j)*dx2
          sqm1(k)=sqlo1(k)+(sqhi1(k)-sqlo1(k))/2.d0
          sqm2(k)=sqlo2(k)+(sqhi2(k)-sqlo2(k))/2.d0
          ddx=sqhi1(k)-sqlo1(k)
          ddy=sqhi2(k)-sqlo2(k)
          dds=dsqrt( ddx*ddx+ddy*ddy )
          sqn1(k)=ddy/dds
          sqn2(k)=-ddx/dds
        enddo
      enddo
c  Calculate Stehfest coefficients v_m
      ns2=ns/2
      do m=1,ns
        if (m.le.ns2) then
          ks=m
        else
          ks=ns2
        endif
        v(m)=.0d0
        do k=int(m+1)/2,ks
          v(m)=v(m)+k**ns2*fact(2*k)
     &      /fact(ns2-k)/fact(k)/fact(k-1)/fact(m-k)/fact(2*k-m)
        enddo
        v(m)=v(m)*(-1)**(ns2+m)
      enddo
      trmse = 0.0d0
      trmse1 = 0.0d0
      trmse2 = 0.0d0
      tmae = 0.0d0
      tmae1 = 0.0d0
      tmae2 = 0.0d0
      tmape = 0.0d0
      tmape1 = 0.0d0
      tmape2 = 0.0d0
      crmse = 0.0d0
      crmse1 = 0.0d0
      crmse2 = 0.0d0
      cmae = 0.0d0
      cmae1 = 0.0d0
      cmae2 = 0.0d0
      cmape = 0.0d0
      cmape1 = 0.0d0
      cmape2 = 0.0d0
c  TIME t ITERATION
      do kkk=1,ntime
        tau=ttt(kkk)
        do ii=1,19
          do jj=1,19
            cci(ii,jj)  = 0.0d0
            cci1(ii,jj) = 0.0d0
            cci2(ii,jj) = 0.0d0
          enddo
        enddo
        do k=1,nn
          ccb(k) = 0.0d0
          ffb(k) = 0.0d0
        enddo
c  STEHFEST ITERATION
        print *,"Stehfest iterations for t =", tau
        do nnnnn=1,ns
          s = real(nnnnn)*dlog(2.d0)/tau
      if (nequ.eq.1) then
          d11 = 1.d0/(1.d0+s)
          d12 = .35d0*d11
          d22 = .15d0*d11
          call findroot(d11,d12,d22)
          reactk = 1.d0/s**2.d0
          rlambda = -.0045d0/(1.d0+s)
          ralpha = (.0045d0*(s+221.2176813d0)*(s+1.004540961d0))
     &             / s**3.d0 / (1.d0+s)
      else
          d11=1.d0/(1.d0+s)
          d12=.25d0*d11
          d22=.25d0*d11
          v1 = 1.d0/s - .8862269255d0*dexp(-.25d0/s)/(s**1.5d0)
          v2 = 2.d0*v1
          reactk = .5d0*dsqrt(pi/s**3.d0)*dexp(.25d0/s)
          rlambda = .0075d0/(1.d0+s) - .1d0/s
     &      + .08862269255d0*dexp(-.25d0/s)/s**1.5d0
          ralpha = -( .8862269255d0*dexp(.25d0/s)*(s+1.d0)
     &    - 0.06d0*s**1.5d0 - .1d0*s**.5d0 +
     &    .08862269255d0*dexp(-.25d0/s)*(s+1.d0)  )
     &    / (s**2.5d0*(1.d0+s))
      endif
          call findroot(d11,d12,d22)
c  THE COEFFICIENT k IN THE GOVERNING EQUATION
      if (nequ.eq.1) then
          SmallK=reactk-(s*ralpha)-rlambda
      else
          SmallK=rlambda+reactk+s*ralpha
      endif
          print *,"Iteration",nnnnn
      if (nequ.eq.1) then
        call fundparDR(d11,d12,d22,SmallK)
      else
        call fundparDCR(d11,d12,d22,SmallK)
      endif
c  THE BOUNDARY DATA
          do k=1,nn
            if (k.le.nsse(1)+nsse(2)+nsse(3)) then
              kode(k)=0
              bc(k)=Fexs(k,sqm1(k),sqm2(k))
            else
              kode(k)=1
              bc(k)=cexs(sqm1(k),sqm2(k))
            endif
          enddo
c  SETUP OF THE LINEAR SYSTEM OF EQUATIONS
          do k=1,nn
            bb(k)=0.d0
            spt1=sqm1(k)
            spt2=sqm2(k)
            do l=1,nn
      if (nequ.eq.1) then
              call bodeDR(99,l,spt1,spt2,smp,smg,pkiphi,smp1,smg1,
     &             pkiphi1,smp2,smg2,pkiphi2)
              if (kode(l).eq.0) then
                aa(k,l) = -(smg-pkiphi)
                bb(k) = bb(k)-smp*bc(l)
              else
                aa(k,l) = smp
                bb(k) = bb(k)+(smg-pkiphi)*bc(l)
              endif
      else
              call bodeDCR(99,l,spt1,spt2,smp,smg,smpv,smp1,smg1,smpv1,
     &               smp2,smg2,smpv2)
              if (kode(l).eq.0) then
                aa(k,l)=-(smpv-smg)
                bb(k)=bb(k)+smp*bc(l)
              else
                aa(k,l)=-smp
                bb(k)=bb(k)+(smpv-smg)*bc(l)
              endif
      endif
            enddo
            if (kode(k).eq.1) then
              bb(k)=bb(k)-0.5d0*bc(k)*h12(sqm1(k),sqm2(k))
            else
              aa(k,k)=aa(k,k)+0.5d0*h12(sqm1(k),sqm2(k))
            endif
          enddo
c  SOLVING THE LINEAR ALGEBRAIC EQUATION SYSTEM
          call laes(nn,aa,bb,xx)
cc  SOLUTIONS ON THE BOUNDARY
          do k=1,nn
            if (kode(k) .eq. 1) then
              cc(k)=bc(k)
              ff(k)=xx(k)
              ccb(k) = ccb(k) + ( cc(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
              ffb(k) = ffb(k) + ( ff(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
            else
              ff(k)=bc(k)
              cc(k)=xx(k)
              ccb(k) = ccb(k) + ( cc(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
              ffb(k) = ffb(k) + ( ff(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
            endif
          enddo
          do ii=1,19
            xint=xxint(ii)
            do jj=1,19
              yint=yyint(jj)
              psi=0.d0
              psi1=0.d0
              psi2=0.d0
              do l=1,nn
      if (nequ.eq.1) then
                call bodeDR(1,l,xint,yint,smp,smg,pkiphi,smp1,smg1,
     &               pkiphi1,smp2,smg2,pkiphi2)
                psi = psi+(smg-pkiphi)*cc(l)-smp*ff(l)
                psi1 = psi1+(smg1-pkiphi1)*cc(l)-smp1*ff(l)
                psi2 = psi2+(smg2-pkiphi2)*cc(l)-smp2*ff(l)
      else
                call bodeDCR(1,l,xint,yint,smp,smg,smpv,smp1,smg1,smpv1,
     &               smp2,smg2,smpv2)
                psi  = psi + smp*ff(l) + (smpv-smg)*cc(l)
                psi1 = psi1 + smp1*ff(l) + (smpv1-smg1)*cc(l)
                psi2 = psi2 + smp2*ff(l) + (smpv2-smg2)*cc(l)
      endif
              enddo
              c=psi/h12(xint,yint)
              dcdx1=(psi1-c*dh12dx1(xint,yint))/h12(xint,yint)
              dcdx2=(psi2-c*dh12dx2(xint,yint))/h12(xint,yint)
              cci(ii,jj)  = cci(ii,jj) + (c *v(nnnnn)
     &                         *dlog(2.0d0)/tau)
              cci1(ii,jj) = cci1(ii,jj)+ (dcdx1*v(nnnnn)
     &                         *dlog(2.0d0)/tau)
              cci2(ii,jj) = cci2(ii,jj)+ (dcdx2*v(nnnnn)
     &                         *dlog(2.0d0)/tau)
            enddo
          enddo
        enddo
        write(6,*)
c  END OF STEHFEST ITERATION
 117    format('#',2x,186('-'))
        write(6,235)
        write(1,235)
 235    format('# Solutions :')
        write(6,117)
        write(1,117)
        write(6,236)
        write(1,236)
 236    format('# Time t (1) | Point (x,y) (2-3) | BEM sol (4-6): c dc/d
     &x dc/dy | Anal sol (7-9): c dc/dx dc/dy | Err Square (10-12): c dc
     &/dx dc/dy | Abs Err (13-15): c dc/dx dc/dy | Abs Rel Err (16-18):
     &c dc/dx dc/dy')
        write(6,117)
        write(1,117)
        ermse = 0.0d0
        ermse1 = 0.0d0
        ermse2 = 0.0d0
        emae = 0.0d0
        emae1 = 0.0d0
        emae2 = 0.0d0
        emape = 0.0d0
        emape1 = 0.0d0
        emape2 = 0.0d0
        do ii=1,19
          xxi=xxint(ii)
          do jj=1,19
            yyi=yyint(jj)
       write(6,237)
     &   tau,xxi,yyi,cci(ii,jj),cci1(ii,jj),cci2(ii,jj),
     &   cext(xxi,yyi,tau), cext1(xxi,yyi,tau), cext2(xxi,yyi,tau),
     &   (cci(ii,jj)-cext(xxi,yyi,tau))**2.d0,
     &   (cci1(ii,jj)-cext1(xxi,yyi,tau))**2.d0,
     &   (cci2(ii,jj)-cext2(xxi,yyi,tau))**2.d0,
     &   dabs(cci(ii,jj)-cext(xxi,yyi,tau)),
     &   dabs(cci1(ii,jj)-cext1(xxi,yyi,tau)),
     &   dabs(cci2(ii,jj)-cext2(xxi,yyi,tau)),
     &   dabs( (cci(ii,jj)-cext(xxi,yyi,tau)) / cext(xxi,yyi,tau) ),
     &   dabs( (cci1(ii,jj)-cext1(xxi,yyi,tau)) / cext1(xxi,yyi,tau) ),
     &   dabs( (cci2(ii,jj)-cext2(xxi,yyi,tau)) / cext2(xxi,yyi,tau) )
       write(1,237)
     &   tau,xxi,yyi,cci(ii,jj),cci1(ii,jj),cci2(ii,jj),
     &   cext(xxi,yyi,tau), cext1(xxi,yyi,tau), cext2(xxi,yyi,tau),
     &   (cci(ii,jj)-cext(xxi,yyi,tau))**2.d0,
     &   (cci1(ii,jj)-cext1(xxi,yyi,tau))**2.d0,
     &   (cci2(ii,jj)-cext2(xxi,yyi,tau))**2.d0,
     &   dabs(cci(ii,jj)-cext(xxi,yyi,tau)),
     &   dabs(cci1(ii,jj)-cext1(xxi,yyi,tau)),
     &   dabs(cci2(ii,jj)-cext2(xxi,yyi,tau)),
     &   dabs( (cci(ii,jj)-cext(xxi,yyi,tau)) / cext(xxi,yyi,tau) ),
     &   dabs( (cci1(ii,jj)-cext1(xxi,yyi,tau)) / cext1(xxi,yyi,tau) ),
     &   dabs( (cci2(ii,jj)-cext2(xxi,yyi,tau)) / cext2(xxi,yyi,tau) )
 237    format(f7.4,1x,f5.2,1x,f5.2,1x,f12.6,1x,f12.6,1x,f12.6,1x,f12.6,
     &   1x,f12.6,1x,f12.6,1x,f9.6,1x,f9.6,1x,f9.6,1x,f9.6,1x,f9.6,
     &   1x,f9.6,1x,f9.6,1x,f9.6,1x,f9.6)
          ermse = ermse
     &                + (cci(ii,jj)-cext(xxi,yyi,tau))**2.d0
          ermse1 = ermse1
     &                 + (cci1(ii,jj)-cext1(xxi,yyi,tau))**2.d0
          ermse2 = ermse2
     &                 + (cci2(ii,jj)-cext2(xxi,yyi,tau))**2.d0
          emae = emae
     &               + dabs(cci(ii,jj)-cext(xxi,yyi,tau))
          emae1 = emae1
     &                + dabs(cci1(ii,jj)-cext1(xxi,yyi,tau))
          emae2 = emae2
     &                + dabs(cci2(ii,jj)-cext2(xxi,yyi,tau))
          emape = emape
     &    + dabs( (cci(ii,jj)-cext(xxi,yyi,tau)) / cext(xxi,yyi,tau) )
          emape1 = emape1
     &    + dabs( (cci1(ii,jj)-cext1(xxi,yyi,tau))/cext1(xxi,yyi,tau) )
          emape2 = emape2
     &    + dabs( (cci2(ii,jj)-cext2(xxi,yyi,tau))/cext2(xxi,yyi,tau) )
          trmse = trmse
     &                + (cci(ii,jj)-cext(xxi,yyi,tau))**2.d0
          trmse1 = trmse1
     &                + (cci1(ii,jj)-cext1(xxi,yyi,tau))**2.d0
          trmse2 = trmse2
     &                + (cci2(ii,jj)-cext2(xxi,yyi,tau))**2.d0
          tmae = tmae
     &               + dabs(cci(ii,jj)-cext(xxi,yyi,tau))
          tmae1 = tmae1
     &               + dabs(cci1(ii,jj)-cext1(xxi,yyi,tau))
          tmae2 = tmae2
     &               + dabs(cci2(ii,jj)-cext2(xxi,yyi,tau))
          tmape = tmape
     &    + dabs( (cci(ii,jj)-cext(xxi,yyi,tau)) / cext(xxi,yyi,tau) )
          tmape1 = tmape1
     &    + dabs( (cci1(ii,jj)-cext1(xxi,yyi,tau))/cext1(xxi,yyi,tau) )
          tmape2 = tmape2
     &    + dabs( (cci2(ii,jj)-cext2(xxi,yyi,tau))/cext2(xxi,yyi,tau) )
      if (ii.eq.ncpx .and. jj.eq.ncpy) then
        crmse = crmse +
     &    (cci(ii,jj)-cext(xxi,yyi,tau))**2.d0
        crmse1 = crmse1 +
     &    (cci1(ii,jj)-cext1(xxi,yyi,tau))**2.d0
        crmse2 = crmse2 +
     &    (cci2(ii,jj)-cext2(xxi,yyi,tau))**2.d0
        cmae = cmae +
     &    dabs(cci(ii,jj)-cext(xxi,yyi,tau))
        cmae1 = cmae1 +
     &    dabs(cci1(ii,jj)-cext1(xxi,yyi,tau))
        cmae2 = cmae2 +
     &    dabs(cci2(ii,jj)-cext2(xxi,yyi,tau))
        cmape = cmape +
     &    dabs( (cci(ii,jj)-cext(xxi,yyi,tau)) / cext(xxi,yyi,tau) )
        cmape1 = cmape1 +
     &    dabs( (cci1(ii,jj)-cext1(xxi,yyi,tau)) / cext1(xxi,yyi,tau) )
        cmape2 = cmape2 +
     &    dabs( (cci2(ii,jj)-cext2(xxi,yyi,tau)) / cext2(xxi,yyi,tau) )
      else
      endif
          enddo
        enddo
        write(6,117)
        write(1,117)
        write(6,444)
        write(1,444)
 444    format('# Time t (1) | RMSE (2-4): c dc/dx dc/dy | MAE (5-7): c
     &dc/dx dc/dy | MAPE (8-10): c dc/dx dc/dy')
        write(1,999) tau,dsqrt(ermse/19.d0/19.d0),
     &                 dsqrt(ermse1/19.d0/19.d0),
     &                 dsqrt(ermse2/19.d0/19.d0),
     &                 emae/19.d0/19.d0,
     &                 emae1/19.d0/19.d0,
     &                 emae2/19.d0/19.d0,
     &                 emape/19.d0/19.d0,
     &                 emape1/19.d0/19.d0,
     &                 emape2/19.d0/19.d0
        write(6,999) tau,dsqrt(ermse/19.d0/19.d0),
     &                 dsqrt(ermse1/19.d0/19.d0),
     &                 dsqrt(ermse2/19.d0/19.d0),
     &                 emae/19.d0/19.d0,
     &                 emae1/19.d0/19.d0,
     &                 emae2/19.d0/19.d0,
     &                 emape/19.d0/19.d0,
     &                 emape1/19.d0/19.d0,
     &                 emape2/19.d0/19.d0
 999    format(f7.4,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,
     &         f10.6,1x,f10.6,1x,f10.6,1x,f10.6)
        write(6,*)
        write(1,*)
      enddo
c  END OF TIME t ITERATION
      call cpu_time ( t2 )
        write(6,234)
        write(1,234)
 234    format('# For all time-steps and all points:')
      write(6,777)
      write(1,777)
 777  format('# N (1) | CPU time (2) | RMSE (3-5): c dc/dx dc/dy | MAE (
     &6-8): c dc/dx dc/dy | MAPE (9-11): c dc/dx dc/dy')
      write(1,998) ns, t2-t1,
     &   dsqrt(trmse/19.d0/19.d0/ntime),
     &   dsqrt(trmse1/19.d0/19.d0/ntime),
     &   dsqrt(trmse2/19.d0/19.d0/ntime),
     &   tmae/19.d0/19.d0/ntime,
     &   tmae1/19.d0/19.d0/ntime,
     &   tmae2/19.d0/19.d0/ntime,
     &   tmape/19.d0/19.d0/ntime,
     &   tmape1/19.d0/19.d0/ntime,
     &   tmape2/19.d0/19.d0/ntime
      write(6,998) ns, t2-t1,
     &   dsqrt(trmse/19.d0/19.d0/ntime),
     &   dsqrt(trmse1/19.d0/19.d0/ntime),
     &   dsqrt(trmse2/19.d0/19.d0/ntime),
     &   tmae/19.d0/19.d0/ntime,
     &   tmae1/19.d0/19.d0/ntime,
     &   tmae2/19.d0/19.d0/ntime,
     &   tmape/19.d0/19.d0/ntime,
     &   tmape1/19.d0/19.d0/ntime,
     &   tmape2/19.d0/19.d0/ntime
 998  format(1x,i2,3x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,
     &          f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6)
        write(6,345) xxint(ncpx),yyint(ncpy)
        write(1,345) xxint(ncpx),yyint(ncpy)
 345    format('# For all time-steps at point (',f5.2,',',f5.2,'):')
      write(6,666)
      write(1,666)
 666  format('# N (1) | RMSE (2-4): c dc/dx dc/dy | MAE (5-7): c dc/dx d
     &c/dy | MAPE (8-10): c dc/dx dc/dy')
      write(1,997) ns,
     &   dsqrt(crmse/ntime),
     &   dsqrt(crmse1/ntime),
     &   dsqrt(crmse2/ntime),
     &   cmae/ntime,
     &   cmae1/ntime,
     &   cmae2/ntime,
     &   cmape/ntime,
     &   cmape1/ntime,
     &   cmape2/ntime
      write(6,997) ns,
     &   dsqrt(crmse/ntime),
     &   dsqrt(crmse1/ntime),
     &   dsqrt(crmse2/ntime),
     &   cmae/ntime,
     &   cmae1/ntime,
     &   cmae2/ntime,
     &   cmape/ntime,
     &   cmape1/ntime,
     &   cmape2/ntime
 997  format(1x,i2,3x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,
     &          f10.6,1x,f10.6,1x,f10.6,1x,f10.6)
        write(6,*)
        write(1,*)
        write(6,567)
        write(1,567)
 567    format('# RMSE = Root Mean Square Error, MAE = Mean Absolute Err
     &or, MAPE = Mean Absolute Percentage Error')
      end
c  ------------------------------------------------------------------------
c  Finding roots
      subroutine findroot(d11,d12,d22)
      implicit real*8 (a-h, o-z)
      complex*8 rho
      common /akar/ rhodot, rhoddot
      if (d11.le.0.d0 .or. d22.le.0.d0) then
        write(6,*) '  k11 and k22 should be both positive. '
        write(1,*) '  k11 and k22 must be both positive. '
        stop
      else
      endif
      det=dble(d11*d22-d12*d12)
      if (det.lt.0.d0) then
        write(6,*) '  d11.d22-d12.d12 =',det
        write(1,*) '  d11.d22-d12.d12 =',det
        write(6,*) '  d11.d22-d12.d12 should be positive.'
        write(1,*) '  d11.d22-d12.d12 must be positive.'
        stop
      else
      endif
      rho=cmplx(-d12/d22,dsqrt(det)/d22)
      rhodot = real(rho)
      rhoddot = aimag(rho)
      return
      end
c  ------------------------------------------------------------------------
c  Parameters necessary for the evaluation of the fund. solution
      subroutine fundparDCR(d11,d12,d22,SmallK)
      implicit real*8 (a-h, o-z)
      complex*8 rho
      dimension VecV(2)
      common /vektorv/ v1, v2
      common /akar/ rhodot, rhoddot
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      common /laplace_constant/ s
      common /alpha/ ralpha
      common /lambda/ rlambda
      CapD = .5d0*(d11+2.d0*rhodot*d12+
     &            (rhodot**2.d0+rhoddot**2.d0)*d22)
      CapK = rhoddot/CapD
      VecV(1) = (v1 + rhodot*v2)
      VecV(2) = (rhoddot*v2)
      Cmu = dsqrt( (VecV(1)**2.d0+VecV(2)**2.d0)/(4.d0*CapD**2.0d0)
     &      + SmallK/CapD )
      return
      end
c  ------------------------------------------------------------------------
c  Parameters necessary for the evaluation of the fund. solution
      subroutine fundparDR(d11,d12,d22,SmallK)
      implicit real*8 (a-h, o-z)
      common /akar/ rhodot, rhoddot
      common /constantAzis/ Azis
      common /constantOmega / Omega
      dummy = d11+2.d0*rhodot*d12+d22*(rhodot**2.d0+
     &        rhoddot**2.d0)
      Omega = dsqrt( 2.d0*dabs(SmallK)/dummy )
      Azis = 1.d0/rhoddot/d22
      return
      end
c  ------------------------------------------------------------------------
c  NUMERICAL INTEGRATION USING BODE'S RULE
      subroutine bodeDCR(ints,li,spt1,spt2,smp,smg,smpv,smp1,smg1,
     &                 smpv1,smp2,smg2,smpv2)
      implicit real*8 (a-h, o-z)
      parameter (nto=1000)
      dimension qlo1(nto),qlo2(nto),qhi1(nto),qhi2(nto),
     &          qn1(nto),qn2(nto)
      dimension w(5)
      complex*8 qkm1,qk,vectx,vecty
      common /vektorv/ v1, v2
      common /koord/ qlo1,qlo2,qhi1,qhi2
      common /normal/ qn1, qn2
      external Phi,Gamma,Phi1,Gamma1,Phi2,Gamma2,Fh
      data w / 2857.0d0, 15741.0d0, 1080.0d0, 19344.0d0, 5778.0d0 /
      smp=0.d0
      smg=0.d0
      smpv=0.d0
      smp1=0.d0
      smg1=0.d0
      smpv1=0.d0
      smp2=0.d0
      smg2=0.d0
      smpv2=0.d0
      qkm1=cmplx(qlo1(li),qlo2(li))
      qk=cmplx(qhi1(li),qhi2(li))
      do kkk=1,5
        vectx=qkm1+dfloat(kkk-1)*(qk-qkm1)/9.d0
        x1=real(vectx)
        x2=aimag(vectx)
        vecty=qkm1+dfloat(10-kkk)*(qk-qkm1)/9.d0
        y1=real(vecty)
        y2=aimag(vecty)
        pkix = Fh(li,x1,x2)
        pkiy = Fh(li,y1,y2)
        smp = smp + ( Phi(x1,x2,spt1,spt2)/h12(x1,x2) +
     &                Phi(y1,y2,spt1,spt2)/h12(y1,y2) ) * w(kkk)
        smg = smg + ( Gamma(li,x1,x2,spt1,spt2)*h12(x1,x2) +
     &                Gamma(li,y1,y2,spt1,spt2)*h12(y1,y2) ) * w(kkk)
        smpv = smpv + ( (pkix - (v1*qn1(li)+v2*qn2(li))*h12(x1,x2)) *
     &                  Phi(x1,x2,spt1,spt2)
     &              +   (pkiy - (v1*qn1(li)+v2*qn2(li))*h12(y1,y2)) *
     &                  Phi(y1,y2,spt1,spt2) ) * w(kkk)
        if (ints.eq.1) then
        smp1 = smp1 + ( Phi1(x1,x2,spt1,spt2)/h12(x1,x2) +
     &                  Phi1(y1,y2,spt1,spt2)/h12(y1,y2) ) * w(kkk)
        smg1 = smg1 + ( Gamma1(li,x1,x2,spt1,spt2)*h12(x1,x2) +
     &                  Gamma1(li,y1,y2,spt1,spt2)*h12(y1,y2) ) * w(kkk)
        smpv1 = smpv1 + ( (pkix - (v1*qn1(li)+v2*qn2(li))*h12(x1,x2)) *
     &                    Phi1(x1,x2,spt1,spt2)
     &                +   (pkiy - (v1*qn1(li)+v2*qn2(li))*h12(y1,y2)) *
     &                    Phi1(y1,y2,spt1,spt2) ) * w(kkk)
        smp2 = smp2 + ( Phi2(x1,x2,spt1,spt2)/h12(x1,x2) +
     &                  Phi2(y1,y2,spt1,spt2)/h12(y1,y2) ) * w(kkk)
        smg2 = smg2 + ( Gamma2(li,x1,x2,spt1,spt2)*h12(x1,x2) +
     &                  Gamma2(li,y1,y2,spt1,spt2)*h12(y1,y2) ) * w(kkk)
        smpv2 = smpv2 + ( (pkix - (v1*qn1(li)+v2*qn2(li))*h12(x1,x2)) *
     &                    Phi2(x1,x2,spt1,spt2)
     &                +   (pkiy - (v1*qn1(li)+v2*qn2(li))*h12(y1,y2)) *
     &                    Phi2(y1,y2,spt1,spt2) ) * w(kkk)
        else
        endif
      enddo
      h=abs(qkm1-qk)
      smp=smp*h/89600.0d0
      smg=smg*h/89600.0d0
      smpv=smpv*h/89600.0d0
      if (ints.eq.1) then
        smp1=smp1*h/89600.0d0
        smg1=smg1*h/89600.0d0
        smpv1=smpv1*h/89600.0d0
        smp2=smp2*h/89600.0d0
        smg2=smg2*h/89600.0d0
        smpv2=smpv2*h/89600.0d0
      else
      endif
      return
      end
c  ------------------------------------------------------------------------
c  NUMERICAL INTEGRATION USING BODE'S RULE
      subroutine bodeDR(ints,li,xi1,xi2,smp,smg,pkiphi,smp1,smg1,
     &                  pkiphi1,smp2,smg2,pkiphi2)
      implicit real*8 (a-h, o-z)
      parameter (nto=1000)
      dimension qlo1(nto),qlo2(nto),qhi1(nto),qhi2(nto)
      dimension w(5)
      complex*8 qkm1,qk,vectx,vecty
      common /koord/ qlo1,qlo2,qhi1,qhi2
      common /constantKsmall/ SmallK
      external h12,Phineg,Phipos,Gammaneg,Gammapos,Phineg1,Phipos1,
     &       Gammaneg1,Gammapos1,Phineg2,Phipos2,Gammaneg2,Gammapos2,
     &       Phizero,Gammazero,Phizero1,Gammazero1,Phizero2,Gammazero2
      data w / 2857.0d0, 15741.0d0, 1080.0d0, 19344.0d0, 5778.0d0 /
      external Fh
      pkiphi = 0.d0
      smp = 0.d0
      smg = 0.d0
      pkiphi1 = 0.d0
      smp1 = 0.d0
      smg1 = 0.d0
      pkiphi2 = 0.d0
      smp2 = 0.d0
      smg2 = 0.d0
      qkm1 = cmplx(qlo1(li),qlo2(li))
      qk = cmplx(qhi1(li),qhi2(li))
      do kkk=1,5
        vectx = qkm1+dfloat(kkk-1)*(qk-qkm1)/9.d0
        x1 = real(vectx)
        x2 = aimag(vectx)
        vecty = qkm1+dfloat(10-kkk)*(qk-qkm1)/9.d0
        y1 = real(vecty)
        y2 = aimag(vecty)
        pkix = Fh(li,x1,x2)
        pkiy = Fh(li,y1,y2)
        if (SmallK.lt.0.d0) then
          pkiphi = pkiphi+(Phineg(x1,x2,xi1,xi2)*pkix+
     &             Phineg(y1,y2,xi1,xi2)*pkiy)*w(kkk)
          smp = smp+(Phineg(x1,x2,xi1,xi2)/h12(x1,x2)+
     &          Phineg(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
          smg = smg+(Gammaneg(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &          Gammaneg(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
          if (ints.eq.1) then
            pkiphi1 = pkiphi1+(Phineg1(x1,x2,xi1,xi2)*pkix+
     &                Phineg1(y1,y2,xi1,xi2)*pkiy)*w(kkk)
            smp1 = smp1+(Phineg1(x1,x2,xi1,xi2)/h12(x1,x2)+
     &             Phineg1(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
            smg1 = smg1+(Gammaneg1(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &             Gammaneg1(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
            pkiphi2 = pkiphi2+(Phineg2(x1,x2,xi1,xi2)*pkix+
     &                Phineg2(y1,y2,xi1,xi2)*pkiy)*w(kkk)
            smp2 = smp2+(Phineg2(x1,x2,xi1,xi2)/h12(x1,x2)+
     &             Phineg2(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
            smg2 = smg2+(Gammaneg2(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &             Gammaneg2(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
          else
          endif
        elseif (SmallK.gt.0.d0) then
          pkiphi = pkiphi+(Phipos(x1,x2,xi1,xi2)*pkix+
     &             Phipos(y1,y2,xi1,xi2)*pkiy)*w(kkk)
          smp = smp+(Phipos(x1,x2,xi1,xi2)/h12(x1,x2)+
     &          Phipos(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
          smg = smg+(Gammapos(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &          Gammapos(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
          if (ints.eq.1) then
            pkiphi1 = pkiphi1+(Phipos1(x1,x2,xi1,xi2)*pkix+
     &                Phipos1(y1,y2,xi1,xi2)*pkiy)*w(kkk)
            smp1 = smp1+(Phipos1(x1,x2,xi1,xi2)/h12(x1,x2)+
     &             Phipos1(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
            smg1 = smg1+(Gammapos1(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &             Gammapos1(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
            pkiphi2 = pkiphi2+(Phipos2(x1,x2,xi1,xi2)*pkix+
     &                Phipos2(y1,y2,xi1,xi2)*pkiy)*w(kkk)
            smp2 = smp2+(Phipos2(x1,x2,xi1,xi2)/h12(x1,x2)+
     &             Phipos2(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
            smg2 = smg2+(Gammapos2(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &             Gammapos2(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
          else
          endif
        else
          pkiphi = pkiphi+(Phizero(x1,x2,xi1,xi2)*pkix+
     &             Phizero(y1,y2,xi1,xi2)*pkiy)*w(kkk)
          smp = smp+(Phizero(x1,x2,xi1,xi2)/h12(x1,x2)+
     &          Phizero(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
          smg = smg+(Gammazero(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &          Gammazero(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
          if (ints.eq.1) then
            pkiphi1 = pkiphi1+(Phizero1(x1,x2,xi1,xi2)*pkix+
     &                Phizero1(y1,y2,xi1,xi2)*pkiy)*w(kkk)
            smp1 = smp1+(Phizero1(x1,x2,xi1,xi2)/h12(x1,x2)+
     &             Phizero1(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
            smg1 = smg1+(Gammazero1(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &             Gammazero1(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
            pkiphi2 = pkiphi2+(Phizero2(x1,x2,xi1,xi2)*pkix+
     &                Phizero2(y1,y2,xi1,xi2)*pkiy)*w(kkk)
            smp2 = smp2+(Phizero2(x1,x2,xi1,xi2)/h12(x1,x2)+
     &             Phizero2(y1,y2,xi1,xi2)/h12(y1,y2))*w(kkk)
            smg2 = smg2+(Gammazero2(li,x1,x2,xi1,xi2)*h12(x1,x2)+
     &             Gammazero2(li,y1,y2,xi1,xi2)*h12(y1,y2))*w(kkk)
          else
          endif
        endif
      enddo
      h = abs(qkm1-qk)
      pkiphi = pkiphi*h/89600.0d0
      smp = smp*h/89600.0d0
      smg = smg*h/89600.0d0
      if (ints.eq.1) then
        pkiphi1 = pkiphi1*h/89600.0d0
        smp1 = smp1*h/89600.0d0
        smg1 = smg1*h/89600.0d0
        pkiphi2 = pkiphi2*h/89600.0d0
        smp2 = smp2*h/89600.0d0
        smg2 = smg2*h/89600.0d0
      else
      endif
      return
      end
c  ------------------------------------------------------------------------
c  The function Fh
      function Fh(kp,x1,x2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension xn(nto),yn(nto)
      common /moduli/ d11,d12,d22
      common /normal/ xn,yn
      external dh12dx1, dh12dx2
      der1 = dh12dx1(x1,x2)
      der2 = dh12dx2(x1,x2)
      Fh = d11*der1*xn(kp) + d12*(der1*yn(kp)+der2*xn(kp))
     &     + d22*der2*yn(kp)
      return
      end
c  ------------------------------------------------------------------------
c  A subroutine for solving a linear algebraic system AX=B
      subroutine laes(nn,aa,bb,xx)
      implicit real*8 (a-h, o-z)
      parameter (nto=1000)
      dimension aa(nto,nto),bb(nto),xx(nto),wk(nto),bp(nto),ap(nto),
     &       cp(nto)
      ig=nn
      jg=nn
      itt=0
  60  do i=1,ig
        ape=0.d0
        do j=1,ig
          ape=aa(i,j)*xx(j)+ape
        enddo
        wk(i)=bb(i)-ape
      enddo
      do i=1,ig
        ape=0.d0
        do j=1,ig
          ape=ape+aa(j,i)*wk(j)
        enddo
        ap(i)=ape
        cp(i)=ape
      enddo
      irr=0
  94  itt=itt+1
      irr=irr+1
      do i=1,ig
        bp(i)=0.d0
      enddo
      do j=1,ig
        do i=1,ig
          bp(i)=aa(i,j)*ap(j)+bp(i)
        enddo
      enddo
      daa=0.d0
      dab=0.d0
      do j=1,ig
        daa=daa+cp(j)*cp(j)
        dab=dab+bp(j)*bp(j)
      enddo
      if(dab.lt.1.0e-37.and.irr.eq.1) go to 115
      if(dab.lt.1.0e-37) go to 60
      dac=daa/dab
      do i=1,ig
        wk(i)=wk(i)-dac*bp(i)
        xx(i)=xx(i)+dac*ap(i)
      enddo
      err=0.d0
      do i=1,ig
        err=err+wk(i)*wk(i)
      enddo
      if(err.lt.1.0e-37.or.itt.eq.jg) go to 115
      if(irr.eq.ig) go to 60
      bong=0.d0
      do i=1,ig
        ape=0.d0
        do j=1,ig
          ape=ape+aa(j,i)*wk(j)
        enddo
        cp(i)=ape
        bong=bong+ape*ape
      enddo
      bong=bong/daa
      do i=1,ig
        ap(i)=cp(i)+bong*ap(i)
      enddo
      go to 94
 115  return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Phi
      function Phi(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2), VecR(2)
      common /akar/ rhodot, rhoddot
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      pi2 = 8.d0*datan(1.d0)
      r = Rcap(x1,x2,xi1,xi2)
      VecR(1) = (x1-xi1) + (x2-xi2)*rhodot
      VecR(2) = (x2-xi2)*rhoddot
      ProductVR = VecV(1)*VecR(1) + VecV(2)*VecR(2)
      x = Cmu*r
      FFK0 = FK0(x)
      Phi = CapK/pi2*dexp(-ProductVR/(2.d0*CapD))*FFK0
      return
      end
c  The derivative d(Phi)/d(a)
      function Phi1(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      vx = VecV(1)
      PPhi = Phi(x1,x2,xi1,xi2)
      dR1 = dRcap1(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi1 = (vx/(2.d0*CapD) - Cmu*dR1*FFK1/FFK0) * PPhi
      return
      end
c  The derivative d(Phi)/d(b)
      function Phi2(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      common /akar/ rhodot, rhoddot
      vx = VecV(1)
      vy = VecV(2)
      PPhi = Phi(x1,x2,xi1,xi2)
      dR2 = dRcap2(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi2 = ((vx*rhodot+vy*rhoddot)/(2.d0*CapD) - Cmu*dR2*FFK1/
     &         FFK0) * PPhi
      return
      end
c  The derivative d^2(Phi)/d(a^2)
      function Phi11(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      vx = VecV(1)
      PPhi = Phi(x1,x2,xi1,xi2)
      PPhi1 = Phi1(x1,x2,xi1,xi2)
      dR1 = dRcap1(x1,x2,xi1,xi2)
      dR11 = dRcap11(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi11 = PPhi1*(vx/(2.d0*CapD)-Cmu*FFK1/FFK0*dR1)
     &   -PPhi*Cmu*(Cmu*dR1**2.d0*(-1.d0-FFK1/FFK0/Cmu/r+
     &   FFK1**2.d0/FFK0**2.d0)+FFK1/FFK0*dR11)
      return
      end
c  The derivative d^2(Phi)/d(a)d(b)
      function Phi12(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      vx = VecV(1)
      PPhi = Phi(x1,x2,xi1,xi2)
      PPhi2 = Phi2(x1,x2,xi1,xi2)
      dR1 = dRcap1(x1,x2,xi1,xi2)
      dR2 = dRcap2(x1,x2,xi1,xi2)
      dR12 = dRcap12(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi12 = PPhi2*(vx/(2.d0*CapD)-Cmu*FFK1/FFK0*dR1)
     &   -PPhi*Cmu*(Cmu*dR1*dR2*(-1.d0-FFK1/FFK0/Cmu/r+
     &   FFK1**2.d0/FFK0**2.d0)+FFK1/FFK0*dR12)
      return
      end
c  The derivative d^2(Phi)/d(b^2)
      function Phi22(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /akar/ rhodot, rhoddot
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      vx = VecV(1)
      vy = VecV(2)
      PPhi = Phi(x1,x2,xi1,xi2)
      PPhi2 = Phi2(x1,x2,xi1,xi2)
      dR2 = dRcap2(x1,x2,xi1,xi2)
      dR22 = dRcap22(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi22 = PPhi2*((vx*rhodot+vy*rhoddot)/(2.d0*CapD)-Cmu*
     &   FFK1/FFK0*dR2)-PPhi*Cmu*(Cmu*dR2**2.d0*(-1.d0-FFK1/
     &   FFK0/Cmu/r+FFK1**2.d0/FFK0**2.d0)+FFK1/FFK0*dR22)
      return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Gamma
      function Gamma(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      dPhi1 = -Phi1(x1,x2,xi1,xi2)
      dPhi2 = -Phi2(x1,x2,xi1,xi2)
      Gamma = d11*dPhi1*qn1(ki) + d12*(dPhi1*qn2(ki)+dPhi2*qn1(ki))
     &        + d22*dPhi2*qn2(ki)
      return
      end
c  The derivative D(Gamma)/D(a)
      function Gamma1(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
c      external Phi11,Phi12
      dPhi11 = -Phi11(x1,x2,xi1,xi2)
      dPhi12 = -Phi12(x1,x2,xi1,xi2)
      Gamma1 = d11*dPhi11*qn1(ki) + d12*(dPhi11*qn2(ki)+
     &         dPhi12*qn1(ki)) + d22*dPhi12*qn2(ki)
      return
      end
c  The derivative D(Gamma)/D(b)
      function Gamma2(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      dPhi12 = -Phi12(x1,x2,xi1,xi2)
      dPhi22 = -Phi22(x1,x2,xi1,xi2)
      Gamma2 = d11*dPhi12*qn1(ki) + d12*(dPhi12*qn2(ki)+
     &         dPhi22*qn1(ki)) + d22*dPhi22*qn2(ki)
      return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Phi, for k<0
      function Phineg(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      common /constantOmega / Omega
      external Rcap,FK0
      pi = 4.d0*datan(1.d0)
      RR = Rcap(x1,x2,xi1,xi2)
      x = Omega*RR
      FFK0 = FK0(x)
      Phineg = -Azis*FFK0/(2.d0*pi)
      return
      end
c  The derivative d(Phi)/d(a), for k<0
      function Phineg1(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      common /constantOmega / Omega
      external Rcap,dRcap1,FK1
      pi = 4.d0*datan(1.d0)
      RRs = Rcap(x1,x2,xi1,xi2)
      dRRs1 = dRcap1(x1,x2,xi1,xi2)
      x = Omega*RRs
      FFK1 = FK1(x)
      Phineg1 = Azis*FFK1*Omega*dRRs1/(2.d0*pi)
      return
      end
c  The derivative d(Phi)/d(b), for k<0
      function Phineg2(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      common /constantOmega / Omega
      external Rcap,dRcap2,FK1
      pi = 4.d0*datan(1.d0)
      RRs = Rcap(x1,x2,xi1,xi2)
      dRRs2 = dRcap2(x1,x2,xi1,xi2)
      x = Omega*RRs
      FFK1 = FK1(x)
      Phineg2 = Azis*FFK1*Omega*dRRs2/(2.*pi)
      return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Phi, for k>0
      function Phipos(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      complex*8 iii, HH20, H20
      common /constantAzis/ Azis
      common /constantOmega / Omega
      external Rcap,H20
      iii = cmplx(0.d0,1.d0)
      RRc = Rcap(x1,x2,xi1,xi2)
      x = Omega*RRc
      HH20 = H20(x)
      Phipos = real( iii/4.d0*Azis*HH20 )
      return
      end
c  The derivative d(Phi)/d(a), for k>0.
c  NOTE :   d(Phi)/d(a) = - d(Phi)/d(x1)
      function Phipos1(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      complex*8 iii,H21,HH21
      common /constantAzis/ Azis
      common /constantOmega / Omega
      external Rcap,H21,dRcap1
      iii = cmplx(0.d0,1.d0)
      RRc = Rcap(x1,x2,xi1,xi2)
      dRRc1 = dRcap1(x1,x2,xi1,xi2)
      x = Omega*RRc
      HH21 = H21(x)
      Phipos1 = real( -iii/4.d0*Azis*HH21*Omega*dRRc1 )
      return
      end
c  The derivative d(Phi)/d(b), for k>0.
c  NOTE :   d(Phi)/d(b) = - d(Phi)/d(x2)
      function Phipos2(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      complex*8 iii,H21,HH21
      common /constantAzis/ Azis
      common /constantOmega / Omega
      external Rcap,H21,dRcap2
      iii = cmplx(0.d0,1.d0)
      RRc = Rcap(x1,x2,xi1,xi2)
      dRRc2 = dRcap2(x1,x2,xi1,xi2)
      x = Omega*RRc
      HH21 = H21(x)
      Phipos2 = real(-iii/4.d0*Azis*HH21*Omega*
     &                dRRc2 )
      return
      end
c  The derivative d^2(Phi)/d(a^2), for k>0.
      function Phipos11(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      complex*8 iii,H20,H21,HH20,HH21
      common /constantAzis/ Azis
      common /constantOmega / Omega
      iii = cmplx(0.d0,1.d0)
      RR = Rcap(x1,x2,xi1,xi2)
      dRR1 = dRcap1(x1,x2,xi1,xi2)
      dRR11 = dRcap11(x1,x2,xi1,xi2)
      x = Omega*RR
      HH20 = H20(x)
      HH21 = H21(x)
      Phipos11 = real( -iii/4.d0*Azis*Omega*((HH20-HH21/x)*Omega*dRR1*
     &                 dRR1+HH21*dRR11) )
      return
      end
c  The derivative d^2(Phi)/d(a)d(b), for k>0.
      function Phipos12(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      complex*8 iii,H20,H21,HH20,HH21
      common /constantAzis/ Azis
      common /constantOmega / Omega
      iii = cmplx(0.d0,1.d0)
      RR = Rcap(x1,x2,xi1,xi2)
      dRR1 = dRcap1(x1,x2,xi1,xi2)
      dRR2 = dRcap2(x1,x2,xi1,xi2)
      dRR12 = dRcap12(x1,x2,xi1,xi2)
      x = Omega*RR
      HH20 = H20(x)
      HH21 = H21(x)
      Phipos12 = real( -iii/4.d0*Azis*Omega*((HH20-HH21/x)*Omega*dRR1*
     &                 dRR2+HH21*dRR12) )
      return
      end
c  The derivative d^2(Phi)/d(b^2), for k>0.
      function Phipos22(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      complex*8 iii,H20,H21,HH20,HH21
      common /constantAzis/ Azis
      common /constantOmega / Omega
      iii = cmplx(0.d0,1.d0)
      RR = Rcap(x1,x2,xi1,xi2)
      dRR2 = dRcap2(x1,x2,xi1,xi2)
      dRR22 = dRcap22(x1,x2,xi1,xi2)
      x = Omega*RR
      HH20 = H20(x)
      HH21 = H21(x)
      Phipos22 = real( -iii/4.d0*Azis*Omega*((HH20-HH21/x)*Omega*dRR2*
     &                 dRR2+HH21*dRR22) )
      return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Phi, for k=0
      function Phizero(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      external Rcap
      pi2 = 8.d0*datan(1.d0)
      RRc = Rcap(x1,x2,xi1,xi2)
      Phizero = Azis*dlog(RRc)/pi2
      return
      end
c  The derivative d(Phi)/d(a), for k=0.
c  NOTE :   d(Phi)/d(a) = - d(Phi)/d(x1)
      function Phizero1(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      external Rcap,dRcap1
      pi2 = 8.d0*datan(1.d0)
      RRc = Rcap(x1,x2,xi1,xi2)
      dRRc1 = dRcap1(x1,x2,xi1,xi2)
      Phizero1 = Azis/RRc/pi2*dRRc1
      return
      end
c  The derivative d(Phi)/d(b), for k=0.
c  NOTE :   d(Phi)/d(b) = - d(Phi)/d(x2)
      function Phizero2(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      external Rcap,dRcap2
      pi2 = 8.d0*datan(1.d0)
      RRc = Rcap(x1,x2,xi1,xi2)
      dRRc2 = dRcap2(x1,x2,xi1,xi2)
      Phizero2 = Azis/RRc/pi2*dRRc2
      return
      end
c  The derivative d^2(Phi)/d(a^2), for k=0.
      function Phizero11(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      external Rcap,dRcap1,dRcap11
      pi2 = 8.d0*datan(1.d0)
      RR = Rcap(x1,x2,xi1,xi2)
      dRR1 = dRcap1(x1,x2,xi1,xi2)
      dRR11 = dRcap11(x1,x2,xi1,xi2)
      Phizero11 = Azis/pi2*(-dRR1**2.d0/RR**2.d0+dRR11/RR)
      return
      end
c  The derivative d^2(Phi)/d(a)d(b), for k=0.
      function Phizero12(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      external Rcap,dRcap1,dRcap2,dRcap12
      pi2 = 8.d0*datan(1.d0)
      RR = Rcap(x1,x2,xi1,xi2)
      dRR1 = dRcap1(x1,x2,xi1,xi2)
      dRR2 = dRcap2(x1,x2,xi1,xi2)
      dRR12 = dRcap12(x1,x2,xi1,xi2)
      Phizero12 = Azis/pi2*(-dRR1*dRR2/RR**2.d0+dRR12/RR)
      return
      end
c  The derivative d^2(Phi)/d(b^2), for k=0.
      function Phizero22(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /constantAzis/ Azis
      external Rcap,dRcap2,dRcap22
      pi2 = 8.d0*datan(1.d0)
      RR = Rcap(x1,x2,xi1,xi2)
      dRR2 = dRcap2(x1,x2,xi1,xi2)
      dRR22 = dRcap22(x1,x2,xi1,xi2)
      Phizero22 = Azis/pi2*(-dRR2**2.d0/RR**2.d0+dRR22/RR)
      return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Gamma, for k<0
      function Gammaneg(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Phineg1,Phineg2
      ddPhineg1 = -Phineg1(x1,x2,xi1,xi2)
      ddPhineg2 = -Phineg2(x1,x2,xi1,xi2)
      Gammaneg = d11*ddPhineg1*qn1(ki)
     &           + d12*(ddPhineg1*qn2(ki)+ddPhineg2*qn1(ki))
     &           + d22*ddPhineg2*qn2(ki)
      return
      end
c  The derivative D(Gamma)/D(a), for k<0
      function Gammaneg1(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /constantAzis/ Azis
      common /constantOmega / Omega
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Rcap,dRcap1,dRcap2,dRcap11,dRcap12,FK0,FK1
      pi = 4.d0*datan(1.d0)
      Rs = Rcap(x1,x2,xi1,xi2)
      Rs1 = dRcap1(x1,x2,xi1,xi2)
      Rs2 = dRcap2(x1,x2,xi1,xi2)
      Rs11 = dRcap11(x1,x2,xi1,xi2)
      Rs12 = dRcap12(x1,x2,xi1,xi2)
      x = Omega*Rs
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      FFK11 = (-FFK0-FFK1/x)*Omega*Rs1
      damn = d11*Rs1*qn1(ki) + d12*(Rs1*qn2(ki)+Rs2*qn1(ki)) +
     &       d22*Rs2*qn2(ki)
      damn1 = d11*Rs11*qn1(ki) + d12*(Rs11*qn2(ki)+Rs12*
     &        qn1(ki)) + d22*Rs12*qn2(ki)
      Gammaneg1 = -Azis*Omega/(2.d0*pi)*(FFK11*damn+FFK1*damn1)
      return
      end
c  The derivative D(Gamma)/D(b), for k<0
      function Gammaneg2(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /constantAzis/ Azis
      common /constantOmega / Omega
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Rcap,dRcap1,dRcap2,dRcap12,dRcap22,FK0,FK1
      pi = 4.d0*datan(1.d0)
      Rs = Rcap(x1,x2,xi1,xi2)
      Rs1 = dRcap1(x1,x2,xi1,xi2)
      Rs2 = dRcap2(x1,x2,xi1,xi2)
      Rs12 = dRcap12(x1,x2,xi1,xi2)
      Rs22 = dRcap22(x1,x2,xi1,xi2)
      x = Omega*Rs
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      FFK11 = (-FFK0-FFK1/x)*Omega*Rs2
      damn = d11*Rs1*qn1(ki) + d12*(Rs1*qn2(ki)+Rs2*qn1(ki)) +
     &       d22*Rs2*qn2(ki)
      damn2 = d11*Rs12*qn1(ki) + d12*(Rs12*qn2(ki)+Rs22*
     &        qn1(ki)) + d22*Rs22*qn2(ki)
      Gammaneg2 = -Azis*Omega/(2.d0*pi)*(FFK11*damn+FFK1*damn2)
      return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Gamma, for k>0
      function Gammapos(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Phipos1,Phipos2
      ddPhipos1 = -Phipos1(x1,x2,xi1,xi2)
      ddPhipos2 = -Phipos2(x1,x2,xi1,xi2)
      Gammapos = d11*ddPhipos1*qn1(ki) +
     &           d12*(ddPhipos1*qn2(ki)+ddPhipos2*
     &           qn1(ki)) + d22*ddPhipos2*qn2(ki)
      return
      end
c  The derivative D(Gamma)/D(a), for k>0
      function Gammapos1(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Phipos11,Phipos12
      ddPhi11 = -Phipos11(x1,x2,xi1,xi2)
      ddPhi12 = -Phipos12(x1,x2,xi1,xi2)
      Gammapos1 = (d11*ddPhi11+d12*ddPhi12)*qn1(ki)
     &            + (d12*ddPhi11+d22*ddPhi12)*qn2(ki)
      return
      end
c  The derivative D(Gamma)/D(b), for k>0
      function Gammapos2(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Phipos22,Phipos12
      ddPhi22 = -Phipos22(x1,x2,xi1,xi2)
      ddPhi12 = -Phipos12(x1,x2,xi1,xi2)
      Gammapos2 = (d11*ddPhi12+d12*ddPhi22)*qn1(ki)
     &            + (d12*ddPhi12+d22*ddPhi22)*qn2(ki)
      return
      end
c  ------------------------------------------------------------------------
c  The fundamental solution Gamma, for k=0
      function Gammazero(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Phizero1,Phizero2
      ddPhi1 = -Phizero1(x1,x2,xi1,xi2)
      ddPhi2 = -Phizero2(x1,x2,xi1,xi2)
      Gammazero = d11*ddPhi1*qn1(ki) + d12*(ddPhi1*qn2(ki)
     &            +ddPhi2*qn1(ki)) + d22*ddPhi2*qn2(ki)
      return
      end
c  The derivative D(Gamma)/D(a), for k=0
      function Gammazero1(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Phizero11,Phizero12
      ddPhi11 = -Phizero11(x1,x2,xi1,xi2)
      ddPhi12 = -Phizero12(x1,x2,xi1,xi2)
      Gammazero1 = (d11*ddPhi11+d12*ddPhi12)*qn1(ki)
     &             + (d12*ddPhi11+d22*ddPhi12)*qn2(ki)
      return
      end
c  The derivative D(Gamma)/D(b), for k=0
      function Gammazero2(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ d11,d12,d22
      common /normal/ qn1,qn2
      external Phizero22,Phizero12
      ddPhi22 = -Phizero22(x1,x2,xi1,xi2)
      ddPhi12 = -Phizero12(x1,x2,xi1,xi2)
      Gammazero2 = (d11*ddPhi12+d12*ddPhi22)*qn1(ki)
     &             + (d12*ddPhi12+d22*ddPhi22)*qn2(ki)
      return
      end
c  ------------------------------------------------------------------------
c  Function K_0, the Modified Bessel function of first kind and order zero
c  see Abramowitz & Stegun formulas 9.8.1, 9.8.5, 9.8.6 pages 378-379
      function FK0(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     &  0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     &  -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.0d0) then
        y=x*x/4.0d0
        FK0=(-dlog(x/2.0d0)*FI0(x))+(p1+y*(p2+y*(p3+
     &  y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        y=(2.0d0/x)
        FK0=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+
     &  y*(q4+y*(q5+y*(q6+y*q7))))))
      endif
      return
      end
c  Function I_0, the modified Bessel function of second kind and order zero
c  see Abramowitz & Stegun formula 9.8.1 page 378
      function FI0(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     &  1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     &  0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     &  0.2635537d-1,-0.1647633d-1,0.392377d-2/
      if (dabs(x).lt.3.75d0) then
        y = (x/3.75d0)**2.d0
        FI0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax = dabs(x)
        y = 3.75d0/ax
        FI0 = (dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     &        +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
      end
c  ------------------------------------------------------------------------
c  Function K_1, the modified Bessel function of first kind and order one
c  see Abramowitz & Stegun formulas 9.8.3, 9.8.7, 9.8.8 pages 378-379
      function FK1(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     &  -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     &  0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      if (x.le.2.0d0) then
        y=x*x/4.0d0
        FK1=(dlog(x/2.0d0)*FI1(x))+(1.0d0/x)*(p1+y*(p2+
     &  y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        y=2.0d0/x
        FK1=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+
     &  y*(q4+y*(q5+y*(q6+y*q7))))))
      endif
      return
      end
c  Function I_1, the modified Bessel function of second kind and order one
c  see Abramowitz & Stegun formula 9.8.3 page 378
      function FI1(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     &  0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     &  -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     &  -0.2895312d-1,0.1787654d-1,-0.420059d-2/
      if (dabs(x).lt.3.75d0) then
        y=(x/3.75d0)**2.d0
        FI1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=dabs(x)
        y=3.75d0/ax
        FI1=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     &  y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
        if (x.lt.0.0d0) FI1=-FI1
      endif
      return
      end
c  ------------------------------------------------------------------------
c  Function H_0^{(2)}, Hankel function of second kind and order zero
      complex*8 function H20(x)
      implicit real*8 (a-h,o-z)
      external FJ0, FY0
      FFJ0 = FJ0(x)
      FFY0 = FY0(x)
      H20 = cmplx(FFJ0,-FFY0)
      return
      end
c  Function J_0, Bessel function of first kind and order zero
c  see Abramowitz & Stegun formulas 9.4.1, 9.4.3 pages 369-370
      function FJ0(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     & -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     & .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     & 651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,
     & s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     & 9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(dabs(x).lt.8.d0)then
        y = x**2.d0
        FJ0 = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     &        /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax = dabs(x)
        z = 8.d0/ax
        y = z**2.d0
        xx = ax-.785398164d0
        FJ0 = sqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     &         *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end
c  Function Y_0, Bessel function of second kind and order zero
c  see Abramowitz & Stegun formulas 9.4.2, 9.4.3 pages 369-370
      function FY0(x)
      implicit real*8 (a-h,o-z)
      external FJ0
      pi = 4.0d0*datan(1.0d0)
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     & -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     & .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/-2957821389.d0,7062834065.d0,
     & -512359803.6d0,10879881.29d0,-86327.92757d0,228.4622733d0/,
     & s1,s2,s3,s4,s5,s6/40076544269.d0,745249964.8d0,
     & 7189466.438d0,47447.26470d0,226.1030244d0,1.d0/
      if (x.lt.8.d0) then
        y = x**2.d0
        FY0 = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y
     &      *(s3+y*(s4+y*(s5+y*s6)))))+.636619772d0*FJ0(x)*dlog(x)
      else
        z = 8.d0/x
        y = z**2.d0
        xx = x-.785398164d0
        FY0 = sqrt(.636619772d0/x)*(dsin(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     &        p5))))+z*dcos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end
c  ------------------------------------------------------------------------
c  Function H_1^{(2)}, Hankel function of second kind and order zero
      complex*8 function H21(x)
      implicit real*8 (a-h,o-z)
      external FJ1, FY1
      FFJ1=FJ1(x)
      FFY1=FY1(x)
      H21=cmplx(FFJ1,-FFY1)
      return
      end
c  Function J_1, Bessel function of first kind and order one
c  see Abramowitz & Stegun formulas 9.4.4, 9.4.6 page 370
      function FJ1(x)
      implicit real*8 (a-h,o-z)
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     & 242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,
     & s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     & 99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     & .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     & -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if (dabs(x).lt.8.d0) then
        y = x**2.d0
        FJ1 = x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     &        /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax = dabs(x)
        z = 8.d0/ax
        y = z**2.d0
        xx = ax-2.356194491d0
        FJ1 = sqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     &        *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     &        *sign(1.d0,x)
      endif
      return
      end
c  Function Y_1, Bessel function of second kind and order one
c  see Abramowitz & Stegun formulas 9.4.4, 9.4.5, 9.4.6 page 370
      function FY1(x)
      implicit real*8 (a-h,o-z)
      external FJ1
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     & .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     & -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      DATA r1,r2,r3,r4,r5,r6/-.4900604943d13,.1275274390d13,
     & -.5153438139d11,.7349264551d9,-.4237922726d7,.8511937935d4/,
     & s1,s2,s3,s4,s5,s6,s7/.2499580570d14,.4244419664d12,
     & .3733650367d10,.2245904002d8,.1020426050d6,.3549632885d3,1.d0/
      if (x.lt.8.d0) then
        y = x**2.d0
        FY1 = x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*
     &        (s3+y*(s4+y*(s5+y*(s6+y*s7))))))+.636619772d0
     &        *(FJ1(x)*dlog(x)-1./x)
      else
        z = 8.d0/x
        y = z**2.d0
        xx = x-2.356194491d0
        FY1 = sqrt(.636619772d0/x)*(dsin(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     &        *p5))))+z*dcos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end
c  ------------------------------------------------------------------------
c  Function R
      function Rcap(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /akar/ rhodot, rhoddot
      x1dot = x1+rhodot*x2
      xi1dot = xi1+rhodot*xi2
      x2dot = x2*rhoddot
      xi2dot = xi2*rhoddot
      Rcap = dsqrt((x1dot-xi1dot)**2.d0+(x2dot-xi2dot)**2.d0)
      return
      end
c  Derivative d(R)/d(a)
      function dRcap1(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /akar/ rhodot, rhoddot
      external Rcap
      RR=Rcap(x1,x2,xi1,xi2)
      if (dabs(x1-xi1).le.1.e-10 .and. dabs(x2-xi2).le.1.e-10) then
        dRcap1 = 0.0d0
      else
        dRcap1 = -(x1+rhodot*x2-xi1-rhodot*xi2)/RR
      endif
      return
      end
c  Derivative d(R)/d(b)
      function dRcap2(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /akar/ rhodot, rhoddot
      external Rcap
      RR = Rcap(x1,x2,xi1,xi2)
      if (dabs(x1-xi1).le.1.e-10 .and. dabs(x2-xi2).le.1.e-10 .and.
     &    dabs(rhodot).lt.1.0d0 .and. dabs(rhoddot).lt.1.0d0) then
        dRcap2 = 0.0d0
      else
        dRcap2 = -( rhodot*(x1+rhodot*x2-xi1-rhodot*xi2)+rhoddot*
     &           (x2*rhoddot-xi2*rhoddot) ) / RR
      endif
      return
      end
c  Derivative d^2(R)/d(a^2)
      function dRcap11(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      external Rcap,dRcap1
      RR = Rcap(x1,x2,xi1,xi2)
      RR1 = dRcap1(x1,x2,xi1,xi2)
      dRcap11 = -(-1.d0+RR1**2.d0)/RR
      return
      end
c  Derivative d^2(R)/d(a)d(b)
      function dRcap12(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /akar/ rhodot, rhoddot
      external Rcap,dRcap1,dRcap2
      RR = Rcap(x1,x2,xi1,xi2)
      RR1 = dRcap1(x1,x2,xi1,xi2)
      RR2 = dRcap2(x1,x2,xi1,xi2)
      dRcap12 = -(-rhodot+RR1*RR2)/RR
      return
      end
c  Derivative d^2(R)/d(b^2)
      function dRcap22(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      common /akar/ rhodot, rhoddot
      external Rcap,dRcap2
      RR = Rcap(x1,x2,xi1,xi2)
      RR2 = dRcap2(x1,x2,xi1,xi2)
      dRcap22 = -(-rhodot**2.d0-rhoddot**2.d0+RR2**2.d0)/RR
      return
      end
c  ------------------------------------------------------------------------
c  The function factorial(n)
      function fact(n)
      implicit real*8 (a-h, o-z)
      fact = 1.d0
      do k=1,n
        fact = fact*real(k)
      enddo
      return
      end
c  ------------------------------------------------------------------------
c  The function h^{1/2}(x1,x2)
      function h12(x,y)
      implicit real*8 (a-h, o-z)
      common /nequ/ nequ
      if (nequ.eq.1) then
      h12 = dcos(.1d0*x-.1d0*y) + dsin(.1d0*x-.1d0*y)
      else
      h12=dexp(.1d0*x-.1d0*y)
      endif
      return
      end
c  The function d(h^{1/2})/d(x1)
      function dh12dx1(x,y)
      implicit real*8 (a-h, o-z)
      common /nequ/ nequ
      if (nequ.eq.1) then
      dh12dx1 = .1d0 * ( dcos(.1d0*x-.1d0*y) - dsin(.1d0*x-.1d0*y) )
      else
      dh12dx1=.1d0*dexp(.1d0*x-.1d0*y)
      endif
      return
      end
c  The function d(h^{1/2})/d(x2)
      function dh12dx2(x,y)
      implicit real*8 (a-h, o-z)
      common /nequ/ nequ
      if (nequ.eq.1) then
      dh12dx2 = -.1d0 * ( dcos(.1d0*x-.1d0*y) - dsin(.1d0*x-.1d0*y) )
      else
      dh12dx2=-.1d0*dexp(.1d0*x-.1d0*y)
      endif
      return
      end
c  The exact solution cext=cext(x,y,t)
      function cext(x,y,t)
      implicit real*8 (a-h, o-z)
      common /nequ/ nequ
      external h12
      if (nequ.eq.1) then
      psiex = -.2d0*x+.13d0*y
      ft = dsin(dsqrt(t))
      cext = psiex/h12(x,y)*ft
      else
      psiex = 10.d0*(dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y))
      ft = t**2.d0
      cext = psiex/h12(x,y)*ft
      endif
      return
      end
c  The exact solution cext1=d(cext(x,y,t))/d(x)
      function cext1(x,y,t)
      implicit real*8 (a-h, o-z)
      common /nequ/ nequ
      external h12, dh12dx1
      if (nequ.eq.1) then
      psiex = -.2d0*x+.13d0*y
      psiex1 = -.2d0
      ft = dsin(dsqrt(t))
      cext1 = ( psiex1*h12(x,y) - psiex*dh12dx1(x,y) )
     &       /h12(x,y)**2.d0*ft
      else
      psiex = 10.d0*( dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y) )
      psiex1 = -2.d0*(dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y))
      ft = t**2.d0
      cext1 = ( psiex1*h12(x,y) - psiex*dh12dx1(x,y) )
     &       /h12(x,y)**2.d0*ft
      endif
      return
      end
c  The exact solution cext2=d(cext(x,y,t))/d(y)
      function cext2(x,y,t)
      implicit real*8 (a-h, o-z)
      common /nequ/ nequ
      external h12, dh12dx2
      if (nequ.eq.1) then
      psiex = -.2d0*x+.13d0*y
      psiex2 = .13d0
      ft = dsin(dsqrt(t))
      cext2 = ( psiex2*h12(x,y) - psiex*dh12dx2(x,y) )
     &       /h12(x,y)**2.d0*ft
      else
      psiex = 10.d0*(dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y))
      psiex2 = 1.d0*(dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y))
      ft = t**2.d0
      cext2 = ( psiex2*h12(x,y) - psiex*dh12dx2(x,y) )
     &       /h12(x,y)**2.d0*ft
      endif
      return
      end
c  The exact solution cexs=cexs(x,y)
      function cexs(x,y)
      implicit real*8 (a-h,o-z)
      common /laplace_constant/ s
      common /nequ/ nequ
      external h12
      if (nequ.eq.1) then
      psiex = -.2d0*x+.13d0*y
      fs = 0.8862269255d0*dexp(-.25d0/s)/s**1.5d0
      cexs = psiex/h12(x,y)*fs
      else
      psiex = 10.d0*(dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y))
      fs = 2.d0/s**3.d0
      cexs = psiex/h12(x,y)*fs
      endif
      return
      end
c  The exact solution Fexs(kp,x,y)
      function Fexs(kp,x,y)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension xn(nto),yn(nto)
      common /moduli/ d11,d12,d22
      common /normal/ xn,yn
      common /laplace_constant/ s
      common /nequ/ nequ
      external h12, Fh
      if (nequ.eq.1) then
      psiex = -.2d0*x+.13d0*y
      psiex1 = -.2d0
      psiex2 = .13d0
      fs = 0.8862269255d0*dexp(-.25d0/s)/s**1.5d0
      psiexs = psiex*fs
      psiexs1 = psiex1*fs
      psiexs2 = psiex2*fs
      Fpsi = d11*psiexs1*xn(kp) + d12*(psiexs1*yn(kp)+
     &       psiexs2*xn(kp)) + d22*psiexs2*yn(kp)
      Fexs = -Fh(kp,x,y)*psiexs + Fpsi*h12(x,y)
      else
      fs = 2.d0/(s**3.d0)
      psiex = 10.d0*( dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y) )
      psiexs = psiex*fs
      psiex1 = -2.d0*(dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y))
      psiexs1 = psiex1*fs
      psiex2 = 1.d0*(dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y))
      psiexs2 = psiex2*fs
      Fpsi = d11*psiexs1*xn(kp) + d12*(psiexs1*yn(kp)+psiexs2*xn(kp))
     &       + d22*psiexs2*yn(kp)
      Fexs = -Fh(kp,x,y)*psiexs + Fpsi*h12(x,y)
      endif
      return
      end
c  ------------------------------------------------------------------------

