program burgers
!implicit none

real :: a_left,a_right,q,r,s,fml,fmr,fml1,fmr1,fml2,fmr2
real :: dx,dt,time
real :: lambda1,lambda2,L,uxb,uxf
real, dimension(:),allocatable :: u0,u1,fl,fr,fm,fm1,fm2
integer :: count,tsteps,jmax,i,j,interval
character*11 ::filename
character *6 :: name
Real , parameter :: PI = 3.1415926
real,dimension(:,:),allocatable :: um



      print*,' This code computes the solution for the burger equation '
      print*,' u_t + (0.5*u^2)_x = 0 ; u(x,0) = q*arctan(s x) + r  '
      print*,' Give the time to run the model '
      read*,time
!      print*,' Give the grid size '
!      read*, dx      
       ! print*,' Give the time step '
       ! read*, dt
       print*,'Give the numerical method option 1-Laxwendoff,2-weno,3-upwind'
       read*,option

   L=2.   
   dx=0.01
   jmax = nint(L/dx)+1

   ! tsteps = nint(time/dt)+1
   tsteps=1000
   dt=time/tsteps
 
   allocate(u0(0:jmax+1))
   allocate(u1(0:jmax+1))
   allocate(fl(0:jmax+1))
   allocate(fr(0:jmax+1))
   
   dx=L/(jmax+2)
   
   do j=0,jmax+1
   u0(j)=sin(pi*L/(jmax+2)*j)+0.5
   u1(j)=u0(j)
   enddo
 
 if (option==1) then   ! Lax_wendoff
                 
!start calculation
do i=1,tsteps


           do j=1,jmax
             
              fb=(u0(j-1)**2)/2
              ff=(u0(j+1)**2)/2
              f=(u0(j)**2)/2

              uxf = (u0(j+1)+u0(j))/2
              uxb = (u0(j)+u0(j-1))/2

              u1(j)=u0(j)-(ff-fb)*(dt/(2*dx))+0.5*(dt/dx)**2*(uxf*(ff-f)-uxb*(f-fb))
             
           end do
             u1(0)=3*u1(1)-3*u1(2)+u1(3)
			 u1(jmax+1)=3*u1(jmax)-3*u1(jmax-1)+u1(jmax-2)

            do j=0,jmax+1
            u0(j)=u1(j)
            end do


 end do
   
  open(6,file="Lax_Wendoff.dat",position='APPEND')
      do j=0,jmax+1
      write (6,*) j,u1(j)
      enddo
  close(6)    

 endif

if(option==2) then  !Weno 3th Order

allocate(um(0:jmax+1,1:2))
allocate(fm(0:jmax+1))
allocate(fm1(0:jmax+1))
allocate(fm2(0:jmax+1))
do i=1,tsteps
!Runge-Kutta-TVD time propagation

call wenoflux(u0(0:jmax+1),fm,jmax)
do j=1,jmax
um(j,1)=u0(j)-dt/dx*(fm(j)-fm(j-1))
enddo
um(0,1)=4*um(1,1)-6*um(2,1)+4*um(3,1)-um(4,1)
um(jmax+1,1)=4*um(jmax,1)-6*um(jmax-1,1)+4*um(jmax-2,1)-um(jmax-3,1)



call wenoflux(um(0:jmax+1,1),fm1,jmax)
do j=1,jmax
um(j,2)=0.75*u0(j)+0.25*um(j,1)-0.25*dt/dx*(fm1(j)-fm1(j-1))
enddo
um(0,2)=4*um(1,2)-6*um(2,2)+4*um(3,2)-um(4,2)
um(jmax+1,2)=4*um(jmax,2)-6*um(jmax-1,2)+4*um(jmax-2,2)-um(jmax-3,2)



call wenoflux(um(0:jmax+1,2),fm2,jmax)
do j=1,jmax
 u1(j)=0.333*u0(j)+0.667*um(j,2)-0.667*dt/dx*(fm2(j)-fm2(j-1))
enddo
u1(0)=4*u1(1)-6*u1(2)+4*u1(3)-u1(4)
u1(jmax+1)=4*u1(jmax)-6*u1(jmax-1)+4*u1(jmax-2)-u1(jmax-3)

     do j=1,jmax
        u0(j)=u1(j)
     end do

enddo


 open(6,file="weno.dat",position='APPEND')
      do j=0,jmax+1
      write (6,*) j,u1(j)
      enddo
  close(6)  

 endif


if(option==3) then

do i=1,tsteps

do j=1,jmax
u1(j)=u0(j)-dt/dx*(u0(j)+abs(u0(j)))/2*(u0(j)-u0(j-1))-dt/dx*(u0(j)-abs(u0(j)))/2*(u0(j+1)-u0(j))
enddo

u1(0)=2*u1(1)-u1(3)
u1(jmax+1)=2*u1(jmax)-u1(jmax-1)

   do j=0,jmax+1
        u0(j)=u1(j)
     end do

enddo


 open(6,file="upwind.dat",position='APPEND')
      do j=0,jmax+1
      write (6,*) j,u1(j)
      enddo
  close(6) 

endif






end program burgers








subroutine wenoflux(u,f,jmax)

real alfa(2),s(2),w(2),v(2)
real favg(0:jmax+1),vl(0:jmax+1),vr(0:jmax+1),f(0:jmax+1),u(0:jmax+1)
integer jmax


eps=1.0d-7

do j=0,jmax+1
favg(j)=u(j)**2/2
enddo

do j=1,jmax
 
   S(1)=(favg(j+1)-favg(j))**2
   S(2)=(favg(j)-favg(j-1))**2
  
   alfa(1)=0.667/((eps+s(1))**2)
   alfa(2)=0.333/((eps+s(2))**2)
   W(1)=alfa(1)/(alfa(1)+alfa(2))
   W(2)=alfa(2)/(alfa(1)+alfa(2))
   v(1)=0.5*favg(j)+0.5*favg(j+1)
   v(2)=-0.5*favg(j-1)+1.5*favg(j)
   vl(j)=v(1)*W(1)+v(2)*W(2)

enddo


do j=0,jmax-1
   S(1)=(favg(j+2)-favg(j+1))**2
   S(2)=(favg(j+1)-favg(j))**2

   alfa(1)=0.333/((eps+s(1))**2)
   alfa(2)=0.667/((eps+s(2))**2)
   W(1)=alfa(1)/(alfa(1)+alfa(2))
   W(2)=alfa(2)/(alfa(1)+alfa(2))
   v(1)=1.5*favg(j+1)-0.5*favg(j+2)
   v(2)=0.5*favg(j)+0.5*favg(j+1)
   vr(j)=v(1)*W(1)+v(2)*W(2)
enddo

 do j=1,jmax
   if (((favg(j+1)-favg(j))/(u(j+1)-u(j)))>0) then
    f(j)=vl(j)
   else
    f(j)=vr(j)
	endif
 enddo

   f(0)=4*f(1)-6*f(2)+4*f(3)-f(4)
   f(jmax+1)=4*f(jmax)-6*f(jmax-1)+4*f(jmax-2)-f(jmax-3)
  ! f(jmax)=4*f(jmax-1)-6*f(jmax-2)+4*f(jmax-3)-f(jmax-4)

end subroutine wenoflux










 