!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Adaptive Numerical Integration with Simpson Base
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 8, 2015
!-----------------------------------------------------------------------------!


program adaptive_integration_simpson
implicit none
integer::np
real*8 ::x,f,a,b,eps,geps,S,S1,S2,h,e1,e2,eI
integer::j

!exact solution:
eI = -0.56681975015d0

geps = 1.0d-4  !error criterion

a = -1.0d0 !lower bound
b =  1.0d0 !upper bound

open(13, file="points.plt")
write(13,*)'variables ="x","f","points"'

!Adaptive integration:
S = 0.0d0
x = a
eps= geps/(b-a)

!span along x direction
200 continue
h  = (b-x)

if (x.ge.b) goto 300

!trial solution 
100 continue
S1 = h/6.0d0*(f(x)+4.0d0*f(x+0.5d0*h)+f(x+h))
S2 = h/12.0d0*(f(x)+4.0d0*f(x+0.25d0*h)+2.0d0*f(x+0.5d0*h)+4.0d0*f(x+0.75d0*h)+f(x+h))

e1 = 1.0d0/15.0d0*dabs(S2-S1)
e2 = h*eps

!error check
if (e1.le.e2) then !accept step
    write(13,*)x+0.5d0*h, f(x+0.5d0*h), 0.0d0
	S = S + (16.0d0*S2-S1)/15.0d0
	x = x + h
	goto 200    
else !reduce step size and perform trial solution
	h = 0.5d0*h
	goto 100
end if
300 continue !done

write(*,19)"exact: ", eI, 100.0d0*dabs((eI-eI)/eI)
write(*,19)"numerical: ", S, 100.0d0*dabs((S-eI)/eI)

19 format(A20,2F20.10)



! Writing given function using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="x","f"'
	do j=0,np
		x = a + dfloat(j)*(b-a)/dfloat(np)
		write(12,*) x,f(x)
	end do
close(12)
  

end

!-----------------------------------------------------------------------------!
!Given function to integrate
!-----------------------------------------------------------------------------!
real*8 function f(x)
implicit none
real*8 :: x
f = 10.0d0*dexp(-50.0d0*dabs(x)) &
  - 0.01d0/((x-0.5d0)*(x-0.5d0) + 0.001d0) &
  + 5.0d0*dsin(5.0d0*x)
end function f





