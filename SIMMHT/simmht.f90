program simmht
  use omp_lib
  use variables
  use functions
  implicit none

  call input
  call cpu_time(ti)

  allocate(ID(m*n))
  allocate(T(m*n), T_ANTES(m*n), T_AGORA(m*n))
  allocate(W(m*n))
  allocate(x(n), y(m))
  allocate(HZERO(Nrea), raio_part(Nrea), omega(Nrea), PHI(Nrea))
  allocate(timesteady(Nrea), Tcenter(Nrea))
  allocate(contagem(Nrea))
 
xc = (xmax-xmin)/2.0 
yc = (ymax-ymin)/2.0

! Definindo elipse

by = raio*((1-eccent**2)**0.25)
ax = (raio**2)/by

b_ymax = yc + by
b_ymin = yc - by
a_xmax = xc + ax
a_xmin = xc - ax

raiomaximo = ((ymax-yc)**2.0 + (xmax-xc)**2.0)**0.5
  
  pi    = acos(-1.0)
  muzero= 4.0*pi*1.0E-07
  HZERO = 0.0; omega = 0.0; PHI = 0.0; raio_part = 0.0
  dt    = 1.0E-03
  npast = int(sim_time/dt)
  npastreal = real(npast)
  
  
print *, '================================================='
print *, '*     LABORATÓRIO DE COMPUTAÇÃO CIENTÍFICA      *' 
print *, '*         EM ESCOAMENTOS COMPLEXOS              *'
print *, '*               LCEC - UNB                      *'
print *, '================================================='
print *, '*   SIMMHT - SIMULADOR DE MAGNETOHIPERTERMIA    *'
print *, '*   Autor: Prof. Rafael Gabler Gontijo, PhD     *'
print *, '================================================='
print *, ''
  

write(*,*)'#################################################'
write(*,*)'#           LENDO DADOS DE ENTRADA              #'
write(*,*)'#################################################'
  if (fileinput) then
    call inputs
  else
    call randomica(1.5E+03,4.0E+03, HZERO,    Nrea, 9)
    call randomica(1.0E+05,3.0E+05, omega,    Nrea, 8)
    call randomica(3.0E-02,5.00E-02,PHI,      Nrea, 1)
    call randomica(5.0E-09,8.0E-09, raio_part,Nrea, 4)
    open(unit=19,file='generated_inputs.dat')
    do p=1,Nrea
      write(19,'(I12,F12.3,F12.1,F12.1,E12.2E2)') p, PHI(p), HZERO(p), omega(p), raio_part(p)
    end do
    close(19)
  end if

  write(*,*)'#################################################'
  write(*,*)'#              MONTANDO A MALHA                 #'
  write(*,*)'#################################################'
  call mesh

  write(*,*)'#################################################'
  write(*,*)'#        IDENTIFICANDO A REGIÃO TUMORAL         #'
  write(*,*)'#################################################'
  call id_vector

  call main
  call cpu_time(tf)
  tpro = tf - ti

  write(*,*) ''
  write(*,*) 'TOTAL SIMULATION TIME:', tpro, 'seconds'

  deallocate(ID,T,T_ANTES,T_AGORA,W,x,y)
  deallocate(HZERO,raio_part,omega,PHI,timesteady,Tcenter,contagem)
end program simmht
