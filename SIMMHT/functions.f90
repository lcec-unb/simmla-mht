module functions
  use omp_lib
  use variables
  implicit none
contains

!---------------- Índice do centro (nos índices de nó) ----------------
pure integer function center_idx() result(idx)
  integer :: jc, kc
  jc = (n+1)/2
  kc = (m+1)/2
  idx = jc + (kc-1)*n
end function center_idx

!---------------- Pós-processamento dos campos ----------------
! Calcula estatísticas no tumor e métrica de "área de ação" no tecido saudável
subroutine compute_metrics(Tmin_tumor, Tavg_tumor, Tstd_tumor, SAR_mean, R_afet)
  implicit none
  real, intent(out) :: Tmin_tumor, Tavg_tumor, Tstd_tumor
  real, intent(out) :: SAR_mean, R_afet

  integer :: j, k, idx
  real :: sumT, sumT2, cnt
  real :: dx, dy, r, thr, Tval
  real :: SAR_local

  ! --- estatísticas no tumor (ID==1) ---
  Tmin_tumor =  1.0e30
  sumT = 0.0
  sumT2= 0.0
  cnt  = 0.0

  do k=1,m
    do j=1,n
      idx = j + (k-1)*n
      if (ID(idx) .eq. 1.0) then
        Tval = T(idx)
        if (Tval < Tmin_tumor) Tmin_tumor = Tval
        sumT  = sumT  + Tval
        sumT2 = sumT2 + Tval*Tval
        cnt   = cnt   + 1.0
        end if
    end do
  end do

  if (cnt > 0.0) then
    Tavg_tumor = sumT / cnt
    Tstd_tumor = sqrt( max(0.0, sumT2/cnt - Tavg_tumor*Tavg_tumor) )
  else
    Tavg_tumor = 0.0
    Tstd_tumor = 0.0
    Tmin_tumor = 0.0
  end if

  ! --- SAR médio no tumor (W/kg) ---
  ! Consistente com o termo de geração que você usa no explícito:
  ! GERACAO = muzero*pi*PHI*MD*HZERO*omega*(-ximag)*cor   (W/m^3)
  ! SAR = GERACAO / RHO2                                  (W/kg)
  SAR_local = cor * muzero * acos(-1.0) * PHI(p) * MD * HZERO(p) * omega(p) * (-ximag)   ! W/m^3
  if (RHO2 > 0.0) then
    SAR_mean = SAR_local / RHO2
  else
    SAR_mean = 0.0
  end if

  ! --- Raio/Área de "afetacao" no tecido saudável (ID==0) ---
  ! Critério: T > T_ini + 1.0 [°C]
  R_afet = 0.0
  thr = T_ini + 1.0

  do k=1,m
    do j=1,n
      idx = j + (k-1)*n
      if (ID(idx) .eq. 0.0) then
        if (T(idx) > thr) then
          dx = x(j) - xc
          dy = y(k) - yc
          r  = sqrt(dx*dx + dy*dy)
          if (r > R_afet) R_afet = r
        end if
      end if
    end do
  end do

 end subroutine compute_metrics


!-------------------------------------------------------------
! LÊ entradas por simulação (PHI, H, omega, raio_part)
!-------------------------------------------------------------
subroutine inputs
  open(15,file='inputs.dat')
  do i=1,Nrea
    read(15,*) iter, PHI(i), HZERO(i), omega(i), raio_part(i), raio(i), eccent(i)
  end do
  close(15)
end subroutine inputs

!-------------------------------------------------------------
! MALHA
!-------------------------------------------------------------
subroutine mesh
  deltax = (xmax - xmin) / real(n-1)
  deltay = (ymax - ymin) / real(m-1)

  x(1) = xmin
  do j=2,n; x(j) = x(j-1) + deltax; end do

  y(1) = ymin
  do k=2,m; y(k) = y(k-1) + deltay; end do
end subroutine mesh

!-------------------------------------------------------------
! ID_VECTOR (0 saudável / 1 tumor)
! índice linear: idx = j + (k-1)*n
!-------------------------------------------------------------
subroutine id_vector(p)
  implicit none
  integer :: j, k, idx
  real :: coef, radx2, arg, denom
  integer, intent(in) :: p
  
  ! zera ID
  do k=1,m
    do j=1,n
      idx     = j + (k-1)*n
      ID(idx) = 0.0
    end do
  end do
    
  ! coeficiente geométrico (evita divisão por zero)
  denom = a_xmax(p) - xc
  if (abs(denom) < 1.0e-12) denom = sign(1.0e-12, denom)
  coef  = (b_ymax(p) - yc) / denom

  ! (a_xmax - xc)^2, usado no argumento da raiz
  radx2 = (a_xmax(p) - xc)**2

  do k=1,m
    do j=1,n
      idx = j + (k-1)*n

      if (x(j) .ge. a_xmin(p) .and. x(j) .le. a_xmax(p)) then
        ! argumento da raiz: garante não-negatividade
        arg = radx2 - (x(j) - xc)**2
        if (arg < 0.0) arg = 0.0

        ! Quadrantes superiores: yc <= y <= b_ymax
        if (y(k) .ge. yc .and. y(k) .le. b_ymax(p)) then
          if ( (y(k) - yc) .le. coef * sqrt(arg) ) ID(idx) = 1.0
        ! Quadrantes inferiores: b_ymin <= y <= yc
        else if (y(k) .le. yc .and. y(k) .ge. b_ymin(p)) then
          if ( (yc - y(k)) .le. coef * sqrt(arg) ) ID(idx) = 1.0
        end if
      end if

    end do
  end do
end subroutine id_vector


!-------------------------------------------------------------
! PERFUSÃO SANGUÍNEA W(T) (corrigido)
!-------------------------------------------------------------
subroutine perfusao_sanguinea
  integer :: j,k,idx

  do k=2,m-1
    do j=2,n-1
      idx = j + (k-1)*n
      if (ID(idx) .eq. 0.0) then
        if      (T(idx) .lt. 37.0) then
          W(idx) = 8.33E-04
        else if (T(idx) .le. 42.0) then
          W(idx) = 8.33E-04 - ((T(idx) - 37.0)**4.8)/(5.438E+06)
        else
          W(idx) = 4.16E-04
        end if
      else
        if (T(idx) .lt. 45.0) then
          W(idx) = 4.5E-04 + 3.55E-03*EXP(-((T(idx)-45.0)**2.0)/12.0)
        else
          W(idx) = 4.00E-03
        end if
      end if
    end do
  end do

  ! Fronteiras: saudável padrão
  do j=1,n
    W(j) = 8.33E-04                      ! k=1
    W((m-1)*n + j) = 8.33E-04            ! k=m
  end do
  do k=2,m-1
    W(1 + (k-1)*n) = 8.33E-04            ! j=1
    W(n + (k-1)*n) = 8.33E-04            ! j=n
  end do
end subroutine perfusao_sanguinea

!-------------------------------------------------------------
! CONTORNOS NEUMANN (grad=0) para o EXPLÍCITO
!-------------------------------------------------------------
subroutine aplicar_contornos_explicito(Tnew, Told)
  real, intent(inout) :: Tnew(:)
  real, intent(in)    :: Told(:)
  integer :: j,k,idx

  ! lower (k=1)
  k = 1
  do j=1,n
    idx = j + (k-1)*n
    Tnew(idx) = 2.0*Told(j + (k  )*n) - Told(j + (k+1)*n)
  end do
  ! upper (k=m)
  k = m
  do j=1,n
    idx = j + (k-1)*n
    Tnew(idx) = 2.0*Told(j + (k-2)*n) - Told(j + (k-3)*n)
  end do
  ! left (j=1)
  do k=2,m-1
    idx = 1 + (k-1)*n
    Tnew(idx) = 2.0*Told(2 + (k-1)*n) - Told(3 + (k-1)*n)
  end do
  ! right (j=n)
  do k=2,m-1
    idx = n + (k-1)*n
    Tnew(idx) = 2.0*Told(n-1 + (k-1)*n) - Told(n-2 + (k-1)*n)
  end do
end subroutine aplicar_contornos_explicito

!-------------------------------------------------------------
! ESQUEMA EXPLÍCITO (Told/Tnew + pronto para OpenMP)
!-------------------------------------------------------------
subroutine solucao_explicita
  integer :: j,k, idx, idxE, idxW, idxN, idxS
  real :: DIFUSAO_local, PERFUSAO_local, GERACAO_local
  real, allocatable :: Told(:), Tnew(:)

  allocate(Told(n*m), Tnew(n*m))
  Told = T

!$omp parallel do collapse(2) private(j,k,idx,idxE,idxW,idxN,idxS,DIFUSAO_local,PERFUSAO_local,GERACAO_local) shared(Told,Tnew)
  do k=2,m-1
    do j=2,n-1
      idx  =  j    + (k-1)*n
      idxE = (j+1) + (k-1)*n
      idxW = (j-1) + (k-1)*n
      idxN =  j    +  k    *n
      idxS =  j    + (k-2)*n

      if (ID(idx) .eq. 0.0) then
        DIFUSAO_local = (K1/deltax**2)*(Told(idxE)+Told(idxW)-2.0*Told(idx)) &
                      + (K1/deltay**2)*(Told(idxN)+Told(idxS)-2.0*Told(idx))
        PERFUSAO_local = RHOB*CB*W(idx)*(TB - Told(idx))
        Tnew(idx) = Told(idx) + (dt/(RHO1*C1))*(DIFUSAO_local + PERFUSAO_local + Q1)
      else
        DIFUSAO_local = (K2/deltax**2)*(Told(idxE)+Told(idxW)-2.0*Told(idx)) &
                      + (K2/deltay**2)*(Told(idxN)+Told(idxS)-2.0*Told(idx))
        PERFUSAO_local = RHOB*CB*W(idx)*(TB - Told(idx))
        GERACAO_local  = muzero*acos(-1.0)*PHI(p)*MD*HZERO(p)*omega(p)*(-ximag)
        Tnew(idx) = Told(idx) + (dt/(RHO2*C2))*(DIFUSAO_local + PERFUSAO_local + Q2 + cor*GERACAO_local)
      end if
    end do
  end do
!$omp end parallel do

  call aplicar_contornos_explicito(Tnew, Told)
  T = Tnew
  deallocate(Told, Tnew)
end subroutine solucao_explicita

!-------------------------------------------------------------
! GS/SOR (dense) – substitui Jacobi (para usos pontuais)
!-------------------------------------------------------------
subroutine gauss_seidel_dense(A, x, b, n, maxit, tol, omegaSOR)
  integer, intent(in) :: n, maxit
  real,    intent(in) :: tol, omegaSOR
  real,    intent(in) :: A(n,n), b(n)
  real,    intent(inout) :: x(n)
  integer :: it, i, j
  real :: sigma, dx, resmax, xnew

  do it=1,maxit
    resmax = 0.0
    do i=1,n
      sigma = 0.0
      do j=1,i-1; sigma = sigma + A(i,j)*x(j); end do
      do j=i+1,n; sigma = sigma + A(i,j)*x(j); end do
      xnew = (b(i) - sigma)/A(i,i)
      dx   = xnew - x(i)
      x(i) = x(i) + omegaSOR*dx
      resmax = max(resmax, abs(dx))
    end do
    if (resmax < tol) exit
  end do
end subroutine gauss_seidel_dense


subroutine tecplot(a)
integer a
character(3) int_char
write(int_char, '(I3)') a
open(a, file="fields"//int_char//".plt")
write(a,*) 'Variables="x","y","ID","W","T"'
write(a,*) 'ZONE F=POINT,I='
write(a,*) n
write(a,*) ',J='
write(a,*) m
do i=1,n
do j=1,m
write(a,'(F12.4,F12.4,F12.4,F12.4)') x(i),y(j),ID(i+((j-1)*n)),W(i+((j-1)*n)),T(i+((j-1)*n))
end do
end do
end subroutine tecplot
!-------------------------------------------------------------
! ROTINA MAIN (varreduras)
!-------------------------------------------------------------
subroutine main
  use variables, only: tecplot_output
  integer :: AUXINTloc, auxiliarinteiro, uout, ios
  real :: AUXREALloc, auxiliarreal
  real :: tmin_t, tmax_t, tavg_t, tstd_t
  real :: sar_mean, r_afet, a_afet
  integer :: ic
  525 FORMAT(F12.4,F12.4,F12.4,F12.4)

  T = T_ini
  lreal = 0.0

  write(*,*)'#################################################'
  write(*,*)'#           EXECUTANDO AS SIMULAÇÕES            #'
  write(*,*)'#################################################'


open(newunit=uout, file='outputs.dat', status='replace', action='write', iostat=ios)
if (ios /= 0) then
  write(*,*) 'ERROR: failed to open outputs.dat, iostat=', ios
  return   
end if

write(uout,*) 'Variables="SIMULATION","PHI","H","W","PARTICLE_RADIUS",'// &
              '"TC","S_TIME",'// &
              '"TMIN","TAVG","TSTD",'// &
              '"SAR","R_AFET"'
call flush(uout)

  do p=1,Nrea
  
  write(*,*)'#################################################'
  write(*,*)'#        IDENTIFICANDO A REGIÃO TUMORAL         #'
  write(*,*)'#################################################'
  
  call id_vector(p)
  
    contagem(p)=0
    call imag
    write(*,'(A,1x)') 'PROGRESSO DA SIMULAÇÃO:',p

    do l=1,npast
      call perfusao_sanguinea
      T_ANTES = T
      tc_antes = T((m*n - m)/2)

      call solucao_explicita

      tc_agora = T((m*n - m)/2)
      T_AGORA  = T

      lreal = l
write(*,'(A,F12.4," s   ",F12.2," °C")',advance='no') achar(27)//'[2K'//achar(13), l*dt, T((m*n-m)/2)

      if (maxval(abs(T_ANTES - T_AGORA)) .lt. 1.0E-06) then
        if (contagem(p) .lt. 1) then
          timesteady(p) = l*dt
          contagem(p)   = contagem(p) + 1
          
        ! índice do centro e métricas
    	ic = center_idx()
    	call compute_metrics(tmin_t, tavg_t, tstd_t, sar_mean, r_afet)

write(uout,'(I5,F8.3,F9.1,F12.1,E10.2E2,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2,E10.2E2)') p, PHI(p), HZERO(p), omega(p), raio_part(p), &
                     T(ic), timesteady(p), &
                     tmin_t, tavg_t, tstd_t, &
                     sar_mean, r_afet     
          call flush(uout)
          if(tecplot_output) then
          call tecplot(p)
          end if
          T = TZERO
          exit
        end if
      end if

    end do
    write(*,*)
    write(*,'(A,1x)') 'FIM DA SIMULAÇÃO', p
  end do

!  close(2112)
  close(uout)
end subroutine main

!-------------------------------------------------------------
! ROTINA QUE CALCULA X'' (mantida, sem alterações conceituais)
!-------------------------------------------------------------
subroutine imag
  implicit none
  double precision :: a1, b1, h1
  integer :: M1,i1,n1,j1
  double precision, allocatable :: y1(:),t1(:)
  double precision :: k11,k22,k33,k44
  double precision :: omega1, lambda1, alpha1

  omega1 = (2.0d0*pi*omega(p)*6.0d0*visc*pi*raio_part(p)**3) / ((1.38E-23)*(T((m*n-m)/2)+273))
  alpha1 = (4.0d0*pi*(raio_part(p)**3)*muzero*MD*HZERO(p)) / (3.0d0*(1.38E-23)*(T((m*n-m)/2)+273))
  lambda1= (muzero*(MD**2)*pi*raio_part(p)**3) / (18.0d0*(1.38E-23)*(T((m*n-m)/2)+273))

  pi = dacos(-1.d0)
  a1 = 0.0d0
  b1 = (2.0d0*pi/omega1)*100
  M1 = 10000.d0
  h1 = (b1-a1)/(1.0d0*M1)

  allocate(y1(0:M1), t1(0:M1))
  y1(0)=0.0000001d0
  t1(0)=a1
  do i1=1,M1
    t1(i1)=a1+i1*h1
  end do

  do i1=0,M1-1
    k11 = h1*f1(t1(i1),y1(i1))
    k22 = h1*f1(t1(i1)+0.5d0*h1,y1(i1)+0.5d0*k11)
    k33 = h1*f1(t1(i1)+0.5d0*h1,y1(i1)+0.5d0*k22)
    k44 = h1*f1(t1(i1)+h1,y1(i1)+k33)
    y1(i1+1) = y1(i1)+(1.0d0/6.0d0)*(k11+2.0d0*k22+2.0d0*k33+k44)
  end do

  ximag = 0.0
  do j1=1,M1-1
    h1 = t1(j1+1)-t1(j1)
    ximag = ximag + (h1/2.0d0)*(omega1/(2.0d0*pi*100))*(((1.0d0/dtanh(y1(j1+1))-1.0d0/y1(j1+1)+ &
      8.0d0*PHI(p)*lambda1*(1.0d0/(dtanh(y1(j1+1)))-1.0d0/(y1(j1+1)))*((1.0d0/y1(j1+1))**2-(1.0d0/dsinh(y1(j1+1)))**2)))*dcos(omega1*t1(j1+1)) + &
      ((1.0d0/dtanh(y1(j1))-1.0d0/y1(j1)+8.0d0*PHI(p)*lambda1*(1.0d0/(dtanh(y1(j1)))-1.0d0/(y1(j1)))*((1.0d0/y1(j1))**2- &
      (1.0d0/dsinh(y1(j1)))**2)))*dcos(omega1*t1(j1)) )
  end do

contains
  function f1(var1,var2) result(res)
    implicit none
    double precision :: res, var1, var2
    res = ((-2.0d0)*(0.75/var2)*(+1.0*var2-alpha1*dsin(omega1*var1))*(1.0d0/DTANH(var2)-1.0d0/var2)) / &
          ((1.00d0)*(1.d0/var2**2 - 1.0d0/DSINH(var2)**2 + 8.0d0*PHI(p)*lambda1*((1.0d0/var2**2-1.0d0/DSINH(var2)**2)**2 + &
           (1.0d0/DSINH(var2)-1.0d0/var2)*(-2.0d0/var2**3+2.0d0*(1.0d0/DTANH(var2))*(1.0d0/DSINH(var2)**2)))))
  end function f1
end subroutine imag

!-------------------------------------------------------------
! GERAÇÃO ALEATÓRIA (mantida)
!-------------------------------------------------------------
subroutine randomica(a,b,c,n,d)
  real a,b       ! range
  integer n, m
  real c(n)
  integer d,i,e
  integer f(8)
  integer, allocatable :: seed(:)

  call random_seed(size=m)
  allocate(seed(m))
  call DATE_AND_TIME(values=f)
  call SYSTEM_CLOCK(count=e)

  do i=1,m
    seed(i) = 47*d + f(8)*i*d*12000 + e*(3*d+i)
  end do
  call random_seed(put=seed)

  call random_number(c)
  c = a + (b-a)*c
end subroutine randomica

end module functions
