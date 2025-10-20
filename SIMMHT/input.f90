subroutine input
  use variables
  implicit none
  character(len=256) :: texto   
  505   FORMAT(1X,A54,1X,E11.4E2)
  507   FORMAT(1X,A54,1X,I6)
  508   FORMAT(1X,A54,L10)

  open(2,file='config.dat')
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,505) texto,xmin
    READ (2,505) texto,xmax
    READ (2,505) texto,ymin
    READ (2,505) texto,ymax
    ! IMPORTANTE: no arquivo, a ordem é X depois Y.
    ! Aqui alinhamos a convenção do código: n = nós em X, m = nós em Y.
    READ (2,507) texto,n   ! NODES QUANTITY IN X DIRECTION
    READ (2,507) texto,m   ! NODES QUANTITY IN Y DIRECTION
    READ (2,508) texto,fileinput
    READ (2,507) texto,Nrea
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,505) texto,raio
    READ (2,505) texto,RHO1
    READ (2,505) texto,RHO2
    READ (2,505) texto,K1
    READ (2,505) texto,K2
    READ (2,505) texto,Q1
    READ (2,505) texto,Q2
    READ (2,505) texto,C1
    READ (2,505) texto,C2
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,505) texto,TZERO
    READ (2,505) texto,T_ini
    READ (2,505) texto,MD
    READ (2,505) texto,cor
    READ (2,505) texto,visc
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,505) texto,sim_time
    READ (2,505) texto,t_dec
    READ (2,505) texto,t_tempr
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,'(A)') texto
    READ (2,505) texto,RHOB
    READ (2,505) texto,CB
    READ (2,505) texto,TB
  close(2)

  t_int   = t_dec
  t_tempi = t_tempr
end subroutine input
