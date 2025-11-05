module variables
!############################################################
!############ VARIABLES DEFINITION MODULE ###################
!############################################################
implicit none

! Convenção usada em todo o código:
! n = número de nós em X (tamanho do vetor x)
! m = número de nós em Y (tamanho do vetor y)

integer :: n, m, i, j, k, l, s, npast, Nrea, p

! Campos principais
real, allocatable :: T(:), T_ANTES(:), T_AGORA(:)
real, allocatable :: ID(:)       ! 0: saudável | 1: tumor
real, allocatable :: W(:)        ! perfusão

! Malha e geometria
real, allocatable :: x(:), y(:)
real :: xmin, xmax, ymin, ymax
real :: deltax, deltay
real :: xc, yc, rmin, rmax, eps, ax, by

! Tempo / controle
real :: dt, sim_time, npastreal, lreal
real :: tc_antes, tc_agora
integer :: t_int, t_tempi
real :: t_dec, t_tempr, tf, ti, tpro
character(3) :: rea_char

! Parâmetros (varreduras)
real, allocatable :: HZERO(:), raio_part(:), omega(:), PHI(:), raio(:), eccent(:)
real, allocatable :: timesteady(:), Tcenter(:), a_xmin(:), a_xmax(:), b_ymin(:), b_ymax(:)
integer, allocatable :: contagem(:)

! Propriedades termofísicas e constantes
real :: RHO1, RHO2, RHOB, C1, C2, CB, Q1, Q2, TB
real :: K1, K2, TZERO, T_ini
real :: muzero, MD, cor, visc, iter
real :: pi, ximag, pimag, omegaestrela, alphazero, aux

! Flags
logical :: decaicampo, tecplot_output, fileinput
logical :: geracaodireta, STEADY, CALCCAMPO, EXPLICIT

! Aux
real :: ylinha, xlinha, raiomaximo
integer :: TIU, AUXINT, auxiliarinteiro
real :: AUXREAL, auxiliarreal

end module variables
