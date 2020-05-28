! 2020 - 28 - 5
! onda.f90
! Félix Cabrera (walberto.cabrera@gmail.com)

! SOLUCION A LA ECUACIÓN DE ONDA
! POR DIFERENCIAS FINITAS

! Codificación del texto: UTF-8
! Compiladores probados: GNU Fortran (Ubuntu 9.2.1-9ubuntu2) 9.2.1 2019008

! Algoritomo tomado del libro de Analisis Numerico de Burden
! Algoritmo 12.1 pgs 566

! Rquiere: Fx.f90, Gx.f90, funciones.f90
! Instrucciones de compilacion:
! gfortran -Wall -pedantic -std=f95 -c FX.f90
! gfortran -Wall -pedantic -std=f95 -c GX.f90
! gfortran -Wall -pedantic -std=f95 -c funciones.f90
! gfortran -Wall -pedantic -std=f95 -c onda.f90
! gfortran -Wall -pedantic -std=f95 -o onda onda.o funciones.o FX.o GX.o
! ./onda


PROGRAM onda
  USE funciones
  IMPLICIT NONE
  
  ! Definimos variables principales
  INTEGER(8):: m,N
  REAL(8):: T,l,alpha
  REAL(8), ALLOCATABLE :: w(:,:)
  ! Definifmos variables auxiliares
  INTEGER(8) :: i,j
  REAL(8):: lambda,h,k,x
  ! Variables de control
  INTEGER(4) :: err
  
    ! Abrimos el archivo de configuración (12)
  OPEN (12, FILE='parametros.config', STATUS='old', IOSTAT=err)
  IF (err .ne. 0) STOP 'parametros.config is missing'

  READ(12,*) l ! limite de  x
  READ(12,*) T ! limite de tiempo
  READ(12,*) alpha ! constante
  READ(12,*) m ! cantidad de puntos
  READ(12,*) N ! pasos
  CLOSE(12)
  
  ALLOCATE(w(m+1,N+1))
  
    ! Abrimos el acribo para los resultados
  OPEN(14, FILE='Resultados.dat', STATUS='new', IOSTAT=err)
  IF (err .ne. 0) STOP 'Resultados.dat already exists'

  PRINT *, '***********************************************************************'
  PRINT *, 'Resolviendo la ecuación de onda'
  PRINT *, 'Método= differencias finitas'
  PRINT *, 'Pasos =',N
  PRINT *, 'Procesando...'
  
  
  ! Paso 1
  h= l/REAL(m,8)
  k= T/REAL(N,8)
  lambda= (k*alpha)/h
  
  ! Paso 2
  w = 0
  
  ! Paso 3
  w(1,1) = Fx(REAL(0.,8))
  w(m+1,1) = Fx(l)
  
  ! Paso 4
  DO i=1, m-1
    w(i+1,1)= Fx(i*h)
    w(i+1,2)= (1-(lambda**2))*Fx(i*h)+((lambda**2)/2)*(Fx((i+1)*h)+Fx((i-1)*h))+k*Gx(i*h)
    IF (MOD(i,((m-1)/5))==0) THEN
      WRITE(*, fmt='(1x,a,i0)', advance='no') '%'
    END IF
  END DO
  
  ! Paso 5
  DO j=1, N-1
    DO i=1, m-1
      w(i+1,j+2)=2*(1-(lambda**2))*w(i+1,j+1)+(lambda**2)*(w(i+2,j+1)+w(i,j+1))-w(i+1,j)
    END DO
    
    IF (MOD(j,((N-1)/5))==0) THEN
      WRITE(*, fmt='(1x,a,i0)', advance='no') '%'
    END IF
  END DO
  
  ! Paso 6
  DO j=1,N+1
    t=(j-1)*k
    DO i=1,m+1
      x = (i-1)*h
      WRITE(14,*) x,t,w(i,j)
    END DO
    IF (MOD(j,((N+1)/10))==0) THEN
      WRITE(*, fmt='(1x,a,i0)', advance='no') '%'
    END IF
  END DO
  
  ! Paso 7 FIN
  DEALLOCATE(w)
  PRINT *, 'Fin'
  PRINT *, '***********************************************************************' 
END PROGRAM onda












