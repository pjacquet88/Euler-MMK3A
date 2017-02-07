PROGRAM main
  USE module
  IMPLICIT NONE
  INTEGER::i,k,a,k1
  REAL*8::time
  REAL*8,DIMENSION(:,:),ALLOCATABLE::U1,U2
  REAL*8,DIMENSION(:),ALLOCATABLE::Entropie


  CALL Initialisation()
  ALLOCATE(U1(3,N),U2(3,N),Entropie(N))
  CALL Condition_Initiale()
  CALL DT_CFL(U0)
  CALL Print_Densite(U0,0)
  CALL Print_Vitesse(U0,0)
  CALL Print_Pression(U0,0)
  U1=U0
  k=0
  k1=0

  OPEN(UNIT=42,file='EntropieRK.dat',action='write')
  DO WHILE (time<=T)
    time=time+dt
    k=k+1
    CALL Time_Step(U1,U2)
    U1=U2

    CALL DT_CFL(U1)
    
    if (modulo(k,100)==0) then
    k1=k1+1
    print*,'TEST'
    CALL Print_Densite(U1,k1)
    ! CALL Print_Vitesse(U1,k)
    ! CALL Print_Pression(U1,k)
    end if
    PRINT*,'Etape :',k,'TEMPS :',time

    ! IMPRESSION DE L'ENTROPIE
    DO i=1,N
      Entropie(i)=8.31/(gamma-1.0)*log(Pression(U1(:,i)/(U1(1,i)**gamma)))
    END DO
    WRITE(42,*),time,maxval(Entropie)
  END DO

  CLOSE(42)
  a=FLOOR(1.0*k1/100)+1

  OPEN(unit=78,file='script.gnuplot',action='write')
  WRITE(78,*),'load "trace1.gnuplot"'
  WRITE(78,*),'n=',k1
  WRITE(78,*),'a=',a
  WRITE(78,*),'load "trace2.gnuplot"'
  CLOSE(78)

  CALL SYSTEM('gnuplot5 script.gnuplot')

   CALL Print_Densite(U1,17)
  ! CALL Print_Vitesse(U1,1)
  ! CALL Print_Pression(U1,1)

END PROGRAM main
