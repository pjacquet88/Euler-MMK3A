MODULE module
  USE limiteurs
  IMPLICIT NONE


  REAL*8,PUBLIC,PARAMETER  ::  gamma = 1.4d0
  INTEGER,PUBLIC :: id_ordre,id_temps,id_problem,N
  REAL*8::L,T,dx,dt,dt_min
  REAL*8,PARAMETER :: cfl=0.5
  REAL*8,DIMENSION(3)::U_Left,U_Right
  REAL*8,DIMENSION(:,:),ALLOCATABLE::U0

  !PUBLIC  :: F,HLL
  !PRIVATE ::


CONTAINS


  FUNCTION Pression(U)
    IMPLICIT NONE
    REAL*8,DIMENSION(3),INTENT(IN)::U
    REAL*8::Pression
    REAL*8::e
    e=U(3)/U(1)-(((U(2)/U(1))**2)/2)
    Pression=U(1)*(gamma-1)*e
  END FUNCTION

  FUNCTION RhoE(rho,u,p)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::rho,u,p
    REAL*8::e
    REAL*8::RhoE
    e=(p/(gamma-1))/rho
    RhoE=e+(u**2)/2
    RhoE=rho*RhoE
  END FUNCTION

  FUNCTION F(U)
    IMPLICIT NONE
    REAL*8,DIMENSION(3),INTENT(IN)::U
    REAL*8,DIMENSION(3)::F
    F(1)=U(2)
    F(2)=U(2)**2/U(1)+Pression(U)
    F(3)=(U(2)/U(1))*(U(3)+Pression(U))
  END FUNCTION F

  FUNCTION Flux_HLL(U1,U2)
    IMPLICIT NONE
    REAL*8,DIMENSION(3),INTENT(IN)::U1,U2
    REAL*8,DIMENSION(3):: Flux_HLL
    REAL*8::vp_min,vp_max,c1,c2
    ! REAL*8::c3,c4,c5,c6
    c1=sqrt(gamma*Pression(U1)/U1(1))
    c2=sqrt(gamma*Pression(U2)/U2(1))

    vp_min=min(U1(2)/U1(1)-c1,U2(2)/U2(1)-c2)
    vp_max=max(U1(2)/U1(1)+c1,U2(2)/U2(1)+c2)


    IF(vp_min>0d0) THEN
      Flux_HLL=F(U1)
    ELSEIF (vp_max<0d0) THEN
      Flux_HLL=F(U2)
    ELSEIF((vp_min<=0.0d0).AND.(vp_max>=0.0d0)) THEN
      Flux_HLL=((vp_max*F(U1)-vp_min*F(U2))+vp_min*vp_max*(U2-U1))/(vp_max-vp_min)
    END IF

  END FUNCTION Flux_HLL

  FUNCTION Flux_MUSCL(U0,U1,U2,U3)
    IMPLICIT NONE
    REAL*8,DIMENSION(3),INTENT(IN)::U0,U1,U2,U3
    REAL*8,DIMENSION(3)::Flux_MUSCL
    REAL*8,DIMENSION(3)::UG,UD
    REAL*8::rd1,rd2,rd3,rg1,rg2,rg3
    ! REAL*8::rd1n,rd2n,rd3n,rg1n,rg2n,rg3n
    ! REAL*8::rd1d,rd2d,rd3d,rg1d,rg2d,rg3d
    REAL*8,PARAMETER::eps=1d-16

    rg1=(U1(1)-U0(1)+eps)/(U2(1)-U1(1)+eps)
    rg2=(U1(2)-U0(2)+eps)/(U2(2)-U1(2)+eps)
    rg3=(U1(3)-U0(3)+eps)/(U2(3)-U1(3)+eps)

    rd1=(U2(1)-U1(1)+eps)/(U3(1)-U2(1)+eps)
    rd2=(U2(2)-U1(2)+eps)/(U3(2)-U2(2)+eps)
    rd3=(U2(3)-U1(3)+eps)/(U3(3)-U2(3)+eps)

    UG(1)=U1(1)+0.5*limiteur(rg1)*(U2(1)-U1(1))
    UG(2)=U1(2)+0.5*limiteur(rg2)*(U2(2)-U1(2))
    UG(3)=U1(3)+0.5*limiteur(rg3)*(U2(3)-U1(3))

    UD(1)=U2(1)-0.5*limiteur(rd1)*(U3(1)-U2(1))
    UD(2)=U2(2)-0.5*limiteur(rd2)*(U3(2)-U2(2))
    UD(3)=U2(3)-0.5*limiteur(rd3)*(U3(3)-U2(3))

    Flux_MUSCL=Flux_HLL(UG,UD)
  END FUNCTION


  SUBROUTINE Initialisation()
    IMPLICIT NONE
    !INTEGER,INTENT(IN) ::  temps,ordre,limiteur
    id_temps=1
    id_ordre=1
    id_limiteur=0
    N=10

    PRINT*,'Nombre de maille pour la discrétisation 1D : '
    READ*, N

    PRINT*, 'Temps total de l experience :'
    READ*, T

    PRINT*,'Pas de temps souhaité :'
    READ*, dt_min

    ALLOCATE(U0(3,N))

    PRINT*,'Quel problème étudier ?'
    PRINT*,'1) = Tube à choc de Sod'
    PRINT*,'2) = 1-détente'
    PRINT*,'3) = Blast Wave'
    READ*, id_problem

    PRINT*,'Quel schéma en temps utiliser ?'
    PRINT*,'1) = EULER  (ordre 1)'
    PRINT*,'2) = RUNGE KUTTA 2 (ordre 2)'
    PRINT*,'3) = Ralston (ordre2)'
    READ*, id_temps

    PRINT*,'A quel ordre approcher les flux numériques ?'
    PRINT*,'1) = Ordre 1 (Méthode HLL)'
    PRINT*,'2) = Ordre 2 (Méthode HLL + MUSCL)'
    READ*, id_ordre

    IF (id_ordre==2) THEN
      PRINT*,'Quel limiteur utiliser pour la méthode MSUCL ?'
      PRINT*,'0) = Sans limiteur'
      PRINT*,'1) = Superbee'
      PRINT*,'2) = Van Leer'
      PRINT*,'3) = Minmod'
      PRINT*,'4) = Ospre'
      READ*,id_limiteur
    END IF

    L=1.0d0
    dx=L/N
  END SUBROUTINE

  SUBROUTINE Condition_Initiale()
    IMPLICIT NONE
    INTEGER::i
    REAL*8,DIMENSION(3,N)::U0bis

    IF (id_problem==1) THEN
      U_Left(1)=1.0d0
      U_Left(2)=0.0d0*U_Left(1)
      U_Left(3)=RhoE(1.0d0,0.0d0,1.0d0)

      U_Right(1)=0.125d0
      U_Right(2)=0.0d0*U_Right(1)
      U_Right(3)=RhoE(0.125d0,0.0d0,0.1d0)

      DO i=1,floor(N/2.0)
        U0(:,i)=U_Left
      END DO
      DO i=floor(N/2.0)+1,N
        U0(:,i)=U_Right
      END DO

    ELSE IF (id_problem==2) THEN
      U_Left(1)=1.0d0
      U_Left(2)=-1.0d0*U_Left(1)
      U_Left(3)=RhoE(1.0d0,-1.0d0,1.5d0)

      U_Right(1)=0.1989d0
      U_Right(2)=1.0d0*U_Right(1)
      U_Right(3)=RhoE(0.1989d0,1.0d0,0.1564d0)

      DO i=1,floor(N/2.0)
        U0(:,i)=U_Left
      END DO
      DO i=floor(N/2.0)+1,N
        U0(:,i)=U_Right
      END DO

    ELSE IF (id_problem==3) THEN
      U_Left(1)=1.0d0
      U_Left(2)=0.0d0*U_Left(1)
      U_Left(3)=RhoE(1.0d0,0.0d0,1000.0d0)

      U_Right(1)=1.0d0
      U_Right(2)=0.0d0*U_Left(1)
      U_Right(3)=RhoE(1.0d0,0.0d0,100.0d0)

      U0(1,:)=1.0d0
      U0(2,:)=0.0d0
      U0(3,:)=RhoE(1.0d0,0.0d0,0.01d0)
      DO i=1,floor(N/10.0)
        U0(:,i)=U_Left
      END DO
      DO i=floor(9*N/10.0),N
        U0(:,i)=U_Right
      END DO
    END IF

    ! U0bis=U0
    ! DO i=1,N
    !   U0(:,i)=U0bis(:,N+1-i)
    ! END DO

  END SUBROUTINE


  SUBROUTINE Time_Step(U1,U2)
    IMPLICIT NONE
    REAL*8,DIMENSION(3,N),INTENT(IN)::U1
    REAL*8,DIMENSION(3,N),INTENT(OUT)::U2
    REAL*8,DIMENSION(3,N)::U
    INTEGER::i
    REAL*8,DIMENSION(3,N)::k1,k2
    IF (id_ordre==1) THEN
      IF (id_temps==1) THEN
        U2(:,1)=U1(:,1)-Space_Term(U1(:,1),U1(:,1),U1(:,2))
        DO i=2,N-1
          U2(:,i)=U1(:,i)-Space_Term(U1(:,i-1),U1(:,i),U1(:,i+1))
        END DO
        U2(:,N)=U1(:,N)-Space_Term(U1(:,N-1),U1(:,N),U1(:,N))
      ELSE IF (id_temps==2) THEN
        !------------------K1------------------------------
        k1(:,1)=-Space_Term(U1(:,1),U1(:,1),U1(:,2))
        DO i=2,N-1
          k1(:,i)=-Space_Term(U1(:,i-1),U1(:,i),U1(:,i+1))
        END DO
        k1(:,N)=-Space_Term(U1(:,N-1),U1(:,N),U1(:,N))
        !------------------K2------------------------------
        U=U1+(2.0/3.0)*k1
        k2(:,1)=-Space_Term(U(:,1),U(:,1),U(:,2))
        DO i=2,N-1
          k2(:,i)=-Space_Term(U(:,i-1),U(:,i),U(:,i+1))
        END DO
        k2(:,N)=-Space_Term(U(:,N-1),U(:,N),U(:,N))
        U2=U1+k2
          ELSE IF (id_temps==3) THEN
        k1(:,1)=-Space_Term(U1(:,1),U1(:,1),U1(:,2))
        DO i=2,N-1
          k1(:,i)=-Space_Term(U1(:,i-1),U1(:,i),U1(:,i+1))
        END DO
        k1(:,N)=-Space_Term(U1(:,N-1),U1(:,N),U1(:,N))
        U=U1+(2.0/3.0)*k1
        k2(:,1)=-Space_Term(U(:,1),U(:,1),U(:,2))
        DO i=2,N-1
          k2(:,i)=-Space_Term(U(:,i-1),U(:,i),U(:,i+1))
        END DO
        k2(:,N)=-Space_Term(U(:,N-1),U(:,N),U(:,N))

        U2=U1+(3.0/4.0)*k2+(1.0/4.0)*k1
      END IF
    ELSE IF (id_ordre==2) THEN
      IF (id_temps==1) THEN
        U2(:,1)=U1(:,1)-Space_Term_2(U1(:,1),U1(:,1),U1(:,1),U1(:,2),U1(:,3))
        U2(:,2)=U1(:,2)-Space_Term_2(U1(:,1),U1(:,1),U1(:,2),U1(:,3),U1(:,4))
        DO i=3,N-2
          U2(:,i)=U1(:,i)-Space_Term_2(U1(:,i-2),U1(:,i-1),U1(:,i),U1(:,i+1),U1(:,i+2))
        END DO
        U2(:,N-1)=U1(:,N-1)-Space_Term_2(U1(:,N-3),U1(:,N-2),U1(:,N-1),U1(:,N),U1(:,N))
        U2(:,N)=U1(:,N)-Space_Term_2(U1(:,N-2),U1(:,N-1),U1(:,N),U1(:,N),U1(:,N))
      ELSE IF (id_temps==2) THEN
        !------------------K1------------------------------
        k1(:,1)=-Space_Term_2(U1(:,1),U1(:,1),U1(:,1),U1(:,2),U1(:,3))
        k1(:,2)=-Space_Term_2(U1(:,1),U1(:,1),U1(:,2),U1(:,3),U1(:,4))
        DO i=3,N-2
          k1(:,i)=-Space_Term_2(U1(:,i-2),U1(:,i-1),U1(:,i),U1(:,i+1),U1(:,i+2))
        END DO
        k1(:,N-1)=-Space_Term_2(U1(:,N-3),U1(:,N-2),U1(:,N-1),U1(:,N),U1(:,N))
        k1(:,N)=-Space_Term_2(U1(:,N-2),U1(:,N-1),U1(:,N),U1(:,N),U1(:,N))
        !------------------K2------------------------------
        U=U1+(2.0/3.0)*k1
        k2(:,1)=-Space_Term_2(U(:,1),U(:,1),U(:,1),U(:,2),U(:,3))
        k2(:,2)=-Space_Term_2(U(:,1),U(:,1),U(:,2),U(:,3),U(:,4))
        DO i=3,N-2
          k2(:,i)=-Space_Term_2(U(:,i-2),U(:,i-1),U(:,i),U(:,i+1),U(:,i+2))
        END DO
        k2(:,N-1)=-Space_Term_2(U(:,N-3),U(:,N-2),U(:,N-1),U(:,N),U(:,N))
        k2(:,N)=-Space_Term_2(U(:,N-2),U(:,N-1),U(:,N),U(:,N),U(:,N))
        U2=U1+k2
        
          ELSE IF (id_temps==3) THEN

          !------------------K1------------------------------
          k1(:,1)=-Space_Term_2(U1(:,1),U1(:,1),U1(:,1),U1(:,2),U1(:,3))
          k1(:,2)=-Space_Term_2(U1(:,1),U1(:,1),U1(:,2),U1(:,3),U1(:,4))
          DO i=3,N-2
            k1(:,i)=-Space_Term_2(U1(:,i-2),U1(:,i-1),U1(:,i),U1(:,i+1),U1(:,i+2))
          END DO
          k1(:,N-1)=-Space_Term_2(U1(:,N-3),U1(:,N-2),U1(:,N-1),U1(:,N),U1(:,N))
          k1(:,N)=-Space_Term_2(U1(:,N-2),U1(:,N-1),U1(:,N),U1(:,N),U1(:,N))
          !------------------K2------------------------------
          U=U1+(2.0/3.0)*k1
          k2(:,1)=-Space_Term_2(U(:,1),U(:,1),U(:,1),U(:,2),U(:,3))
          k2(:,2)=-Space_Term_2(U(:,1),U(:,1),U(:,2),U(:,3),U(:,4))
          DO i=3,N-2
            k2(:,i)=-Space_Term_2(U(:,i-2),U(:,i-1),U(:,i),U(:,i+1),U(:,i+2))
          END DO
          k2(:,N-1)=-Space_Term_2(U(:,N-3),U(:,N-2),U(:,N-1),U(:,N),U(:,N))
          k2(:,N)=-Space_Term_2(U(:,N-2),U(:,N-1),U(:,N),U(:,N),U(:,N))
          U2=U1+k2*(3.0/4.0)+(1.0/4.0)*k1



      END IF
    END IF
  END SUBROUTINE Time_Step

  FUNCTION Space_Term(U0,U1,U2)
    REAL*8,DIMENSION(3),INTENT(IN)::U0,U1,U2
    REAL*8,DIMENSION(3)::Space_Term
    Space_Term=Flux_HLL(U1,U2)-Flux_HLL(U0,U1)
    Space_Term=Space_Term/dx*dt
  END FUNCTION

  FUNCTION Space_Term_2(U_,U0,U1,U2,U3)
    REAL*8,DIMENSION(3),INTENT(IN)::U_,U0,U1,U2,U3
    REAL*8,DIMENSION(3)::Space_Term_2
    Space_Term_2=Flux_MUSCL(U0,U1,U2,U3)-Flux_MUSCL(U_,U0,U1,U2)
    Space_Term_2=Space_Term_2/dx*dt
  END FUNCTION

  SUBROUTINE DT_CFL(U)
    IMPLICIT NONE
    REAL*8,DIMENSION(3,N),INTENT(IN)::U
    REAL*8,DIMENSION(N)::VALEURS
    REAL*8::valeur,dt_intermediaire
    INTEGER::i
    DO i=1,N
      VALEURS(i)=abs(U(2,i)/U(1,i))+sqrt(gamma*Pression(U(:,i))/U(1,i))
    END DO
    valeur=maxval(VALEURS)
    dt_intermediaire=cfl*dx/valeur
    dt=min(dt_min,dt_intermediaire)
    ! PRINT*,'TEST EVALUATION DU PAS DE TEMPS'
    ! PRINT*,'dt_intermediaire',dt_intermediaire,'dt',dt
  END SUBROUTINE

  SUBROUTINE Print_Densite(U,Nb)
    IMPLICIT NONE
    REAL*8,DIMENSION(3,N),INTENT(IN)::U
    INTEGER,INTENT(IN)::Nb
    INTEGER::i
    CHARACTER(len=36)::F_NAME
    WRITE(F_NAME,"(A,I0,'.dat')") 'Result/Densite',Nb
    OPEN(UNIT=2,file=F_NAME,action='write')
    DO i=1,N
      WRITE(2,*),(i-0.5)*dx, U(1,i)
    END DO
    CLOSE(2)
  END SUBROUTINE Print_Densite

  SUBROUTINE Print_Vitesse(U,Nb)
    IMPLICIT NONE
    REAL*8,DIMENSION(3,N),INTENT(IN)::U
    INTEGER,INTENT(IN)::Nb
    INTEGER::i
    CHARACTER(len=36)::F_NAME

    WRITE(F_NAME,"(A,I0,'.dat')") 'Result/Vitesse',Nb
    OPEN(UNIT=2,file=F_NAME,action='write')
    DO i=1,N
      WRITE(2,*),(i-0.5)*dx, U(2,i)/U(1,i)
    END DO
    CLOSE(2)
  END SUBROUTINE Print_Vitesse

  SUBROUTINE Print_Pression(U,Nb)
    IMPLICIT NONE
    REAL*8,DIMENSION(3,N),INTENT(IN)::U
    INTEGER,INTENT(IN)::Nb
    INTEGER::i
    CHARACTER(len=36)::F_NAME

    WRITE(F_NAME,"(A,I0,'.dat')") 'Result/Pression',Nb
    OPEN(UNIT=2,file=F_NAME,action='write')
    DO i=1,N
      WRITE(2,*),(i-0.5)*dx, Pression(U(:,i))
    END DO
    CLOSE(2)
  END SUBROUTINE Print_Pression

END MODULE module
