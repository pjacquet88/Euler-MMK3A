MODULE limiteurs
  IMPLICIT NONE

  INTEGER,PUBLIC :: id_limiteur

  PUBLIC :: limiteur
  PRIVATE :: Superbee,Minmod,VanLeer,Ospre

CONTAINS

  FUNCTION Limiteur(r)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::r
    REAL*8::Limiteur


    IF (id_limiteur==0) THEN
      limiteur=1.0
    ELSE IF (id_limiteur==1) THEN
      limiteur=Superbee(r)
    ELSE IF (id_limiteur==2) THEN
      limiteur=VanLeer(r)
    ELSE IF (id_limiteur==3) THEN
      limiteur=Minmod(r)
    ELSE IF (id_limiteur==4) THEN
      limiteur=Ospre(r)
    END IF
  END FUNCTION


  !---------------------PRIVATE---------------------------------------
  FUNCTION Superbee(r)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::r
    REAL*8::Superbee
    Superbee=max(0.0d0,min(2*r,1.0d0),min(r,2.0d0))
  !  print*,'Superbee r : ',r ,Superbee
  END FUNCTION

  FUNCTION Minmod(r)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::r
    REAL*8::Minmod
    Minmod=max(0.0d0,min(1.0d0,r))
  !  print*,'Minmod r : ',r ,Minmod
  END FUNCTION


  FUNCTION VanLeer(r)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::r
    REAL*8::VanLeer
    VanLeer=(r+abs(r))/(1+abs(r))
    !print*,'VanLeer r : ',r ,VanLeer
  END FUNCTION

  FUNCTION Ospre(r)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::r
    REAL*8::Ospre
    Ospre=1.5*(r**2+r)/(r**2+r+1)
    !print*,'Ospre r : ',r ,Ospre
  END FUNCTION

END MODULE limiteurs
