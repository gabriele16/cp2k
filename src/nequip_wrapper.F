!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the NequIP.
!> \par History
!>      07.2019 created [Yongbin Zhuang]
!>      06.2021 refactored [Yunpei Liu]
!>      11.2022 adapted to NequIP [Gabriele Tocci]
!> \author Yongbin Zhuang
! **************************************************************************************************

MODULE nequip_wrapper

   USE ISO_C_BINDING,                   ONLY: C_CHAR,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_NULL_CHAR,&
                                              C_PTR
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: nequip_nnp, create_nequip, delete_nequip_c, compute_nequip

#if(__NEQUIP)
   INTERFACE
      FUNCTION create_nequip_c(model) BIND(C, name="create_nequip")
         USE ISO_C_BINDING, ONLY: C_CHAR, C_PTR 
      CHARACTER(KIND=C_CHAR)                      :: model(*)
      TYPE(C_PTR)                                 :: create_nequip_c

      END FUNCTION
      SUBROUTINE delete_nequip_c(nequip_nnp) BIND(C, name="delete_nequip")
         USE ISO_C_BINDING, ONLY: C_PTR
      TYPE(C_PTR), INTENT(IN), VALUE                     :: nequip_nnp

      END SUBROUTINE
      SUBROUTINE compute_nequip_c(nequip_nnp, vecsize, &
                               dener, dforce, datom_ener, &
                               dcoord, datype, dbox) BIND(C, name="compute_nequip")
         USE ISO_C_BINDING, ONLY : C_PTR
      TYPE(C_PTR), VALUE                     :: nequip_nnp, vecsize, dener, dforce, &
                                                            datom_ener, dcoord, &
                                                            datype, dbox

      END SUBROUTINE
   END INTERFACE
#endif

   TYPE nequip_nnp
      TYPE(C_PTR) :: ptr
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param model ...
!> \return ...
! **************************************************************************************************
   FUNCTION create_nequip(model)
      CHARACTER(len=*), INTENT(IN), TARGET               :: model
      TYPE(nequip_nnp)                                   :: create_nequip
#if(__NEQUIP)
      create_nequip%ptr = create_nequip_c(TRIM(model)//C_NULL_CHAR)
#endif
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param pot ...
!> \param vecsize ...
!> \param dener ...
!> \param dforce ...
!> \param datom_ener ...
!> \param dcoord ...
!> \param datype ...
!> \param dbox ...
! **************************************************************************************************
   SUBROUTINE compute_nequip(pot, vecsize, dener, dforce, datom_ener, dcoord, datype, dbox)
      IMPLICIT NONE
      TYPE(C_PTR) :: pot
      INTEGER(C_INT), TARGET  :: vecsize
      REAL(C_DOUBLE), POINTER :: dener
      REAL(C_DOUBLE), POINTER :: dforce(:)
      REAL(C_DOUBLE), POINTER :: datom_ener(:)
      REAL(C_DOUBLE), POINTER :: dcoord(:)
      INTEGER(C_INT), POINTER :: datype(:)
      REAL(C_DOUBLE), POINTER :: dbox(:)
#if(__NEQUIP)
      CALL compute_nequip_c(pot, C_LOC(vecsize), C_LOC(dener), C_LOC(dforce(1)), &
                         C_LOC(datom_ener(1)), C_LOC(dcoord(1)), C_LOC(datype(1)), &
                         C_LOC(dbox(1)))
                         
#endif
!             print*, "this is energy"
!             print*, dener
!             print*, "this is force"
!             print*, dforce
!             print*, "this is atom_energy"
!             print*, datom_ener
!             print*, 'after call'
   END SUBROUTINE
END MODULE nequip_wrapper
