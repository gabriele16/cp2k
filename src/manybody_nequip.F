!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE manybody_nequip

   USE cp_log_handling, ONLY: cp_logger_get_default_io_unit
   USE atomic_kind_types, ONLY: atomic_kind_type

   USE cell_types, ONLY: cell_type
   USE cp_para_types, ONLY: cp_para_env_type
   USE fist_nonbond_env_types, ONLY: fist_nonbond_env_get, &
                                     fist_nonbond_env_set, &
                                     fist_nonbond_env_type, &
                                     nequip_data_type
   USE kinds, ONLY: dp
   USE pair_potential_types, ONLY: pair_potential_pp_type, &
                                   pair_potential_single_type, &
                                   nequip_pot_type, &
                                   nequip_type
   USE particle_types, ONLY: particle_type
   USE physcon, ONLY: angstrom, &
                      evolt
#ifdef __NEQUIP
   USE nequip_wrapper, ONLY: nequip_nnp, &
                             create_nequip, &
                             delete_nequip_c, &
                             compute_nequip
#endif

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PUBLIC nequip_energy_store_force_virial, nequip_add_force

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_nequip'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param cell ...
!> \param atomic_kind_set ...
!> \param potparm ...
!> \param fist_nonbond_env ...
!> \param pot_nequip ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE nequip_energy_store_force(particle_set, cell, atomic_kind_set, potparm, fist_nonbond_env, &
                                               pot_nequip, para_env)
      TYPE(particle_type), POINTER             :: particle_set(:)
      TYPE(cell_type), POINTER                 :: cell
      TYPE(atomic_kind_type), POINTER          :: atomic_kind_set(:)
      TYPE(pair_potential_pp_type), POINTER    :: potparm
      TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
      REAL(kind=dp), TARGET                    :: pot_nequip
      TYPE(cp_para_env_type), OPTIONAL, &
         POINTER                                :: para_env

      CHARACTER(LEN=*), PARAMETER :: &
         routineN = 'nequip_energy_store_force', &
         routineP = moduleN//':'//routineN

#ifdef __NEQUIP
      INTEGER                                  :: i, iat, iat_use, ikind, &
                                                  jkind, n_atoms, n_atoms_use, &
                                                  output_unit
      LOGICAL                                  :: do_parallel
      LOGICAL, ALLOCATABLE                     :: use_atom(:)
      REAL(kind=dp)                            :: lattice(3, 3)
      TYPE(pair_potential_single_type), &
         POINTER                                :: pot
      TYPE(nequip_data_type), POINTER          :: nequip_data
      TYPE(nequip_pot_type), POINTER           :: nequip
      REAL(kind=dp), ALLOCATABLE, TARGET       :: dpmd_force(:), &
                                                  dpmd_atom_ener(:), &
                                                  dpmd_box(:), dpmd_coord(:)
      INTEGER, ALLOCATABLE, TARGET             :: dpmd_atype(:), use_atom_type(:)
      INTEGER, TARGET                          :: vecsize
      REAL(kind=dp), POINTER                   :: dforce(:), datom_ener(:), &
                                                  dbox(:), dcoord(:), dener
      INTEGER, POINTER                         :: datype(:)

#endif
#ifndef __NEQUIP
      MARK_USED(particle_set)
      MARK_USED(cell)
      MARK_USED(atomic_kind_set)
      MARK_USED(potparm)
      MARK_USED(fist_nonbond_env)
      MARK_USED(pot_nequip)
      MARK_USED(para_env)
      CALL cp_abort(__LOCATION__, "In order to use NequIP you need to download "// &
                    "and install PyTorch (check CP2K manual for details)")
#else
      n_atoms = SIZE(particle_set)
      ALLOCATE (use_atom(n_atoms))
      ALLOCATE (use_atom_type(n_atoms))
      use_atom = .FALSE.
      use_atom_type = 100 ! just initial the array

      NULLIFY (nequip)

      DO ikind = 1, SIZE(atomic_kind_set)
      DO jkind = 1, SIZE(atomic_kind_set)
         pot => potparm%pot(ikind, jkind)%pot
         DO i = 1, SIZE(pot%type)
            IF (pot%type(i) /= nequip_type) CYCLE
            IF (.NOT. ASSOCIATED(nequip)) nequip => pot%set(i)%nequip
            DO iat = 1, n_atoms
               IF (particle_set(iat)%atomic_kind%kind_number == ikind .OR. &
                   particle_set(iat)%atomic_kind%kind_number == jkind) THEN
                  use_atom(iat) = .TRUE.
                  use_atom_type(iat) = pot%set(i)%nequip%atom_nequip_type
               END IF
            END DO ! iat
         END DO ! i
      END DO ! jkind
      END DO ! ikind
!     print*, use_atom_type
      n_atoms_use = COUNT(use_atom)
      !      ALLOCATE (pos(3, n_atoms_use), force(3, n_atoms_use), elem_symbol(n_atoms_use))
      vecsize = n_atoms_use
      ALLOCATE (dpmd_box(9))
      ALLOCATE (dpmd_atype(vecsize))
      ALLOCATE (dpmd_coord(vecsize*3))
      ALLOCATE (dpmd_force(vecsize*3))
      ALLOCATE (dpmd_atom_ener(vecsize))

! initial the properties
      pot_nequip = 0.0
      dpmd_force = 0.0
      dpmd_atom_ener = 0.0
      iat_use = 0
      DO iat = 1, n_atoms
         IF (.NOT. use_atom(iat)) CYCLE
         iat_use = iat_use + 1
         dpmd_coord((iat_use - 1)*3 + 1:(iat_use - 1)*3 + 3) = particle_set(iat)%r*angstrom
         dpmd_atype(iat_use) = use_atom_type(iat)
      END DO
      IF (iat_use > 0) 
      output_unit = cp_logger_get_default_io_unit()
      lattice = cell%hmat*angstrom
! change matrix to one d array
      DO i = 1, 3
         dpmd_box((i - 1)*3 + 1:(i - 1)*3 + 3) = lattice(:, i)
      END DO
      do_parallel = .FALSE.
      IF (PRESENT(para_env)) THEN
         do_parallel = para_env%num_pe > 1
      END IF
! change to pointer type
      dener => pot_nequip
      dforce => dpmd_force
      datom_ener => dpmd_atom_ener
      dcoord => dpmd_coord
      dbox => dpmd_box
      datype => dpmd_atype
! get nequip_data to save shared state as well as force, virial info
      CALL fist_nonbond_env_get(fist_nonbond_env, nequip_data=nequip_data)
      IF (.NOT. ASSOCIATED(nequip_data)) THEN
         ALLOCATE (nequip_data)
         CALL fist_nonbond_env_set(fist_nonbond_env, nequip_data=nequip_data)
         NULLIFY (nequip_data%use_indices, nequip_data%force)
      END IF
      IF (ASSOCIATED(nequip_data%force)) THEN
         IF (SIZE(nequip_data%force, 2) /= n_atoms_use) THEN
            DEALLOCATE (nequip_data%force, nequip_data%use_indices)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(nequip_data%force)) THEN
         ALLOCATE (nequip_data%force(3, n_atoms_use))
         ALLOCATE (nequip_data%use_indices(n_atoms_use))
      END IF
! initial the two different potential if needed
      IF (nequip_data%do_count == 0) THEN
         IF (nequip_data%do_create(1)) THEN
            nequip_data%nequip_tmp_pot1 = create_nequip(TRIM(nequip%nequip_file_name))
            nequip_data%do_create(1) = .FALSE.
         END IF
         CALL compute_nequip( &
            pot=nequip_data%nequip_tmp_pot1%ptr, &
            vecsize=vecsize, &
            dener=dener, &
            dforce=dforce, &
            datom_ener=datom_ener, &
            dcoord=dcoord, &
            datype=datype, &
            dbox=dbox)
        nequip_data%do_count = 1
      ELSE IF (nequip_data%do_count == 1) THEN
         IF (nequip_data%do_create(1)) THEN
            nequip_data%nequip_tmp_pot2 = create_nequip(TRIM(nequip%nequip_file_name))
            nequip_data%do_create(2) = .FALSE.
         END IF
         CALL compute_nequip( &
            pot=nequip_data%nequip_tmp_pot2%ptr, &
            vecsize=vecsize, &
            dener=dener, &
            dforce=dforce, &
            datom_ener=datom_ener, &
            dcoord=dcoord, &
            datype=datype, &
            dbox=dbox)
         nequip_data%do_count = 0
      END IF

!      ENDIF
      ! convert units
      pot_nequip = pot_nequip/evolt
      dpmd_force = dpmd_force/(evolt/angstrom)
      ! save force, virial info
      iat_use = 0
      DO iat = 1, n_atoms
         IF (use_atom(iat)) THEN
            iat_use = iat_use + 1
            nequip_data%use_indices(iat_use) = iat
            nequip_data%force(1:3, iat_use) = dpmd_force((iat_use - 1)*3 + 1:(iat_use - 1)*3 + 3)
         END IF
      END DO
      DEALLOCATE (use_atom, use_atom_type, dpmd_coord, dpmd_force, &
                   dpmd_atom_ener,  &
                  dpmd_box, dpmd_atype)
#endif
   END SUBROUTINE nequip_energy_store_force

! **************************************************************************************************
!> \brief ...
!> \param fist_nonbond_env ...
!> \param force ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE nequip_add_force(fist_nonbond_env, force)
      TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
      REAL(KIND=dp)                            :: force(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nequip_add_force', &
                                     routineP = moduleN//':'//routineN

#ifdef __NEQUIP
      INTEGER                                  :: iat, iat_use
      TYPE(nequip_data_type), POINTER            :: nequip_data
#endif

#ifndef __NEQUIP
      MARK_USED(fist_nonbond_env)
      MARK_USED(force)
      RETURN
#else
      CALL fist_nonbond_env_get(fist_nonbond_env, nequip_data=nequip_data)
      IF (.NOT. ASSOCIATED(nequip_data)) RETURN

      DO iat_use = 1, SIZE(nequip_data%use_indices)
         iat = nequip_data%use_indices(iat_use)
         CPASSERT(iat >= 1 .AND. iat <= SIZE(force, 2))
         force(1:3, iat) = force(1:3, iat) + nequip_data%force(1:3, iat_use)
      END DO
#endif
   END SUBROUTINE nequip_add_force

END MODULE nequip_force
