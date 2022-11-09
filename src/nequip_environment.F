!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with Neural Network potentials
!> \author Gabriele Tocci (gabriele.tocci@uzh.ch)
!> \date   2022-09-11
! **************************************************************************************************
MODULE nequip_environment

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_methods,                    ONLY: read_cell,&
                                              write_cell
   USE cell_types,                      ONLY: cell_release,&
                                              cell_type,&
                                              get_cell
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_parser_methods,               ONLY: parser_read_line,&
                                              parser_search_string
   USE cp_parser_types,                 ONLY: cp_parser_type,&
                                              parser_create,&
                                              parser_release,&
                                              parser_reset
   USE cp_subsys_methods,               ONLY: cp_subsys_create
   USE cp_subsys_types,                 ONLY: cp_subsys_release,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                              write_molecule_kind_set
   USE molecule_types,                  ONLY: molecule_type
   USE nnp_model,                       ONLY: nnp_write_arc
   USE particle_methods,                ONLY: write_fist_particle_coordinates,&
                                              write_particle_distances,&
                                              write_structure_data
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nequip_environment'

   PUBLIC :: nequip_init

CONTAINS

! **************************************************************************************************
!> \brief Read and initialize all the information for neural network potentials
!> \param nnp_env ...
!> \param root_section ...
!> \param para_env ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nequip_init(nequip_env, root_section, para_env, force_env_section, subsys_section, &
                       use_motion_section)
      TYPE(nequip_type), INTENT(INOUT), POINTER             :: nequip_env
      TYPE(section_vals_type), INTENT(IN), POINTER       :: root_section
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(section_vals_type), INTENT(INOUT), POINTER    :: force_env_section, subsys_section
      LOGICAL, INTENT(IN)                                :: use_motion_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nequip_init'

      INTEGER                                            :: handle
      LOGICAL                                            :: explicit, use_ref_cell
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: cell_section, nnp_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(nnp_env))

      NULLIFY (cell_section, nnp_section, cell, cell_ref, subsys)

      IF (.NOT. ASSOCIATED(subsys_section)) THEN
         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      END IF
      cell_section => section_vals_get_subs_vals(subsys_section, "CELL")
      nnp_section => section_vals_get_subs_vals(force_env_section, "NEQUIP")
      CALL section_vals_get(nnp_section, explicit=explicit)
      IF (.NOT. explicit) THEN
         CPWARN("NEQUIP section not explicitly stated. Using default file names.")
      END IF

      CALL nnp_env_set(nnp_env=nnp_env, nnp_input=nnp_section, &
                       force_env_input=force_env_section)

      CALL read_cell(cell=cell, cell_ref=cell_ref, use_ref_cell=use_ref_cell, cell_section=cell_section, &
                     para_env=para_env)
      CALL get_cell(cell=cell, abc=abc)
      CALL write_cell(cell=cell, subsys_section=subsys_section)

      CALL cp_subsys_create(subsys, para_env, root_section, &
                            force_env_section=force_env_section, subsys_section=subsys_section, &
                            use_motion_section=use_motion_section)

      CALL nnp_init_subsys(nnp_env=nnp_env, subsys=subsys, cell=cell, &
                           cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                           subsys_section=subsys_section)

      CALL cell_release(cell)
      CALL cell_release(cell_ref)
      CALL cp_subsys_release(subsys)

      CALL timestop(handle)

   END SUBROUTINE nequip_init

! **************************************************************************************************
!> \brief Read and initialize all the information for neural network potentials
!> \param nnp_env ...
!> \param subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param subsys_section ...
!> \date   2022-9-11
!> \author Gabriele Tocci (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nequip_init_subsys(nnp_env, subsys, cell, cell_ref, use_ref_cell, subsys_section)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env
      TYPE(cp_subsys_type), INTENT(IN), POINTER          :: subsys
      TYPE(cell_type), INTENT(INOUT), POINTER            :: cell, cell_ref
      LOGICAL, INTENT(IN)                                :: use_ref_cell
      TYPE(section_vals_type), INTENT(IN), POINTER       :: subsys_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nnp_init_subsys'

      INTEGER                                            :: handle, natom
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, molecule_kind_set, particle_set, molecule_set, &
               local_molecules, local_particles)

      particle_set => subsys%particles%els
      atomic_kind_set => subsys%atomic_kinds%els
      molecule_kind_set => subsys%molecule_kinds%els
      molecule_set => subsys%molecules%els

      !Print the molecule kind set
      CALL write_molecule_kind_set(molecule_kind_set, subsys_section)

      !Print the atomic coordinates
      CALL write_fist_particle_coordinates(particle_set, subsys_section)
      CALL write_particle_distances(particle_set, cell=cell, &
                                    subsys_section=subsys_section)
      CALL write_structure_data(particle_set, cell=cell, &
                                input_section=subsys_section)

      !Distribute molecules and atoms using the new data structures
      CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   local_particles=local_particles, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   local_molecules=local_molecules, &
                                   force_env_section=nnp_env%force_env_input)

      natom = SIZE(particle_set)

      ALLOCATE (nnp_env%nnp_forces(3, natom))

      nnp_env%nnp_forces(:, :) = 0.0_dp

      nnp_env%nnp_potential_energy = 0.0_dp

      ! Set up arrays for calculation:
      nnp_env%num_atoms = natom
      ALLOCATE (nnp_env%ele_ind(natom))
      ALLOCATE (nnp_env%nuc_atoms(natom))
      ALLOCATE (nnp_env%coord(3, natom))
      ALLOCATE (nnp_env%atoms(natom))
      ALLOCATE (nnp_env%sort(natom))
      ALLOCATE (nnp_env%sort_inv(natom))

      CALL cp_subsys_set(subsys, cell=cell)

      CALL nnp_env_set(nnp_env=nnp_env, subsys=subsys, &
                       cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                       local_molecules=local_molecules, &
                       local_particles=local_particles)

      CALL distribution_1d_release(local_particles)
      CALL distribution_1d_release(local_molecules)

      CALL nnp_init_model(nnp_env=nnp_env)

      CALL timestop(handle)

   END SUBROUTINE nequip_init_subsys

! **************************************************************************************************
!> \brief Initialize the Neural Network Potential
!> \param nnp_env ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nequip_init_model(nnp_env)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'nequip_init_model'
      INTEGER, PARAMETER                                 :: def_str_len = 256, &
                                                            default_path_length = 256

      CHARACTER(len=1), ALLOCATABLE, DIMENSION(:)        :: cactfnct
      CHARACTER(len=2)                                   :: ele
      CHARACTER(len=def_str_len)                         :: dummy, line
      CHARACTER(len=default_path_length)                 :: base_name, file_name
      INTEGER                                            :: handle, i, i_com, io, iweight, j, k, l, &
                                                            n_weight, nele, nuc_ele, symfnct_type, &
                                                            unit_nr
      LOGICAL                                            :: at_end, atom_e_found, explicit, first, &
                                                            found
      REAL(KIND=dp)                                      :: energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: weights
      REAL(KIND=dp), DIMENSION(7)                        :: test_array
      REAL(KIND=dp), DIMENSION(:), POINTER               :: work
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_parser_type)                               :: parser
      TYPE(section_vals_type), POINTER                   :: bias_section, model_section

      CALL timeset(routineN, handle)

      NULLIFY (logger)

      logger => cp_get_default_logger()

      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger)
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "NEQUIP| Neural Network Potential Force Environment"
      END IF

      model_section => section_vals_get_subs_vals(nnp_env%nnp_input, "MODEL")
      CALL section_vals_get(model_section, n_repetition=nnp_env%n_committee)
      ALLOCATE (nnp_env%atomic_energy(nnp_env%num_atoms, nnp_env%n_committee))
      ALLOCATE (nnp_env%myforce(3, nnp_env%num_atoms, nnp_env%n_committee))

      CALL parser_create(parser, file_name, para_env=logger%para_env)

      nele = nnp_env%n_ele

      ALLOCATE (nnp_env%atom_energies(nele))
      nnp_env%atom_energies = 0.0_dp


      IF (logger%para_env%ionode) THEN
         WRITE (unit_nr, *) "NEQUIP| NEQUIP force environment initialized"
      END IF

      CALL timestop(handle)

   END SUBROUTINE nequip_init_model

END MODULE nequip_environment
