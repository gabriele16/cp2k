!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Creates the NEQUIP section of the input
!> \author Gabriele Tocci (gabriele.tocci@uzh.chem)
!> \date   2022-06-10
! **************************************************************************************************
MODULE input_cp2k_nequip

   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              medium_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              real_t
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_nequip'

   PUBLIC :: create_nequip_section

CONTAINS

! **************************************************************************************************
!> \brief Create the input section for NEQUIP
!> \param section the section to create
!> \date   2022-06-10
!> \author Gabriele Tocci (gabriele.tocci@uzh.ch)
! **************************************************************************************************
   SUBROUTINE create_nequip_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_nequip_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection, subsubsection

!   ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="NEQUIP", &
                          description="This section contains all information to run a "// &
                          "Nequip calculation.", &
                          n_keywords=0, n_subsections=2, repeats=.FALSE.)

      NULLIFY (subsection, subsubsection, keyword)

      CALL section_create(subsection, __LOCATION__, name="MODEL", &
                          description="Section for a single nequip model. "// &
                          "If this section is repeated, a committee model (C-Nequip)"// &
                          "is used where the nequip members have different weights.", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)
      CALL keyword_create(keyword, __LOCATION__, name="NEQUIP_MODEL_FILE_NAME", &
                          description="File containing the weights and model definition "// &
                          "artificial neural networks of nequip. "// &
                          "The extension is a PyTorch file .XXX.pth", &
                          repeats=.FALSE., default_lc_val="nequip_model")
      CALL section_add_keyword(subsection, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_nequip_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_nequip_section

! **************************************************************************************************
!> \brief Creates the print section for the nequip subsection
!> \param section the section to create
!> \date   2022-10-06
!> \author Gabriele Tocci (gabriele.tocci@uzh.ch)
! **************************************************************************************************
   SUBROUTINE create_nequip_print_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_nequip_print_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

!   ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options in NEQUIP code.", &
                          n_keywords=0, n_subsections=5, repeats=.FALSE.)

      NULLIFY (print_key, keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "ENERGIES", &
                                       description="Controls the printing of the NEQUIP energies.", &
                                       print_level=medium_print_level, common_iter_levels=1)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_nequip_print_section

END MODULE input_cp2k_nequip
