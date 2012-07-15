!> \file
!> \author Chris Bradley
!> \brief This module contains all interface conditions routines.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>   
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>This module contains all interface conditions routines.
MODULE INTERFACE_CONDITIONS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE COMP_ENVIRONMENT
  USE DATA_POINT_ROUTINES
  USE DATA_PROJECTION_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_EQUATIONS_ROUTINES
  USE INTERFACE_MAPPING_ROUTINES
  USE INTERFACE_MATRICES_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  !Module types

  !Module variables

  !Interfaces

  PUBLIC INTERFACE_CONDITION_CREATE_FINISH,INTERFACE_CONDITION_CREATE_START

  PUBLIC INTERFACE_CONDITION_DEPENDENT_VARIABLE_ADD

  PUBLIC INTERFACE_CONDITION_DESTROY

  PUBLIC INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH,INTERFACE_CONDITION_EQUATIONS_CREATE_START

  PUBLIC INTERFACE_CONDITION_EQUATIONS_DESTROY

  PUBLIC INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH,INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START

  PUBLIC INTERFACE_CONDITION_METHOD_GET,INTERFACE_CONDITION_METHOD_SET

  PUBLIC INTERFACE_CONDITION_OPERATOR_GET,INTERFACE_CONDITION_OPERATOR_SET
  
  PUBLIC INTERFACE_CONDITION_TRANSLATION_INCREMENT_APPLY

  PUBLIC INTERFACE_CONDITION_USER_NUMBER_FIND

  PUBLIC INTERFACE_CONDITIONS_FINALISE,INTERFACE_CONDITIONS_INITIALISE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assembles the equations for an interface condition.
  SUBROUTINE INTERFACE_CONDITION_ASSEMBLE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            IF(INTERFACE_CONDITION%OPERATOR==INTERFACE_CONDITION_FRICTIONLESS_CONTACT_OPERATOR) THEN
              IF (INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD%VARIABLES(1)%COMPONENTS(1)%INTERPOLATION_TYPE== &
                  & FIELD_DATA_POINT_BASED_INTERPOLATION) THEN
                CALL INTERFACE_CONDITION_DATA_REPROJECTION(INTERFACE_CONDITION,ERR,ERROR,*999)
                CALL INTERFACE_CONDITION_DATA_POINTS_NORMAL_CALCULATE(INTERFACE_CONDITION,ERR,ERROR,*999)
              ENDIF
            ENDIF
            CALL INTERFACE_CONDITION_ASSEMBLE_FEM(INTERFACE_CONDITION,ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Interface equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_ASSEMBLE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_ASSEMBLE

  !
  !================================================================================================================================
  !
  
  !>Assembles the interface matricesand rhs for using the finite element method.
  SUBROUTINE INTERFACE_CONDITION_ASSEMBLE_FEM(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    
!#ifdef TAUPROF
!    CHARACTER(28) :: CVAR
!    INTEGER :: PHASE(2) = (/ 0, 0 /)
!    SAVE PHASE
!#endif

    CALL ENTERS("INTERFACE_CONDITION_ASSEMBLE_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
        LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
        IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
            IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              !Initialise the matrices and rhs vector
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_VALUES_INITIALISE()")
#endif
              CALL INTERFACE_MATRICES_VALUES_INITIALISE(INTERFACE_MATRICES,0.0_DP,ERR,ERROR,*999)
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_VALUES_INITIALISE()")
#endif
              !Assemble the elements
              !Allocate the element matrices 
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_INITIALISE()")
#endif
              CALL INTERFACE_MATRICES_ELEMENT_INITIALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
              ELEMENTS_MAPPING=>LAGRANGE_FIELD%DECOMPOSITION%DOMAIN(LAGRANGE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
                & MAPPINGS%ELEMENTS
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_INITIALISE()")
#endif
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for interface equations setup and initialisation = ", &
                  & USER_ELAPSED,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for interface equations setup and initialisation = ", &
                  & SYSTEM_ELAPSED,ERR,ERROR,*999)
                ELEMENT_USER_ELAPSED=0.0_SP
                ELEMENT_SYSTEM_ELAPSED=0.0_SP
              ENDIF
              NUMBER_OF_TIMES=0
              !Loop over the internal elements

#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("Internal Elements Loop")
#endif
              DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
!#ifdef TAUPROF
!              WRITE (CVAR,'(a23,i3)') 'Internal Elements Loop ',element_idx
!              CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
!              CALL TAU_PHASE_START(PHASE)
!#endif
                ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
                NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
                CALL INTERFACE_MATRICES_ELEMENT_CALCULATE(INTERFACE_MATRICES,ne,ERR,ERROR,*999)
                CALL INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ne,ERR,ERROR,*999)
                CALL INTERFACE_MATRICES_ELEMENT_ADD(INTERFACE_MATRICES,ERR,ERROR,*999)
!#ifdef TAUPROF
!              CALL TAU_PHASE_STOP(PHASE)
!#endif
              ENDDO !element_idx
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("Internal Elements Loop")
#endif

              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
                SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
                ELEMENT_USER_ELAPSED=USER_ELAPSED
                ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal interface equations assembly = ", &
                  & USER_ELAPSED, ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal interface equations assembly = ", &
                  & SYSTEM_ELAPSED,ERR,ERROR,*999)
              ENDIF
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
                SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
              ENDIF
              !Loop over the boundary and ghost elements
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("Boundary and Ghost Elements Loop")
#endif
              DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
                ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
                NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
                CALL INTERFACE_MATRICES_ELEMENT_CALCULATE(INTERFACE_MATRICES,ne,ERR,ERROR,*999)
                CALL INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ne,ERR,ERROR,*999)
                CALL INTERFACE_MATRICES_ELEMENT_ADD(INTERFACE_MATRICES,ERR,ERROR,*999)
              ENDDO !element_idx
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("Boundary and Ghost Elements Loop")
#endif
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
                SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
                ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
                ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                IF(NUMBER_OF_TIMES>0) THEN
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                    & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                    & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                ENDIF
              ENDIF
              !Finalise the element matrices
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_FINALISE()")
#endif
              CALL INTERFACE_MATRICES_ELEMENT_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_FINALISE()")
#endif
              !Output equations matrices and vector if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_MATRIX_OUTPUT) THEN
                CALL INTERFACE_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,INTERFACE_MATRICES,ERR,ERROR,*999)
              ENDIF
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Lagrange field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE_FEM")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_ASSEMBLE_FEM",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE_FEM")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_ASSEMBLE_FEM

  !
  !==================================================================================================================================
  !

  !>Finishes the process of creating an interface condition. \see OPENCMISS::CMISSInterfaceConditionCreateStart
  SUBROUTINE INTERFACE_CONDITION_CREATE_FINISH(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_idx_count,NUMBER_OF_COMPONENTS,variable_idx
    INTEGER(INTG), POINTER :: NEW_VARIABLE_MESH_INDICES(:)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_VARIABLE_PTR_TYPE), POINTER :: NEW_FIELD_VARIABLES(:)
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    NULLIFY(NEW_FIELD_VARIABLES)
    NULLIFY(NEW_VARIABLE_MESH_INDICES)
    
    CALL ENTERS("INTERFACE_CONDITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Interface condition has already been finished.",ERR,ERROR,*999)
      ELSE
        INTERFACE=>INTERFACE_CONDITION%INTERFACE
        IF(ASSOCIATED(INTERFACE)) THEN
          !Test various inputs have been set up.
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
            IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
              !Check the dependent field variables have been set.
              IF(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES<2) THEN
                LOCAL_ERROR="The number of added dependent variables of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES,"*",ERR,ERROR))// &
                  & " is invalid. The number must be >= 2."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              !\todo check if interface mesh connectivity basis has same number of gauss points as interface geometric field IF(INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY%BASIS%QUADRATURE%NUMBER_OF_GAUSS_XI/=)
              !Note There is no need to check that the dependent variables have the same number of components. 
              !The user will need to set a fixed BC on the interface dof relating to the field components 
              !not present in each of the coupled bodies, eliminating this dof from the solver matrices
              !Reorder the dependent variables based on mesh index order
              ALLOCATE(NEW_FIELD_VARIABLES(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new field variables.",ERR,ERROR,*999)
              ALLOCATE(NEW_VARIABLE_MESH_INDICES(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new variable mesh indices.",ERR,ERROR,*999)
              NEW_VARIABLE_MESH_INDICES=0
              mesh_idx_count=0
              DO mesh_idx=1,INTERFACE%NUMBER_OF_COUPLED_MESHES
                DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                  IF(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)==mesh_idx) THEN
                    mesh_idx_count=mesh_idx_count+1
                    NEW_FIELD_VARIABLES(mesh_idx_count)%PTR=>INTERFACE_DEPENDENT%FIELD_VARIABLES(variable_idx)%PTR
                    NEW_VARIABLE_MESH_INDICES(mesh_idx_count)=mesh_idx
                  ENDIF
                ENDDO !variable_idx
              ENDDO !mesh_idx
              IF(mesh_idx_count/=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES) &
                & CALL FLAG_ERROR("Invalid dependent variable mesh index setup.",ERR,ERROR,*999)
              IF(ASSOCIATED(INTERFACE_DEPENDENT%FIELD_VARIABLES)) DEALLOCATE(INTERFACE_DEPENDENT%FIELD_VARIABLES)
              IF(ASSOCIATED(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES)) DEALLOCATE(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES)
              INTERFACE_DEPENDENT%FIELD_VARIABLES=>NEW_FIELD_VARIABLES
              INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES=>NEW_VARIABLE_MESH_INDICES
            ELSE
              CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "//TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Finish the interface condition creation
          INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED=.TRUE.
        ELSE
          CALL FLAG_ERROR("Interface condition interface is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(NEW_FIELD_VARIABLES)) DEALLOCATE(NEW_FIELD_VARIABLES)
    IF(ASSOCIATED(NEW_VARIABLE_MESH_INDICES)) DEALLOCATE(NEW_VARIABLE_MESH_INDICES)
    CALL ERRORS("INTERFACE_CONDITION_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("INTERFACE_CONDITION_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_CONDITION_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating an interface condition on an interface. \see OPENCMISS::CMISSInterfaceConditionCreateStart
  SUBROUTINE INTERFACE_CONDITION_CREATE_START(USER_NUMBER,INTERFACE,GEOMETRIC_FIELD,INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the interface condition
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the interface condition on
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field for the interface condition.
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<On return, a pointer to the interface condition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,interface_conditions_idx
    TYPE(INTERFACE_TYPE), POINTER :: GEOMETRIC_INTERFACE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: NEW_INTERFACE_CONDITION
    TYPE(INTERFACE_CONDITION_PTR_TYPE), POINTER :: NEW_INTERFACE_CONDITIONS(:)
    TYPE(REGION_TYPE), POINTER :: GEOMETRIC_REGION,GEOMETRIC_INTERFACE_PARENT_REGION,INTERFACE_PARENT_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    NULLIFY(NEW_INTERFACE_CONDITION)
    NULLIFY(NEW_INTERFACE_CONDITIONS)

    CALL ENTERS("INTERFACE_CONDITION_CREATE_START",ERR,ERROR,*997)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS)) THEN
        CALL INTERFACE_CONDITION_USER_NUMBER_FIND(USER_NUMBER,INTERFACE,NEW_INTERFACE_CONDITION,ERR,ERROR,*997)
        IF(ASSOCIATED(NEW_INTERFACE_CONDITION)) THEN
          LOCAL_ERROR="Interface condition user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
        ELSE
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
            IF(GEOMETRIC_FIELD%FIELD_FINISHED) THEN
              !Check the geometric field is defined on the interface
              GEOMETRIC_INTERFACE=>GEOMETRIC_FIELD%INTERFACE
              IF(ASSOCIATED(GEOMETRIC_INTERFACE)) THEN
                IF(ASSOCIATED(GEOMETRIC_INTERFACE,INTERFACE)) THEN
                  NULLIFY(NEW_INTERFACE_CONDITION)
                  !Initialise the new interface condition
                  CALL INTERFACE_CONDITION_INITIALISE(NEW_INTERFACE_CONDITION,ERR,ERROR,*999)
                  !Set default interface condition values
                  NEW_INTERFACE_CONDITION%USER_NUMBER=USER_NUMBER
                  NEW_INTERFACE_CONDITION%GLOBAL_NUMBER=INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS+1
                  NEW_INTERFACE_CONDITION%INTERFACE_CONDITIONS=>INTERFACE%INTERFACE_CONDITIONS
                  NEW_INTERFACE_CONDITION%INTERFACE=>INTERFACE
                  !Default attributes
                  NEW_INTERFACE_CONDITION%GEOMETRY%GEOMETRIC_FIELD=>GEOMETRIC_FIELD
                  NEW_INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD
                  NEW_INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_FIELD_GAUSS_CONTINUITY_OPERATOR
                  CALL INTERFACE_CONDITION_DEPENDENT_INITIALISE(NEW_INTERFACE_CONDITION,ERR,ERROR,*999)
                  !Add new interface condition into list of interface conditions in the interface
                  ALLOCATE(NEW_INTERFACE_CONDITIONS(INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS+1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface conditions.",ERR,ERROR,*999)
                  DO interface_conditions_idx=1,INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS
                    NEW_INTERFACE_CONDITIONS(interface_conditions_idx)%PTR=>INTERFACE%INTERFACE_CONDITIONS% &
                      & INTERFACE_CONDITIONS(interface_conditions_idx)%PTR
                  ENDDO !interface_conditions_idx
                  NEW_INTERFACE_CONDITIONS(INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS+1)%PTR=> &
                    & NEW_INTERFACE_CONDITION
                  IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)) DEALLOCATE(INTERFACE%INTERFACE_CONDITIONS% &
                    & INTERFACE_CONDITIONS)
                  INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS=>NEW_INTERFACE_CONDITIONS
                  INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=INTERFACE%INTERFACE_CONDITIONS% &
                    NUMBER_OF_INTERFACE_CONDITIONS+1
                  !Return the pointer
                  INTERFACE_CONDITION=>NEW_INTERFACE_CONDITION
                ELSE
                  INTERFACE_PARENT_REGION=>INTERFACE%PARENT_REGION
                  IF(ASSOCIATED(INTERFACE_PARENT_REGION)) THEN
                    GEOMETRIC_INTERFACE_PARENT_REGION=>GEOMETRIC_INTERFACE%PARENT_REGION
                    IF(ASSOCIATED(GEOMETRIC_INTERFACE_PARENT_REGION)) THEN
                      LOCAL_ERROR="Geometric field interface does not match specified interface. "// &
                        "The geometric field was created on interface number "// &
                        & TRIM(NUMBER_TO_VSTRING(GEOMETRIC_INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
                        & " of parent region number "// &
                        & TRIM(NUMBER_TO_VSTRING(GEOMETRIC_INTERFACE_PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                        & " and the specified interface was created as number "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" on parent region number "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE_PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Geometric interface parent region is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface parent region is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                GEOMETRIC_REGION=>GEOMETRIC_FIELD%REGION
                IF(ASSOCIATED(GEOMETRIC_REGION)) THEN
                  LOCAL_ERROR="The geometric field was created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(GEOMETRIC_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and not on the specified interface."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("The geometric field does not have a region or interface created.",ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Geometric field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Geometric field is not finished.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The interface conditions on interface number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*997)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_INTERFACE_CONDITION)) CALL INTERFACE_CONDITION_FINALISE(NEW_INTERFACE_CONDITION,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_INTERFACE_CONDITIONS)) DEALLOCATE(NEW_INTERFACE_CONDITIONS)
997 CALL ERRORS("INTERFACE_CONDITION_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_CREATE_START")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Reprojecting contact points for frictionless contact
  SUBROUTINE INTERFACE_CONDITION_DATA_REPROJECTION(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CONNECTIVITY
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection to store the xi locations and element number for the data points
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_1,DEPENDENT_FIELD_PROJECTION!< Dependent field to interpolate
    TYPE(DATA_PROJECTION_RESULT_TYPE), POINTER :: DATA_PROJECTION_RESULT
    INTEGER(INTG) :: data_point_idx,coupled_mesh_idx,xi_idx,MESH_COMPONENT_NUMBER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
 
    CALL ENTERS("INTERFACE_CONDITION_DATA_REPROJECTION",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      POINTS_CONNECTIVITY=>INTERFACE_CONDITION%INTERFACE%POINTS_CONNECTIVITY
      IF(ASSOCIATED(POINTS_CONNECTIVITY)) THEN
        DATA_POINTS=>INTERFACE_CONDITION%INTERFACE%DATA_POINTS
        !TODO: Default first region data projection to be unchanged throughout the simulation, i.e. data embedded in the first region mesh
        DEPENDENT_FIELD_1=>INTERFACE_CONDITION%DEPENDENT%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
        !Evaluate the data points positions according to the xi locations of first mesh and its dependent field.
        CALL DATA_POINTS_VALUES_POINTS_CONNECTIVITY_FIELD_EVALUATE(DATA_POINTS,POINTS_CONNECTIVITY,1,DEPENDENT_FIELD_1, &
          & ERR,ERROR,*999)
        !Get the mesh component number for this field
        !MESH_COMPONENT_NUMBER=DEPENDENT_FIELD_1%DECOMPOSITION%MESH_COMPONENT_NUMBER
        MESH_COMPONENT_NUMBER=1
        DO coupled_mesh_idx=2,INTERFACE_CONDITION%INTERFACE%NUMBER_OF_COUPLED_MESHES
          DATA_PROJECTION=>DATA_POINTS%DATA_PROJECTIONS(coupled_mesh_idx)%PTR
          DEPENDENT_FIELD_PROJECTION=>INTERFACE_CONDITION%DEPENDENT%EQUATIONS_SETS(coupled_mesh_idx)%PTR%DEPENDENT%DEPENDENT_FIELD
          !Projection the data points (with know spatial positions) on the dependent field of the second (or more) mesh
          CALL DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE(DATA_PROJECTION,DEPENDENT_FIELD_PROJECTION,err,error,*999)
          !Loop through data points and record the new element number and xi location of the data points in the second mesh
          DO data_point_idx=1,POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS
            DATA_PROJECTION_RESULT=>DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)
            POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%COUPLED_MESH_ELEMENT_NUMBER= &
              & DATA_PROJECTION_RESULT%ELEMENT_NUMBER
            IF(DATA_PROJECTION_RESULT%ELEMENT_LINE_NUMBER/=0) THEN
              POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%COUPLED_MESH_CONTACT_NUMBER= &
                & DATA_PROJECTION_RESULT%ELEMENT_LINE_NUMBER
              SELECT CASE(ABS(POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)% &
                & COUPLED_MESH_CONTACT_XI_NORMAL))
              CASE(1)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(2,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(1)
              CASE(2)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(1,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(1)
              END SELECT
            ELSEIF(DATA_PROJECTION_RESULT%ELEMENT_FACE_NUMBER/=0) THEN
              POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%COUPLED_MESH_CONTACT_NUMBER= &
                & DATA_PROJECTION_RESULT%ELEMENT_FACE_NUMBER
              SELECT CASE(ABS(POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)% &
                & COUPLED_MESH_CONTACT_XI_NORMAL))
              CASE(1)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(2,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(1)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(3,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(2)
              CASE(2)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(1,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(1)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(3,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(2)
              CASE(3)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(1,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(1)
                POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(data_point_idx,coupled_mesh_idx)%XI(2,MESH_COMPONENT_NUMBER)= &
                  & DATA_PROJECTION_RESULT%XI(2)
              END SELECT
            ENDIF
          ENDDO !data_point_idx      
        ENDDO
        POINTS_CONNECTIVITY%DATA_POINT_PROJECTED=.TRUE.
        ! Check if the data has been orthogonally projected
        DO data_point_idx=1,POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS
          DO xi_idx=1,SIZE(DATA_POINTS%DATA_PROJECTIONS(2)%PTR%DATA_PROJECTION_RESULTS(data_point_idx)%XI)
              IF(DATA_POINTS%DATA_PROJECTIONS(2)%PTR%DATA_PROJECTION_RESULTS(data_point_idx)%XI(xi_idx)==0.0_DP .OR.  &
                  & DATA_POINTS%DATA_PROJECTIONS(2)%PTR%DATA_PROJECTION_RESULTS(data_point_idx)%XI(xi_idx)==1.0_DP) THEN
                POINTS_CONNECTIVITY%DATA_POINT_PROJECTED(data_point_idx)=.FALSE.
              ENDIF
            ENDDO !xi_idx
        ENDDO !data_point_idx
      ELSE
        CALL FLAG_ERROR("Interface condition points connectivity is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_DATA_REPROJECTION")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_DATA_REPROJECTION",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_DATA_REPROJECTION")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_DATA_REPROJECTION
  
  !
  !================================================================================================================================
  !

  !>Reprojecting contact points for frictionless contact
  SUBROUTINE INTERFACE_CONDITION_DATA_POINTS_NORMAL_CALCULATE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CONNECTIVITY
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection to store the xi locations and element number for the data points
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_PROJECTION!< Dependent field to interpolate
    TYPE(DATA_PROJECTION_RESULT_TYPE), POINTER :: DATA_PROJECTION_RESULT
     TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: INTERPOLATED_POINT_METRICS
    INTEGER(INTG) :: dataPointIdx,projectionMeshIdx,xi_idx,MESH_COMPONENT_NUMBER,coupledMeshNumberOfXi,localContactXiNormal, &
      & globalContactNumber,coupledMeshElementNumber,localContactNumber,numberOfDimensions
    REAL(DP) :: XI_POSITION(2),XI(3),POSITION(3),NORMAL(3),TANGENTS(3,3)
    LOGICAL :: REVERSE_NORMAL
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
 
    CALL ENTERS("INTERFACE_CONDITION_DATA_POINTS_NORMAL_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      POINTS_CONNECTIVITY=>INTERFACE_CONDITION%INTERFACE%POINTS_CONNECTIVITY
      IF(ASSOCIATED(POINTS_CONNECTIVITY)) THEN
        DATA_POINTS=>INTERFACE_CONDITION%INTERFACE%DATA_POINTS
        !TODO: Default normal to be on second mesh surface 
        DEPENDENT_FIELD_PROJECTION=>INTERFACE_CONDITION%DEPENDENT%EQUATIONS_SETS(2)%PTR%DEPENDENT%DEPENDENT_FIELD
        MESH_COMPONENT_NUMBER=1
        coupledMeshNumberOfXi=INTERFACE_CONDITION%INTERFACE%MESHES%MESHES(1)%PTR%NUMBER_OF_DIMENSIONS+1
        projectionMeshIdx=2
        numberOfDimensions=INTERFACE_CONDITION%INTERFACE%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
        INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
        DO dataPointIdx=1,POINTS_CONNECTIVITY%NUMBER_OF_DATA_POINTS
          IF(POINTS_CONNECTIVITY%DATA_POINT_PROJECTED(dataPointIdx)) THEN
            !The xi direction of the surface/line normal
            localContactXiNormal=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointIdx,projectionMeshIdx)% &
              & COUPLED_MESH_CONTACT_XI_NORMAL 
            coupledMeshElementNumber=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointIdx,projectionMeshIdx)% &
              & COUPLED_MESH_ELEMENT_NUMBER
            !Local number of the line in contact
            localContactNumber=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointIdx,projectionMeshIdx)% &
              & COUPLED_MESH_CONTACT_NUMBER
            !Decide if the normal vector should be reversed to be outward  
            IF(ABS(localContactXiNormal)==localContactXiNormal) THEN
              REVERSE_NORMAL=.FALSE.
            ELSE
              REVERSE_NORMAL=.TRUE. 
            ENDIF   
            
            XI(1:coupledMeshNumberOfXi)=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointIdx, &
              & projectionMeshIdx)%XI(1:coupledMeshNumberOfXi,MESH_COMPONENT_NUMBER)   
            
            
            SELECT CASE(coupledMeshNumberOfXi)
            CASE(2)
              !Global number of the line in contact
              globalContactNumber=POINTS_CONNECTIVITY%INTERFACE%COUPLED_MESHES(projectionMeshIdx)%PTR% &
                & DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)% &
                & ELEMENT_LINES(localContactNumber)   
              !Compute interpolation parameters for the line for the couple mesh field
              CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,globalContactNumber, &
                & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)% &
                & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)!Needed for computing metrics
              !Compute the other xi direction of the line normal, i.e. xi position input for calculating interpolated point
              XI_POSITION(1)=XI(OTHER_XI_DIRECTIONS2(ABS(localContactXiNormal)))
            CASE(3)
              !Global number of the line in contact
              globalContactNumber=POINTS_CONNECTIVITY%INTERFACE%COUPLED_MESHES(projectionMeshIdx)%PTR% &
                & DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)% &
                & ELEMENT_FACES(localContactNumber)   
              !Compute interpolation parameters for the face for the couple mesh field
              CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,globalContactNumber, &
                & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)% &
                & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              SELECT CASE(ABS(localContactXiNormal))
              CASE(1)
              !Compute the other xi directions of the face normal, i.e. xi positions input for calculating interpolated point
                XI_POSITION(1)=XI(2)
                XI_POSITION(2)=XI(3)
              CASE(2)
                XI_POSITION(1)=XI(1)
                XI_POSITION(2)=XI(3)
              CASE(3)
                XI_POSITION(1)=XI(1)
                XI_POSITION(2)=XI(2)
              END SELECT
            CASE DEFAULT
                CALL FLAG_ERROR("Xi dimension of coupled mesh should be <= 3 and >=0.",ERR,ERROR,*999)
            END SELECT
            CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,XI_POSITION(1:coupledMeshNumberOfXi-1),INTERFACE_EQUATIONS% &
              & INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
              & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)!Needed for computing metrics
            ALLOCATE(INTERFACE_EQUATIONS% &
              & INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
              & INTERPOLATED_POINT_METRICS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated points metrics.",ERR,ERROR,*999)
            NULLIFY(INTERFACE_EQUATIONS% &
              & INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
              & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR)
            CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(INTERFACE_EQUATIONS% &
              & INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
              & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,INTERFACE_EQUATIONS% &
              & INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
              & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(coupledMeshNumberOfXi,INTERFACE_EQUATIONS% &
              & INTERPOLATION%VARIABLE_INTERPOLATION(projectionMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
              & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)                      
            INTERPOLATED_POINT_METRICS=>INTERFACE_EQUATIONS%INTERPOLATION% &
              & VARIABLE_INTERPOLATION(projectionMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
              & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR                 
            !Calculate outward normal of the contact surface in the coupled mesh
            CALL FIELD_POSITION_NORMAL_TANGENTS_CALCULATE_INT_PT_METRIC(INTERPOLATED_POINT_METRICS, &
              & REVERSE_NORMAL,POSITION(1:numberOfDimensions), &
              & NORMAL(1:numberOfDimensions),TANGENTS(1:numberOfDimensions, &
              & 1:coupledMeshNumberOfXi),ERR,ERROR,*999)  
            POINTS_CONNECTIVITY%NORMAL(1:numberOfDimensions,dataPointIdx)=NORMAL

            
          ENDIF
        ENDDO !dataPointIdx
      ELSE
        CALL FLAG_ERROR("Interface condition points connectivity is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_DATA_POINTS_NORMAL_CALCULATE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_DATA_POINTS_NORMAL_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_DATA_POINTS_NORMAL_CALCULATE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_DATA_POINTS_NORMAL_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition dependent field information and deallocate all memory.
  SUBROUTINE INTERFACE_CONDITION_DEPENDENT_FINALISE(INTERFACE_DEPENDENT,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT !<A pointer to the interface condition dependent field information to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_DEPENDENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
      IF(ASSOCIATED(INTERFACE_DEPENDENT%EQUATIONS_SETS)) DEALLOCATE(INTERFACE_DEPENDENT%EQUATIONS_SETS)
      IF(ASSOCIATED(INTERFACE_DEPENDENT%FIELD_VARIABLES)) DEALLOCATE(INTERFACE_DEPENDENT%FIELD_VARIABLES)
      IF(ASSOCIATED(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES)) DEALLOCATE(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES)
      DEALLOCATE(INTERFACE_DEPENDENT)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_DEPENDENT_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_DEPENDENT_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_DEPENDENT_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_DEPENDENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition dependent field information.
  SUBROUTINE INTERFACE_CONDITION_DEPENDENT_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<The pointer to the interface condition to initialise to initialise the dependent field information for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_DEPENDENT_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%DEPENDENT)) THEN
        CALL FLAG_ERROR("Interface condition dependent is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE_CONDITION%DEPENDENT,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface condition dependent.",ERR,ERROR,*999)
        INTERFACE_CONDITION%DEPENDENT%INTERFACE_CONDITION=>INTERFACE_CONDITION
        INTERFACE_CONDITION%DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES=0
        NULLIFY(INTERFACE_CONDITION%DEPENDENT%EQUATIONS_SETS)
        NULLIFY(INTERFACE_CONDITION%DEPENDENT%FIELD_VARIABLES)
        NULLIFY(INTERFACE_CONDITION%DEPENDENT%VARIABLE_MESH_INDICES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_DEPENDENT_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITION_DEPENDENT_FINALISE(INTERFACE_CONDITION%DEPENDENT,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITION_DEPENDENT_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_DEPENDENT_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_DEPENDENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds an equations set to an interface condition. \see OPENCMISS::CMISSInterfaceConditionEquationsSetAdd
  SUBROUTINE INTERFACE_CONDITION_DEPENDENT_VARIABLE_ADD(INTERFACE_CONDITION,MESH_INDEX,EQUATIONS_SET,VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to add the dependent variable to
    INTEGER(INTG), INTENT(IN) :: MESH_INDEX !<The mesh index in the interface conditions interface that the dependent variable corresponds to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set containing the dependent field to add the variable from.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The variable type of the dependent field to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    INTEGER(INTG), POINTER :: NEW_VARIABLE_MESH_INDICES(:)
    LOGICAL :: FOUND_MESH_INDEX
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,INTERFACE_VARIABLE
    TYPE(FIELD_VARIABLE_PTR_TYPE), POINTER :: NEW_FIELD_VARIABLES(:)
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(MESH_TYPE), POINTER :: DEPENDENT_MESH,INTERFACE_MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_DEPENDENT_VARIABLE_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
      IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
        INTERFACE=>INTERFACE_CONDITION%INTERFACE
        IF(ASSOCIATED(INTERFACE)) THEN
          IF(MESH_INDEX>0.AND.MESH_INDEX<=INTERFACE%NUMBER_OF_COUPLED_MESHES) THEN
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                  FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                    !Check that the field variable hasn't already been added.
                    variable_idx=1
                    NULLIFY(INTERFACE_VARIABLE)                   
                    DO WHILE(variable_idx<=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES.AND. &
                      & .NOT.ASSOCIATED(INTERFACE_VARIABLE))
                      !To check if the field variable is the variable associated with the interface
                      IF(ASSOCIATED(FIELD_VARIABLE,INTERFACE_DEPENDENT%FIELD_VARIABLES(variable_idx)%PTR)) THEN
                        INTERFACE_VARIABLE=>INTERFACE_DEPENDENT%FIELD_VARIABLES(variable_idx)%PTR
                      ELSE
                        variable_idx=variable_idx+1
                      ENDIF
                    ENDDO
                    IF(ASSOCIATED(INTERFACE_VARIABLE)) THEN
                      !Check if we are dealing with the same mesh index.
                      IF(MESH_INDEX/=INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)) THEN
                        LOCAL_ERROR="The dependent variable has already been added to the interface condition at "// &
                          & "position index "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      !Check the dependent variable and the mesh index match.
                      INTERFACE_MESH=>INTERFACE%COUPLED_MESHES(MESH_INDEX)%PTR
                      IF(ASSOCIATED(INTERFACE_MESH)) THEN
                        DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
                        IF(ASSOCIATED(DECOMPOSITION)) THEN
                          DEPENDENT_MESH=>DECOMPOSITION%MESH
                          IF(ASSOCIATED(DEPENDENT_MESH)) THEN
                            IF(ASSOCIATED(INTERFACE_MESH,DEPENDENT_MESH)) THEN
                              !The meshes match. Check if the dependent variable has already been added for the mesh index.
                              FOUND_MESH_INDEX=.FALSE.
                              DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                                IF(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)==MESH_INDEX) THEN
                                  FOUND_MESH_INDEX=.TRUE.
                                  EXIT
                                ENDIF
                              ENDDO !variable_idx
                              IF(FOUND_MESH_INDEX) THEN
                                !The mesh index has already been added to replace the dependent variable with the specified variable
                                INTERFACE_DEPENDENT%FIELD_VARIABLES(variable_idx)%PTR=>DEPENDENT_FIELD% &
                                  & VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
                              ELSE
                                !The mesh index has not been found so add a new dependent variable.
                                ALLOCATE(NEW_EQUATIONS_SETS(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations sets.",ERR,ERROR,*999)
                                ALLOCATE(NEW_FIELD_VARIABLES(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new field variables.",ERR,ERROR,*999)
                                ALLOCATE(NEW_VARIABLE_MESH_INDICES(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new variable mesh indices.",ERR,ERROR,*999)
                                DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                                  NEW_EQUATIONS_SETS(variable_idx)%PTR=>INTERFACE_DEPENDENT%EQUATIONS_SETS(variable_idx)%PTR
                                  NEW_FIELD_VARIABLES(variable_idx)%PTR=>INTERFACE_DEPENDENT%FIELD_VARIABLES(variable_idx)%PTR
                                  NEW_VARIABLE_MESH_INDICES(variable_idx)=INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES(variable_idx)
                                ENDDO !variable_idx
                                NEW_EQUATIONS_SETS(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1)%PTR=>EQUATIONS_SET
                                NEW_FIELD_VARIABLES(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1)%PTR=>DEPENDENT_FIELD% &
                                  & VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
                                NEW_VARIABLE_MESH_INDICES(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1)=MESH_INDEX
                                IF(ASSOCIATED(INTERFACE_DEPENDENT%EQUATIONS_SETS)) DEALLOCATE(INTERFACE_DEPENDENT%EQUATIONS_SETS)
                                IF(ASSOCIATED(INTERFACE_DEPENDENT%FIELD_VARIABLES)) DEALLOCATE(INTERFACE_DEPENDENT%FIELD_VARIABLES)
                                IF(ASSOCIATED(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES)) &
                                  & DEALLOCATE(INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES)
                                INTERFACE_DEPENDENT%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
                                INTERFACE_DEPENDENT%FIELD_VARIABLES=>NEW_FIELD_VARIABLES
                                INTERFACE_DEPENDENT%VARIABLE_MESH_INDICES=>NEW_VARIABLE_MESH_INDICES
                                INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES= &
                                  & INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("The dependent field mesh does not match the interface mesh.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("The dependent field decomposition mesh is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("The dependent field decomposition is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="The interface mesh for mesh index "//TRIM(NUMBER_TO_VSTRING(MESH_INDEX,"*",ERR,ERROR))// &
                          & " is not associated."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " has not been created on field number "// &
                      & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The variable type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specificed mesh index of "//TRIM(NUMBER_TO_VSTRING(MESH_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The mesh index must be > 0 and <= "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE%NUMBER_OF_COUPLED_MESHES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface condition interface is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_DEPENDENT_VARIABLE_ADD")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_DEPENDENT_VARIABLE_ADD",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_DEPENDENT_VARIABLE_ADD")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_DEPENDENT_VARIABLE_ADD
  
  !
  !================================================================================================================================
  !

  !>Destroys an interface condition. \see OPENCMISS::CMISSInterfaceConditionDestroy
  SUBROUTINE INTERFACE_CONDITION_DESTROY(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_condition_idx,interface_condition_position
    TYPE(INTERFACE_CONDITION_PTR_TYPE), POINTER :: NEW_INTERFACE_CONDITIONS(:)
    TYPE(INTERFACE_CONDITIONS_TYPE), POINTER :: INTERFACE_CONDITIONS

    NULLIFY(NEW_INTERFACE_CONDITIONS)

    CALL ENTERS("INTERFACE_CONDITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_CONDITIONS=>INTERFACE_CONDITION%INTERFACE_CONDITIONS
      IF(ASSOCIATED(INTERFACE_CONDITIONS)) THEN
        interface_condition_position=INTERFACE_CONDITION%GLOBAL_NUMBER

        !Destroy all the interface condition components
        CALL INTERFACE_CONDITION_FINALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
        
        !Remove the interface condition from the list of interface conditions
        IF(INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS>1) THEN
          ALLOCATE(NEW_INTERFACE_CONDITIONS(INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface conditions.",ERR,ERROR,*999)
          DO interface_condition_idx=1,INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS
            IF(interface_condition_idx<interface_condition_position) THEN
              NEW_INTERFACE_CONDITIONS(interface_condition_idx)%PTR=>INTERFACE_CONDITIONS% &
                & INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            ELSE IF(interface_condition_idx>interface_condition_position) THEN
              INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interface_condition_idx)%PTR%GLOBAL_NUMBER=INTERFACE_CONDITIONS% &
                & INTERFACE_CONDITIONS(interface_condition_idx)%PTR%GLOBAL_NUMBER-1
              NEW_INTERFACE_CONDITIONS(interface_condition_idx-1)%PTR=>INTERFACE_CONDITIONS% &
                & INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            ENDIF
          ENDDO !interface_conditions_idx
          IF(ASSOCIATED(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)) DEALLOCATE(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)
          INTERFACE_CONDITIONS%INTERFACE_CONDITIONS=>NEW_INTERFACE_CONDITIONS
          INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS-1
        ELSE
          DEALLOCATE(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)
          INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Interface conditions interface conditions is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("INTERFACE_CONDITIONS_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_INTERFACE_CONDITIONS)) DEALLOCATE(NEW_INTERFACE_CONDITIONS)
998 CALL ERRORS("INTERFACE_CONDITION_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_DESTROY")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of interface equations for the interface condition. \see OPENCMISS::CMISSInterfaceConditionEquationsCreateFinish
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finish the creation of the interface equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: STORAGE_TYPE(:),STRUCTURE_TYPE(:)
    LOGICAL, ALLOCATABLE :: MATRICES_TRANSPOSE(:)
    INTEGER(INTG) :: number_of_dependent_variables
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_CONDITIONS_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      SELECT CASE(INTERFACE_CONDITION%METHOD)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
        !Finish the interface equations creation
        NULLIFY(INTERFACE_EQUATIONS)
        CALL INTERFACE_CONDITION_EQUATIONS_GET(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*999)
        IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
          CALL FLAG_ERROR("Interface condition equations have already been finished.",ERR,ERROR,*999)
        ELSE
          CALL INTERFACE_EQUATIONS_CREATE_FINISH(INTERFACE_EQUATIONS,ERR,ERROR,*999)
          INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
          IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
            !Create the interface mapping.
            NULLIFY(INTERFACE_MAPPING)
            CALL INTERFACE_MAPPING_CREATE_START(INTERFACE_EQUATIONS,INTERFACE_MAPPING,ERR,ERROR,*999)
            CALL INTERFACE_MAPPING_LAGRANGE_VARIABLE_TYPE_SET(INTERFACE_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
            SELECT CASE(INTERFACE_CONDITION%METHOD)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
              number_of_dependent_variables=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
            CASE(INTERFACE_CONDITION_PENALTY_METHOD)
              number_of_dependent_variables=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES+1
            ENDSELECT
            CALL INTERFACE_MAPPING_MATRICES_NUMBER_SET(INTERFACE_MAPPING,number_of_dependent_variables,ERR,ERROR,*999)
            ALLOCATE(MATRICES_TRANSPOSE(number_of_dependent_variables),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrices transpose.",ERR,ERROR,*999)
            MATRICES_TRANSPOSE=.TRUE.
            SELECT CASE(INTERFACE_CONDITION%METHOD)
            CASE(INTERFACE_CONDITION_PENALTY_METHOD)
              !Set the last interface matrix to have no transpose
              MATRICES_TRANSPOSE(number_of_dependent_variables)=.FALSE.
            ENDSELECT
            CALL INTERFACE_MAPPING_MATRICES_TRANSPOSE_SET(INTERFACE_MAPPING,MATRICES_TRANSPOSE,ERR,ERROR,*999)
            IF(ALLOCATED(MATRICES_TRANSPOSE)) DEALLOCATE(MATRICES_TRANSPOSE)
            CALL INTERFACE_MAPPING_RHS_VARIABLE_TYPE_SET(INTERFACE_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
            CALL INTERFACE_MAPPING_CREATE_FINISH(INTERFACE_MAPPING,ERR,ERROR,*999)
            !Create the interface matrices
            NULLIFY(INTERFACE_MATRICES)
            CALL INTERFACE_MATRICES_CREATE_START(INTERFACE_EQUATIONS,INTERFACE_MATRICES,ERR,ERROR,*999)
            ALLOCATE(STORAGE_TYPE(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate storage type.",ERR,ERROR,*999)
            SELECT CASE(INTERFACE_EQUATIONS%SPARSITY_TYPE)
            CASE(INTERFACE_MATRICES_FULL_MATRICES) 
              STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
              CALL INTERFACE_MATRICES_STORAGE_TYPE_SET(INTERFACE_MATRICES,STORAGE_TYPE,ERR,ERROR,*999)
            CASE(INTERFACE_MATRICES_SPARSE_MATRICES) 
              ALLOCATE(STRUCTURE_TYPE(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate structure type.",ERR,ERROR,*999)
              STORAGE_TYPE=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
              STRUCTURE_TYPE=INTERFACE_MATRIX_FEM_STRUCTURE
              CALL INTERFACE_MATRICES_STORAGE_TYPE_SET(INTERFACE_MATRICES,STORAGE_TYPE,ERR,ERROR,*999)
              CALL INTERFACE_MATRICES_STRUCTURE_TYPE_SET(INTERFACE_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*999)
              IF(ALLOCATED(STRUCTURE_TYPE)) DEALLOCATE(STRUCTURE_TYPE)
            CASE DEFAULT
              LOCAL_ERROR="The interface equations sparsity type of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            IF(ALLOCATED(STORAGE_TYPE)) DEALLOCATE(STORAGE_TYPE)
            CALL INTERFACE_MATRICES_CREATE_FINISH(INTERFACE_MATRICES,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The interface condition method of "//TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH")
    RETURN
999 IF(ALLOCATED(MATRICES_TRANSPOSE)) DEALLOCATE(MATRICES_TRANSPOSE)
    IF(ALLOCATED(STORAGE_TYPE)) DEALLOCATE(STORAGE_TYPE)
    IF(ALLOCATED(STRUCTURE_TYPE)) DEALLOCATE(STRUCTURE_TYPE)
    CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of interface equations for the interface condition. \see CMISSInterfaceConditionEquationsCreateStart
  !>Default values set for the INTERFACE_EQUATIONS's attributes are:
  !>- OUTPUT_TYPE: 0 (INTERFACE_EQUATIONS_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (INTERFACE_EQUATIONS_SPARSE_MATRICES)
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_START(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to create the interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<On exit, a pointer to the created interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        CALL FLAG_ERROR("Interface equations is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(INTERFACE_EQUATIONS)
        SELECT CASE(INTERFACE_CONDITION%METHOD)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
            IF(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FINISHED) THEN
              INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
              IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                !Initialise the setup
                CALL INTERFACE_EQUATIONS_CREATE_START(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*999)
                !Set the number of interpolation sets
                CALL INTERFACE_EQUATIONS_INTERFACE_INTERP_SETS_NUMBER_SET(INTERFACE_EQUATIONS,1,1,1,ERR,ERROR,*999)
                DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                  CALL INTERFACE_EQUATIONS_VARIABLE_INTERP_SETS_NUMBER_SET(INTERFACE_EQUATIONS,variable_idx,1,1,0, &
                    & ERR,ERROR,*999)
                ENDDO !variable_idx
              ELSE
                CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
              ENDIF
              !Return the pointer
              INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
            ELSE
              CALL FLAG_ERROR("Interface condition Lagrange field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The interface condition method of "//TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the interface equations for an interface condition. \see OPENCMISS::CMISSInterfaceConditionEquationsDestroy
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_DESTROY(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface conditions to destroy the interface equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) THEN
        CALL INTERFACE_EQUATIONS_DESTROY(INTERFACE_CONDITION%INTERFACE_EQUATIONS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_DESTROY")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the interface condition and deallocate all memory.
  SUBROUTINE INTERFACE_CONDITION_FINALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      CALL INTERFACE_CONDITION_GEOMETRY_FINALISE(INTERFACE_CONDITION%GEOMETRY,ERR,ERROR,*999)
      CALL INTERFACE_CONDITION_LAGRANGE_FINALISE(INTERFACE_CONDITION%LAGRANGE,ERR,ERROR,*999)
      CALL INTERFACE_CONDITION_PENALTY_FINALISE(INTERFACE_CONDITION%PENALTY,ERR,ERROR,*999)
      CALL INTERFACE_CONDITION_DEPENDENT_FINALISE(INTERFACE_CONDITION%DEPENDENT,ERR,ERROR,*999)
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) &
        & CALL INTERFACE_EQUATIONS_DESTROY(INTERFACE_CONDITION%INTERFACE_EQUATIONS,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_CONDITION)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_FINALISE

  !
  !================================================================================================================================
  !


  !
  !================================================================================================================================
  !

  !>Finalise the interface condition geometry information and deallocate all memory.
  SUBROUTINE INTERFACE_CONDITION_GEOMETRY_FINALISE(INTERFACE_GEOMETRY,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_GEOMETRY_TYPE) :: INTERFACE_GEOMETRY !<The interface condition geometry information to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_GEOMETRY_FINALISE",ERR,ERROR,*999)

    NULLIFY(INTERFACE_GEOMETRY%INTERFACE_CONDITION)
    NULLIFY(INTERFACE_GEOMETRY%GEOMETRIC_FIELD)
       
    CALL EXITS("INTERFACE_CONDITION_GEOMETRY_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_GEOMETRY_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_GEOMETRY_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_GEOMETRY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition geometry information.
  SUBROUTINE INTERFACE_CONDITION_GEOMETRY_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<The pointer to the interface condition to initialise to initialise the geometry information for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_GEOMETRY_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_CONDITION%GEOMETRY%INTERFACE_CONDITION=>INTERFACE_CONDITION
      NULLIFY(INTERFACE_CONDITION%GEOMETRY%GEOMETRIC_FIELD)
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_GEOMETRY_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITION_GEOMETRY_FINALISE(INTERFACE_CONDITION%GEOMETRY,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITION_GEOMETRY_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_GEOMETRY_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_GEOMETRY_INITIALISE

   !
  !================================================================================================================================
  !

  !>Initialises an interface condition.
  SUBROUTINE INTERFACE_CONDITION_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<The pointer to the interface condition to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      CALL FLAG_ERROR("Interface condition is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(INTERFACE_CONDITION,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface condition.",ERR,ERROR,*999)
      INTERFACE_CONDITION%USER_NUMBER=0
      INTERFACE_CONDITION%GLOBAL_NUMBER=0
      INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED=.FALSE.
      NULLIFY(INTERFACE_CONDITION%INTERFACE_CONDITIONS)
      NULLIFY(INTERFACE_CONDITION%INTERFACE)
      INTERFACE_CONDITION%METHOD=0
      INTERFACE_CONDITION%OPERATOR=0
      NULLIFY(INTERFACE_CONDITION%LAGRANGE)
      NULLIFY(INTERFACE_CONDITION%PENALTY)
      NULLIFY(INTERFACE_CONDITION%DEPENDENT)
      NULLIFY(INTERFACE_CONDITION%INTERFACE_EQUATIONS)
      CALL INTERFACE_CONDITION_GEOMETRY_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
      NULLIFY(INTERFACE_CONDITION%BOUNDARY_CONDITIONS)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITION_FINALISE(INTERFACE_CONDITION,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITION_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an interface condition's Lagrange multiplier field \see OPENCMISS::CMISSInterfaceConditionLagrangeConditionCreateFinish
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finish creating the Lagrange field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
        IF(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FINISHED) THEN
          CALL FLAG_ERROR("Interface condition Lagrange field has already been finished.",ERR,ERROR,*999)
        ELSE
          !Finish the Lagrange field creation
          IF(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD_AUTO_CREATED) THEN
            CALL FIELD_CREATE_FINISH(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,ERR,ERROR,*999)
          ENDIF
          INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating the Lagrange multiplyer field for interface condition. \see OPENCMISS::CMISSInterfaceConditionLagrangeFieldCreateStart
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START(INTERFACE_CONDITION,LAGRANGE_FIELD_USER_NUMBER,LAGRANGE_FIELD, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to create the Lagrange field on
    INTEGER(INTG), INTENT(IN) :: LAGRANGE_FIELD_USER_NUMBER !<The user specified Lagrange field number
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD !<If associated on entry, a pointer to the user created Lagrange field which has the same user number as the specified Lagrange field user number. If not associated on entry, on exit, a pointer to the created Lagrange field for the interface condition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,interpolation_type,GEOMETRIC_SCALING_TYPE,dependent_variable_number
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(REGION_TYPE), POINTER :: INTERFACE_REGION,LAGRANGE_FIELD_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
        CALL FLAG_ERROR("Interface condition Lagrange is already associated.",ERR,ERROR,*999)
      ELSE
        INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
        IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
          INTERFACE=>INTERFACE_CONDITION%INTERFACE
          IF(ASSOCIATED(INTERFACE)) THEN
            INTERFACE_REGION=>INTERFACE%PARENT_REGION
            IF(ASSOCIATED(INTERFACE_REGION)) THEN
              IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                !Check the Lagrange field has been finished
                IF(LAGRANGE_FIELD%FIELD_FINISHED) THEN
                  !Check the user numbers match
                  IF(LAGRANGE_FIELD_USER_NUMBER/=LAGRANGE_FIELD%USER_NUMBER) THEN
                    LOCAL_ERROR="The specified Lagrange field user number of "// &
                      & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                      & " does not match the user number of the specified Lagrange field of "// &
                      & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  LAGRANGE_FIELD_REGION=>LAGRANGE_FIELD%REGION
                  IF(ASSOCIATED(LAGRANGE_FIELD_REGION)) THEN
                    !Check the field is defined on the same region as the interface
                    IF(LAGRANGE_FIELD_REGION%USER_NUMBER/=INTERFACE_REGION%USER_NUMBER) THEN
                      LOCAL_ERROR="Invalid region setup. The specified Lagrange field has been created on interface number "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" in parent region number "// &
                        & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                        & " and the specified interface has been created in parent region number "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The Lagrange field region is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The specified Lagrange field has not been finished.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !Check the user number has not already been used for a field in this region.
                NULLIFY(FIELD)
                CALL FIELD_USER_NUMBER_FIND(LAGRANGE_FIELD_USER_NUMBER,INTERFACE,FIELD,ERR,ERROR,*999)
                IF(ASSOCIATED(FIELD)) THEN
                  LOCAL_ERROR="The specified Lagrange field user number of "// &
                    & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                    & " has already been used to create a field on interface number "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
              CALL INTERFACE_CONDITION_LAGRANGE_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
              IF(.NOT.ASSOCIATED(LAGRANGE_FIELD)) THEN
                !Create the Lagrange field
                INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD_AUTO_CREATED=.TRUE.
                CALL FIELD_CREATE_START(LAGRANGE_FIELD_USER_NUMBER,INTERFACE_CONDITION%INTERFACE,INTERFACE_CONDITION%LAGRANGE% &
                  & LAGRANGE_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,"Lagrange Multipliers Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_DEPENDENT_TYPE, &
                  & ERR,ERROR,*999)
                NULLIFY(GEOMETRIC_DECOMPOSITION)
                CALL FIELD_MESH_DECOMPOSITION_GET(INTERFACE_CONDITION%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,INTERFACE_CONDITION%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)              
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,2,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_U_VARIABLE_TYPE,"Lambda", &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & "Lambda RHS",ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                !Note that only components present in both the coupled meshes interface dependent fields can be coupled
                !Default the number of component to be the minimum number of components across all the coupled dependent variables
                !\todo Check ordering of variable components which are coupled and uncoupled are handled correctly to ensure that
                !coupled variable components don't have to always come before the uncoupled variable components
                INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS=0
                DO dependent_variable_number=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                  IF (INTERFACE_DEPENDENT%FIELD_VARIABLES(dependent_variable_number)%PTR%NUMBER_OF_COMPONENTS< &
                    & INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS) THEN
                    INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS= &
                      & INTERFACE_DEPENDENT%FIELD_VARIABLES(dependent_variable_number)%PTR%NUMBER_OF_COMPONENTS-1
                  ELSEIF (INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS==0) THEN
                    INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS= &
                      & INTERFACE_DEPENDENT%FIELD_VARIABLES(dependent_variable_number)%PTR%NUMBER_OF_COMPONENTS-1
                  ENDIF
                ENDDO
                CALL FIELD_NUMBER_OF_COMPONENTS_SET(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)                
                DO component_idx=1,INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS
                  !CALL FIELD_COMPONENT_INTERPOLATION_GET(INTERFACE_DEPENDENT%FIELD_VARIABLES(1)%PTR%FIELD,FIELD_U_VARIABLE_TYPE, &
                  !  & component_idx,interpolation_type,ERR,ERROR,*999)
                  !CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD, &
                  !  & FIELD_U_VARIABLE_TYPE,component_idx,interpolation_type,ERR,ERROR,*999)
                  !CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD, &
                  !  & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,interpolation_type,ERR,ERROR,*999)
                ENDDO !component_idx
                CALL FIELD_SCALING_TYPE_GET(INTERFACE_CONDITION%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                !Check the Lagrange field
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              ENDIF
              !Set pointers
              IF(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD_AUTO_CREATED) THEN
                LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
              ELSE
                INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD=>LAGRANGE_FIELD
              ENDIF
            ELSE
              CALL FLAG_ERROR("The interface parent region is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The interface interface conditions is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition Lagrange information and deallocate all memory.
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FINALISE(INTERFACE_LAGRANGE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: INTERFACE_LAGRANGE !<A pointer to the interface condition Lagrange information to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_LAGRANGE)) THEN
      DEALLOCATE(INTERFACE_LAGRANGE)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition Lagrange information.
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<The pointer to the interface condition to initialise to initialise the Lagrange information for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
        CALL FLAG_ERROR("Interface condition Lagrange is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE_CONDITION%LAGRANGE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface condition Lagrange.",ERR,ERROR,*999)
        INTERFACE_CONDITION%LAGRANGE%INTERFACE_CONDITION=>INTERFACE_CONDITION
        INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FINISHED=.FALSE.
        INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD)
        INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITION_LAGRANGE_FINALISE(INTERFACE_CONDITION%LAGRANGE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an interface condition's penalty field'. \see OPENCMISS::CMISSInterfaceConditionPenaltyConditionCreateFinish
  SUBROUTINE INTERFACE_CONDITION_PENALTY_FIELD_CREATE_FINISH(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finish creating the penalty field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%PENALTY)) THEN
        IF(INTERFACE_CONDITION%PENALTY%PENALTY_FINISHED) THEN
          CALL FLAG_ERROR("Interface condition penalty field has already been finished.",ERR,ERROR,*999)
        ELSE
          !Finish the penalty field creation
          IF(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD_AUTO_CREATED) THEN
            CALL FIELD_CREATE_FINISH(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,ERR,ERROR,*999)
          ENDIF
          INTERFACE_CONDITION%PENALTY%PENALTY_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition penalty is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_CONDITION_PENALTY_FIELD_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the process of creating the penalty field for interface condition. \see OPENCMISS::CMISSInterfaceConditionPenaltyFieldCreateStart
  SUBROUTINE INTERFACE_CONDITION_PENALTY_FIELD_CREATE_START(INTERFACE_CONDITION,PENALTY_FIELD_USER_NUMBER,PENALTY_FIELD, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to create the penalty field on
    INTEGER(INTG), INTENT(IN) :: PENALTY_FIELD_USER_NUMBER !<The user specified penalty field number
    TYPE(FIELD_TYPE), POINTER :: PENALTY_FIELD !<If associated on entry, a pointer to the user created penalty field which has the same user number as the specified penalty field user number. If not associated on entry, on exit, a pointer to the created penalty field for the interface condition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_SCALING_TYPE,INTERPOLATION_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(REGION_TYPE), POINTER :: INTERFACE_REGION,PENALTY_FIELD_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%PENALTY)) THEN
        CALL FLAG_ERROR("Interface condition penalty is already associated.",ERR,ERROR,*999)
      ELSE
        INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
        IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
          INTERFACE=>INTERFACE_CONDITION%INTERFACE
          IF(ASSOCIATED(INTERFACE)) THEN
            INTERFACE_REGION=>INTERFACE%PARENT_REGION
            IF(ASSOCIATED(INTERFACE_REGION)) THEN
              IF(ASSOCIATED(PENALTY_FIELD)) THEN
                !Check the penalty field has been finished
                IF(PENALTY_FIELD%FIELD_FINISHED) THEN
                  !Check the user numbers match
                  IF(PENALTY_FIELD_USER_NUMBER/=PENALTY_FIELD%USER_NUMBER) THEN
                    LOCAL_ERROR="The specified penalty field user number of "// &
                      & TRIM(NUMBER_TO_VSTRING(PENALTY_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                      & " does not match the user number of the specified penalty field of "// &
                      & TRIM(NUMBER_TO_VSTRING(PENALTY_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  PENALTY_FIELD_REGION=>PENALTY_FIELD%REGION
                  IF(ASSOCIATED(PENALTY_FIELD_REGION)) THEN
                    !Check the field is defined on the same region as the interface
                    IF(PENALTY_FIELD_REGION%USER_NUMBER/=INTERFACE_REGION%USER_NUMBER) THEN
                      LOCAL_ERROR="Invalid region setup. The specified penalty field has been created on interface number "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" in parent region number "// &
                        & TRIM(NUMBER_TO_VSTRING(PENALTY_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                        & " and the specified interface has been created in parent region number "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The penalty field region is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The specified penalty field has not been finished.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !Check the user number has not already been used for a field in this region.
                NULLIFY(FIELD)
                CALL FIELD_USER_NUMBER_FIND(PENALTY_FIELD_USER_NUMBER,INTERFACE,FIELD,ERR,ERROR,*999)
                IF(ASSOCIATED(FIELD)) THEN
                  LOCAL_ERROR="The specified penalty field user number of "// &
                    & TRIM(NUMBER_TO_VSTRING(PENALTY_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                    & " has already been used to create a field on interface number "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
              CALL INTERFACE_CONDITION_PENALTY_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
              IF(.NOT.ASSOCIATED(PENALTY_FIELD)) THEN
                !Create the penalty field
                INTERFACE_CONDITION%PENALTY%PENALTY_FIELD_AUTO_CREATED=.TRUE.
                CALL FIELD_CREATE_START(PENALTY_FIELD_USER_NUMBER,INTERFACE_CONDITION%INTERFACE,INTERFACE_CONDITION%PENALTY% &
                  & PENALTY_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,"Penalty Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,FIELD_DEPENDENT_TYPE, &
                  & ERR,ERROR,*999)
                NULLIFY(GEOMETRIC_DECOMPOSITION)
                CALL FIELD_MESH_DECOMPOSITION_GET(INTERFACE_CONDITION%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,INTERFACE_CONDITION%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,FIELD_U_VARIABLE_TYPE,"Alpha", &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                !Default the number of component to the first variable of the interface dependent field's number of components
                !CALL FIELD_NUMBER_OF_COMPONENTS_SET(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,FIELD_U_VARIABLE_TYPE, &
                !  & INTERFACE_DEPENDENT%FIELD_VARIABLES(1)%PTR%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                DO component_idx=1,INTERFACE_CONDITION%LAGRANGE%NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_GET(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,INTERPOLATION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,INTERPOLATION_TYPE,ERR,ERROR,*999)
                  !CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD, &
                  !  & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                CALL FIELD_SCALING_TYPE_GET(INTERFACE_CONDITION%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                !Check the penalty field
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              ENDIF
              !Set pointers
              IF(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD_AUTO_CREATED) THEN
                PENALTY_FIELD=>INTERFACE_CONDITION%PENALTY%PENALTY_FIELD
              ELSE
                INTERFACE_CONDITION%PENALTY%PENALTY_FIELD=>PENALTY_FIELD
              ENDIF
            ELSE
              CALL FLAG_ERROR("The interface parent region is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The interface interface conditions is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface condition dependent is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_PENALTY_FIELD_CREATE_START")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_PENALTY_FIELD_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition penalty information and deallocate all memory.
  SUBROUTINE INTERFACE_CONDITION_PENALTY_FINALISE(INTERFACE_PENALTY,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_PENALTY_TYPE), POINTER :: INTERFACE_PENALTY !<A pointer to the interface condition penalty information to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_PENALTY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_PENALTY)) THEN
      DEALLOCATE(INTERFACE_PENALTY)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_PENALTY_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_PENALTY_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_PENALTY_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_PENALTY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition penalty information.
  SUBROUTINE INTERFACE_CONDITION_PENALTY_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<The pointer to the interface condition to initialise to initialise the penalty information for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_PENALTY_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%PENALTY)) THEN
        CALL FLAG_ERROR("Interface condition penalty is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE_CONDITION%PENALTY,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface condition penalty.",ERR,ERROR,*999)
        INTERFACE_CONDITION%PENALTY%INTERFACE_CONDITION=>INTERFACE_CONDITION
        INTERFACE_CONDITION%PENALTY%PENALTY_FINISHED=.FALSE.
        INTERFACE_CONDITION%PENALTY%PENALTY_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(INTERFACE_CONDITION%PENALTY%PENALTY_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_PENALTY_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITION_PENALTY_FINALISE(INTERFACE_CONDITION%PENALTY,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITION_PENALTY_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_PENALTY_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_PENALTY_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the interface condition method \see OPENCMISS::CMISSInterfaceConditionMethodGet
  SUBROUTINE INTERFACE_CONDITION_METHOD_GET(INTERFACE_CONDITION,INTERFACE_CONDITION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to get the method for
    INTEGER(INTG), INTENT(OUT) :: INTERFACE_CONDITION_METHOD !<On return, the interface condition method. \see INTERFACE_CONDITIONS_Methods,INTERFACE_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_METHOD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        INTERFACE_CONDITION_METHOD=INTERFACE_CONDITION%METHOD
      ELSE
        CALL FLAG_ERROR("Interface condition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_METHOD_GET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_METHOD_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_METHOD_GET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_METHOD_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the interface condition method \see OPENCMISS::CMISSInterfaceConditionMethodSet
  SUBROUTINE INTERFACE_CONDITION_METHOD_SET(INTERFACE_CONDITION,INTERFACE_CONDITION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to set the method for
    INTEGER(INTG), INTENT(IN) :: INTERFACE_CONDITION_METHOD !<The interface condition method to set. \see INTERFACE_CONDITIONS_Methods,INTERFACE_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Interface condition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(INTERFACE_CONDITION_METHOD)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_POINT_TO_POINT_METHOD
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD
         CASE(INTERFACE_CONDITION_PENALTY_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_PENALTY_METHOD
       CASE DEFAULT
          LOCAL_ERROR="The specified interface condition method of "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION_METHOD,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_METHOD_SET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_METHOD_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_METHOD_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_METHOD_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the interface condition operator \see OPENCMISS::CMISSInterfaceConditionOperatorGet
  SUBROUTINE INTERFACE_CONDITION_OPERATOR_GET(INTERFACE_CONDITION,INTERFACE_CONDITION_OPERATOR,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to get the operator for
    INTEGER(INTG), INTENT(OUT) :: INTERFACE_CONDITION_OPERATOR !<On return, the interface condition operator. \see INTERFACE_CONDITIONS_Operators,INTERFACE_CONDITIONS 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_OPERATOR_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        INTERFACE_CONDITION_OPERATOR=INTERFACE_CONDITION%OPERATOR
      ELSE
        CALL FLAG_ERROR("Interface condition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_GET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_OPERATOR_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_GET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_OPERATOR_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the interface condition operator \see OPENCMISS::CMISSInterfaceConditionOperatorSet
  SUBROUTINE INTERFACE_CONDITION_OPERATOR_SET(INTERFACE_CONDITION,INTERFACE_CONDITION_OPERATOR,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to set the operator for
    INTEGER(INTG), INTENT(IN) :: INTERFACE_CONDITION_OPERATOR !<The interface condition operator to set. \see INTERFACE_CONDITIONS_Operators,INTERFACE_CONDITIONS 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_OPERATOR_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Interface condition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(INTERFACE_CONDITION_OPERATOR)
        CASE(INTERFACE_CONDITION_FIELD_GAUSS_CONTINUITY_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_FIELD_GAUSS_CONTINUITY_OPERATOR
        CASE(INTERFACE_CONDITION_FIELD_NODE_CONTINUITY_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_FIELD_NODE_CONTINUITY_OPERATOR
        CASE(INTERFACE_CONDITION_FIELD_GAUSS_NODE_CONTINUITY_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_FIELD_GAUSS_NODE_CONTINUITY_OPERATOR
        CASE(INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR
        CASE(INTERFACE_CONDITION_SOLID_FLUID_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_SOLID_FLUID_OPERATOR
        CASE(INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR
        CASE(INTERFACE_CONDITION_FRICTIONLESS_CONTACT_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_FRICTIONLESS_CONTACT_OPERATOR
        CASE DEFAULT
          LOCAL_ERROR="The specified interface condition operator of "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION_OPERATOR,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_SET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_OPERATOR_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_OPERATOR_SET
  
  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an interface condition.
  SUBROUTINE INTERFACE_CONDITION_RESIDUAL_EVALUATE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            IF(INTERFACE_CONDITION%OPERATOR==INTERFACE_CONDITION_FRICTIONLESS_CONTACT_OPERATOR) THEN
              IF (INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD%VARIABLES(1)%COMPONENTS(1)%INTERPOLATION_TYPE== &
                  & FIELD_DATA_POINT_BASED_INTERPOLATION) THEN
                CALL INTERFACE_CONDITION_DATA_REPROJECTION(INTERFACE_CONDITION,ERR,ERROR,*999)
              ENDIF
            ENDIF
            CALL INTERFACE_CONDITION_RESIDUAL_EVALUATE_FEM(INTERFACE_CONDITION,ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Interface equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_RESIDUAL_EVALUATE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_CONDITION_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an interface condition using the finite element method
  SUBROUTINE INTERFACE_CONDITION_RESIDUAL_EVALUATE_FEM(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: LAGRANGE
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
 
    CALL ENTERS("INTERFACE_CONDITION_RESIDUAL_EVALUATE_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
      IF(ASSOCIATED(LAGRANGE)) THEN
        LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
        IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
            IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
!!Do we need to transfer parameter sets???
              !Initialise the matrices and rhs vector
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_VALUES_INITIALISE()")
#endif
              CALL INTERFACE_MATRICES_VALUES_INITIALISE(INTERFACE_MATRICES,0.0_DP,ERR,ERROR,*999)
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_VALUES_INITIALISE()")
#endif
              !Assemble the elements
              !Allocate the element matrices 
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_INITIALISE()")
#endif
              CALL INTERFACE_MATRICES_ELEMENT_INITIALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
              ELEMENTS_MAPPING=>LAGRANGE_FIELD%DECOMPOSITION%DOMAIN(LAGRANGE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
                & MAPPINGS%ELEMENTS
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_INITIALISE()")
#endif
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for interface setup and initialisation = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for interface setup and initialisation = ", &
                  & SYSTEM_ELAPSED,ERR,ERROR,*999)
                ELEMENT_USER_ELAPSED=0.0_SP
                ELEMENT_SYSTEM_ELAPSED=0.0_SP
              ENDIF
              NUMBER_OF_TIMES=0
              !Loop over the internal elements
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("Internal Elements Loop")
#endif
              DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
                ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
                NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
                CALL INTERFACE_MATRICES_ELEMENT_CALCULATE(INTERFACE_MATRICES,ne,ERR,ERROR,*999)
                CALL INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ne,ERR,ERROR,*999)
                CALL INTERFACE_MATRICES_ELEMENT_ADD(INTERFACE_MATRICES,ERR,ERROR,*999)
              ENDDO !element_idx                  
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("Internal Elements Loop")
#endif
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
                SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
                ELEMENT_USER_ELAPSED=USER_ELAPSED
                ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal interface assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal interface assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
              ENDIF
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
                SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)              
              ENDIF
              !Loop over the boundary and ghost elements
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("Boundary and Ghost Elements Loop")
#endif
              DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
                ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
                NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
                CALL INTERFACE_MATRICES_ELEMENT_CALCULATE(INTERFACE_MATRICES,ne,ERR,ERROR,*999)
                CALL INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ne,ERR,ERROR,*999)
                CALL INTERFACE_MATRICES_ELEMENT_ADD(INTERFACE_MATRICES,ERR,ERROR,*999)
              ENDDO !element_idx
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("Boundary and Ghost Elements Loop")
#endif
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
                SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
                ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
                ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost interface assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE, &
                  & "System time for boundary+ghINTERFACE_CONDITION_PENALTY_METHODost interface assembly = ", &
                  & SYSTEM_ELAPSED,ERR,ERROR,*999)
                IF(NUMBER_OF_TIMES>0) THEN
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for interface assembly = ", &
                    & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for interface assembly = ", &
                    & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                ENDIF
              ENDIF
              !Finalise the element matrices
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_FINALISE()")
#endif
              CALL INTERFACE_MATRICES_ELEMENT_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_FINALISE()")
#endif
              !Output equations matrices and RHS vector if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_MATRIX_OUTPUT) THEN
                CALL INTERFACE_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,INTERFACE_MATRICES,ERR,ERROR,*999)
              ENDIF
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for interface equations assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for interface equations assembly = ", &
                  & SYSTEM_ELAPSED,ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface condition Lagrange field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition Lagrange is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_RESIDUAL_EVALUATE_FEM")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_RESIDUAL_EVALUATE_FEM",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_RESIDUAL_EVALUATE_FEM")
    RETURN 1
    
  END SUBROUTINE INTERFACE_CONDITION_RESIDUAL_EVALUATE_FEM
  
  !
  !================================================================================================================================
  !

   !>Calculates the element stiffness matries for the given element number for a finite element interface equations.
  SUBROUTINE INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    INTEGER(INTG) :: ng, mh, mhs, ms, nh, nhs, ns, mhc, dir,derivativeIdx,derivative,localNode,localNodeIdx,localElementNode
    INTEGER(INTG) :: connectedLine,lineNodeIdx,decompositionLineNumber,localLineNode,localLineNodeIdx
    INTEGER(INTG) :: connectedFace,faceNodeIdx,decompositionFaceNumber,localFaceNode,localFaceNodeIdx
    INTEGER(INTG) :: interfaceNode,interfaceDerivative,coupledMeshElementNumber,elementNumberGap
    REAL(DP) :: XI(3),RWG,PGMSI,PGNSI,MATRIX_COEFFICIENT,ROWBASIS,XIGap(3),gap(3),penetration,junk
    INTEGER(INTG) :: interface_matrix_idx,interface_matrix_gap_idx
    INTEGER(INTG) :: INTERFACE_MATRIX_VARIABLE_TYPE,LAGRANGE_VARIABLE_TYPE
    TYPE(BASIS_TYPE), POINTER :: INTERFACE_DEPENDENT_BASIS,COUPLED_MESH_BASIS,INTERFACE_GEOMETRIC_BASIS, &
      & INTERFACE_PENALTY_BASIS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: INTERFACE_MATRIX_VARIABLE,LAGRANGE_VARIABLE
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX,ELEMENT_MATRIX,INTERFACE_ELEMENT_MATRIX
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE), POINTER :: INTERFACE_INTERPOLATION
    TYPE(FIELD_TYPE), POINTER :: INTERFACE_MATRIX_DEPENDENT_FIELD,INTERFACE_DEPENDENT_FIELD,INTERFACE_GEOMETRIC_FIELD, &
      & INTERFACE_PENALTY_FIELD,COUPLED_MESH_GEOMETRIC_FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: INTERFACE_QUADRATURE_SCHEME
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    TYPE(DOMAIN_LINE_TYPE), POINTER :: COUPLED_MESH_DOMAIN_LINE
    TYPE(DOMAIN_FACE_TYPE), POINTER :: COUPLED_MESH_DOMAIN_FACE
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VEC
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: INTERPOLATED_POINT_METRICS !<A pointer to the interpolated point metric information to calculate the position etc. for
    TYPE(INTERFACE_POINTS_CONNECTIVITY_TYPE), POINTER :: POINTS_CONNECTIVITY !<A pointer to the points connectivity
    TYPE(DOMAIN_DATA_POINTS_TYPE), POINTER :: DATA_POINTS_TOPOLOGY
    LOGICAL :: REVERSE_NORMAL !<logical to decide if the normal needs to be reversed to get the outward normal if .TRUE. then reverse 
    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    

    INTEGER(INTG) :: coupledMeshNumberOfXi,ELEM
    REAL(DP) :: POSITION(3),NORMAL(3),TANGENTS(3,3),XI_POSITION(2)
    REAL(DP) :: dofValue, weight
    INTEGER(INTG) :: rowGlobalDofNumber
    INTEGER(INTG) :: localContactNumber,localContactXiNormal

    
    INTEGER(INTG) :: dataPointIdx,rowComponentIdx,rowParameterIdx,rowIdx,colComponentIdx,colIdx, &
      & colParameterIdx,coupledMeshElementIdx,dataPointNumber,rowMeshComponentNumber, &
      & globalContactNumber,coupledMeshNumberOfGeometricComponents,InterfaceNumberOfGeometricComponent,gaussPointIdx, &
      & numberOfCoupledMesh
    
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        SELECT CASE(INTERFACE_CONDITION%METHOD)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
          !Pointers to interface variables (columns of interface element matrix)
          INTERFACE_INTERPOLATION=>INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION
          INTERFACE_GEOMETRIC_FIELD=>INTERFACE_INTERPOLATION%GEOMETRIC_FIELD
          INTERFACE_GEOMETRIC_BASIS=>INTERFACE_GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(INTERFACE_GEOMETRIC_FIELD% &
            & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          INTERFACE_DEPENDENT_FIELD=>INTERFACE_INTERPOLATION%DEPENDENT_FIELD
          INTERFACE_DEPENDENT_BASIS=>INTERFACE_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(INTERFACE_DEPENDENT_FIELD% &
            & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            INTERFACE_PENALTY_FIELD=>INTERFACE_INTERPOLATION%PENALTY_FIELD
            SELECT CASE(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD%VARIABLES(1)%COMPONENTS(1)%INTERPOLATION_TYPE)
            CASE(FIELD_NODE_BASED_INTERPOLATION) !Mesh connectivity
              INTERFACE_PENALTY_BASIS=>INTERFACE_PENALTY_FIELD%DECOMPOSITION%DOMAIN(INTERFACE_PENALTY_FIELD% &
                & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,INTERFACE_INTERPOLATION% &
                & PENALTY_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            END SELECT
          END SELECT
          !Integrate using the interface quadrature scheme
          INTERFACE_QUADRATURE_SCHEME=>INTERFACE_GEOMETRIC_BASIS%QUADRATURE% &
            & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          LAGRANGE_VARIABLE=>INTERFACE_EQUATIONS%INTERFACE_MAPPING%LAGRANGE_VARIABLE
          LAGRANGE_VARIABLE_TYPE=LAGRANGE_VARIABLE%VARIABLE_TYPE
          !Get element interpolation parameters from the first geometric interpolation set (to get Jacobian for interface surface integral)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,INTERFACE_INTERPOLATION% &
            & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
          MATRIX_COEFFICIENT=1.0_DP
          DO interface_matrix_idx=1,INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
            IF(INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%UPDATE_MATRIX) THEN
              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
              IF(interface_matrix_idx>1) THEN
                MATRIX_COEFFICIENT=-1.0_DP
              ENDIF
              !Pointers to the interface_matrix_idx'th coupled mesh variables (rows of interface element matrix)
              INTERFACE_MATRIX_VARIABLE=>INTERFACE_EQUATIONS%INTERFACE_MAPPING% & 
                & INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE
              INTERFACE_MATRIX_VARIABLE_TYPE=INTERFACE_MATRIX_VARIABLE%VARIABLE_TYPE
              INTERFACE_ELEMENT_MATRIX=>INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%ELEMENT_MATRIX
              !Hard code-interpolation type for the first component of the first variable   
              SELECT CASE(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD%VARIABLES(1)%COMPONENTS(1)%INTERPOLATION_TYPE)
              CASE(FIELD_NODE_BASED_INTERPOLATION) !Mesh connectivity
                ELEMENT_CONNECTIVITY=>INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY% &
                  & ELEMENT_CONNECTIVITY(ELEMENT_NUMBER,interface_matrix_idx)
                coupledMeshElementNumber=ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER
              CASE(FIELD_DATA_POINT_BASED_INTERPOLATION) !Points connecitivity
                INTERFACE_GEOMETRIC_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%GEOMETRIC_FIELD      
                POINTS_CONNECTIVITY=>INTERFACE_CONDITION%INTERFACE%POINTS_CONNECTIVITY
                DATA_POINTS_TOPOLOGY=>POINTS_CONNECTIVITY%INTERFACE_MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR% &
                  & DOMAIN(INTERFACE_GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%DATA_POINTS 
              END SELECT
              SELECT CASE(INTERFACE_CONDITION%OPERATOR)
              CASE(INTERFACE_CONDITION_FIELD_GAUSS_CONTINUITY_OPERATOR)
                !Loop over gauss points
                DO ng=1,INTERFACE_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
                  CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,INTERFACE_INTERPOLATION% &
                    & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(INTERFACE_GEOMETRIC_BASIS%NUMBER_OF_XI,INTERFACE_INTERPOLATION% &
                    & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  RWG=INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                    & JACOBIAN*INTERFACE_QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
                  IF(INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                      & interface_matrix_idx==INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
                    CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,INTERFACE_INTERPOLATION% &
                      & PENALTY_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                    mhs=0
                    DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                      !Loop over the Lagrange variable matrix rows
                      DO ms=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        mhs=mhs+1
                        INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs)- &
                          & (1.0_DP/INTERFACE_INTERPOLATION%PENALTY_INTERPOLATION(1)% &
                          & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1))*PGNSI**2.0_DP*RWG
                      ENDDO !ms
                    ENDDO !mh
                  ELSE
                    !\todo defaults to first mesh component, Generalise
                    XI(1:INTERFACE_DEPENDENT_BASIS%NUMBER_OF_XI)=INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM( &
                      & INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY,ELEMENT_NUMBER,interface_matrix_idx,ng,ERR,ERROR)
                    ! Loop over number of Lagrange variable components as not all components in the dependent field variable may be coupled
                    !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                    DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                      mhc=INTERFACE_MATRIX_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                      COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% & 
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                      SELECT CASE(INTERFACE_DEPENDENT_BASIS%NUMBER_OF_XI)
                      CASE(1) !1D interface (line)
                        connectedLine = ELEMENT_CONNECTIVITY%CONNECTED_LINE
                        decompositionLineNumber=INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY% &
                          & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(connectedLine)
                        COUPLED_MESH_DOMAIN_LINE=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% &
                          & LINES%LINES(decompositionLineNumber)
                        DO localLineNodeIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(connectedLine)
                          localElementNode=COUPLED_MESH_BASIS%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                          DO derivativeIdx=1,COUPLED_MESH_DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES!COUPLED_MESH_BASIS%NUMBER_OF_DERIVATIVES(localNode)
                            derivative=COUPLED_MESH_BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                            derivative=COUPLED_MESH_DOMAIN_LINE%DERIVATIVES_IN_LINE(1,derivativeIdx,localLineNodeIdx)
                            ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                            IF (mh==4) THEN
                              PGMSI=1.0_DP
                            ELSE
                              PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,ms,NO_PART_DERIV,XI,ERR,ERROR)
                            ENDIF
                            mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                            DO interfaceNode=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_NODES
                              DO interfaceDerivative=1,INTERFACE_DEPENDENT_BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                                ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                                PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                                nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                                !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                                INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                                  & PGNSI*PGMSI*RWG*MATRIX_COEFFICIENT
                              ENDDO !interfaceDerivative
                            ENDDO !interfaceNode
                          ENDDO !derivativeIdx
                        ENDDO !lineNodeIdx
                      CASE(2) !2D interface (face)
                        SELECT CASE(COUPLED_MESH_BASIS%NUMBER_OF_XI)
                        CASE(2) !Coupled Mesh has 2 xi directions
                          DO localElementNode=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES
                            DO derivative=1,COUPLED_MESH_DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES!COUPLED_MESH_BASIS%NUMBER_OF_DERIVATIVES(localNode)
                              ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                              IF (mh==4) THEN
                                PGMSI=1.0_DP
                              ELSE
                                PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,ms,NO_PART_DERIV, &
                                  & XI(1:COUPLED_MESH_BASIS%NUMBER_OF_XI),ERR,ERROR)
                              ENDIF
                              mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                              DO interfaceNode=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_NODES
                                DO interfaceDerivative=1,INTERFACE_DEPENDENT_BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                  !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                                  ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                                  PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                                  nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                                  !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                                  INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                                    & PGNSI*PGMSI*RWG*MATRIX_COEFFICIENT
                                ENDDO !interfaceDerivative
                              ENDDO !interfaceNode
                            ENDDO !derivative
                          ENDDO !localElementNode
                        CASE(3) !Coupled Mesh has 3 xi directions
                          connectedFace = ELEMENT_CONNECTIVITY%CONNECTED_FACE
                          decompositionFaceNumber=INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY% &
                            & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(connectedFace)
                          COUPLED_MESH_DOMAIN_FACE=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% &
                            & FACES%FACES(decompositionFaceNumber)
                          DO localFaceNodeIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                            localElementNode=COUPLED_MESH_BASIS%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                            DO derivativeIdx=1,COUPLED_MESH_DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                              derivative=COUPLED_MESH_BASIS% &
                                & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(derivativeIdx,localFaceNodeIdx,connectedFace)
                              ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                              IF (mh==4) THEN
                                PGMSI=1.0_DP
                              ELSE
                                PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,ms,NO_PART_DERIV, &
                                  & XI(1:COUPLED_MESH_BASIS%NUMBER_OF_XI),ERR,ERROR)
                              ENDIF
                              mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                              DO interfaceNode=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_NODES
                                DO interfaceDerivative=1,INTERFACE_DEPENDENT_BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                  !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                                  ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                                  PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                                  nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                                  !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                                  INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                                    & PGNSI*PGMSI*RWG*MATRIX_COEFFICIENT
                                ENDDO !interfaceDerivative
                              ENDDO !interfaceNode
                            ENDDO !derivativeIdx
                          ENDDO !FaceNodeIdx
                        END SELECT
                      END SELECT
                    ENDDO !mh
                  ENDIF
                ENDDO !ng
                !Scale factor adjustment
                !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix contribution to the residual for non-linear problems
                !\todo update looping of variables/components for non-zero matrix elements as done above 
                IF(INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                  & interface_matrix_idx==INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
                  !Scale factor adjustment for the Lagrange Variable (columns)
                  IF(INTERFACE_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER, &
                      & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR, &
                      & ERR,ERROR,*999)
                    mhs=0
                    !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                    !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                    DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                      !Loop over element Lagrange variable rows
                      DO ms=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        mhs=mhs+1
                        INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs) * &
                          & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
                          & INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR%SCALE_FACTORS(ms,mh)**2
                      ENDDO !ms
                    ENDDO !mh
                  ENDIF
                ELSE
                  !Scale factor adjustment for the Lagrange Variable (columns)
                  IF(INTERFACE_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER, &
                      & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR, &
                      & ERR,ERROR,*999)
                    mhs=0
                    !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                    !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                    DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                      mhc=INTERFACE_MATRIX_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                      COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% & 
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                      !Loop over element rows
                      DO ms=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        mhs=mhs+1
                        nhs=0
                        !Loop over element columns
                        DO nh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                          DO ns=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            nhs=nhs+1
                            INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs) * &
                            & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                          ENDDO !ns
                        ENDDO !nh
                      ENDDO !ms
                    ENDDO !mh
                  ENDIF
                  !Scale factor adjustment for the row dependent variable
                  IF(INTERFACE_MATRIX_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(coupledMeshElementNumber, &
                      & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                      & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                    mhs=0
                    DO mh=1,INTERFACE_MATRIX_VARIABLE%NUMBER_OF_COMPONENTS
                      !Loop over element rows
                      DO ms=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        mhs=mhs+1
                        nhs=0
                        !Loop over element columns
                        DO nh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                          DO ns=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            nhs=nhs+1
                            INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                            & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                            & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR% &
                            & SCALE_FACTORS(ms,mh)
                          ENDDO !ns
                        ENDDO !nh
                      ENDDO !ms
                    ENDDO !mh
                  ENDIF
                ENDIF
              CASE(INTERFACE_CONDITION_FIELD_NODE_CONTINUITY_OPERATOR)
                ! Loop over number of Lagrange variable components as not all components in the dependent field variable may be coupled
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                  mhc=INTERFACE_MATRIX_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% & 
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                  SELECT CASE(INTERFACE_DEPENDENT_BASIS%NUMBER_OF_XI)
                  CASE(1) !1D interface (line)
                    connectedLine = ELEMENT_CONNECTIVITY%CONNECTED_LINE
                    decompositionLineNumber=INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(connectedLine)
                    COUPLED_MESH_DOMAIN_LINE=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% &
                      & LINES%LINES(decompositionLineNumber)
                    DO localLineNodeIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(connectedLine)
                      localElementNode=COUPLED_MESH_BASIS%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                      IF (mh==4) THEN
                        PGMSI=1.0_DP
                      ELSE
                        !Calculate PGMSI based on node NO_PART_DERIV
                        !\todo defaults to first mesh component, Generalise
                        PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,COUPLED_MESH_BASIS% &
                          & ELEMENT_PARAMETER_INDEX(1,localElementNode),NO_PART_DERIV,INTERFACE_CONDITION%INTERFACE% &
                          & MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ELEMENT_NUMBER,interface_matrix_idx)% &
                          & XI(:,1,localLineNodeIdx),ERR,ERROR)
                      ENDIF
                      IF (PGMSI<1.0_DP .AND. PGMSI >ZERO_TOLERANCE)THEN
                        PGMSI=PGMSI*2.0_DP
                      ENDIF
                      DO derivativeIdx=1,COUPLED_MESH_DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES!COUPLED_MESH_BASIS%NUMBER_OF_DERIVATIVES(localNode)
                        derivative=COUPLED_MESH_BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                        derivative=COUPLED_MESH_DOMAIN_LINE%DERIVATIVES_IN_LINE(1,derivativeIdx,localLineNodeIdx)
                        ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                        mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                        !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                        ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(derivativeIdx,localLineNodeIdx)
                        nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                        !Direct map between nodal lagrange parameter and corresponding coupled mesh parameter therefore PGNSI=1, so not required below
                        !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                        INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & PGMSI*MATRIX_COEFFICIENT
                      ENDDO !derivativeIdx
                    ENDDO !lineNodeIdx
                  CASE(2) !2D interface (face)
                    SELECT CASE(COUPLED_MESH_BASIS%NUMBER_OF_XI)
                    CASE(2) !Coupled Mesh has 2 xi directions
                      DO localElementNode=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES
                        DO derivative=1,COUPLED_MESH_DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES!COUPLED_MESH_BASIS%NUMBER_OF_DERIVATIVES(localNode)
                          ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                          !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                          ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                          !Direct map between nodal lagrange parameter and corresponding coupled mesh parameter therefore PGNSI=1, so not required below
                          !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                          INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & MATRIX_COEFFICIENT
                        ENDDO !derivative
                      ENDDO !localElementNode
                    CASE(3) !Coupled Mesh has 3 xi directions
                      connectedFace = ELEMENT_CONNECTIVITY%CONNECTED_FACE
                      decompositionFaceNumber=INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(connectedFace)
                      COUPLED_MESH_DOMAIN_FACE=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% &
                        & FACES%FACES(decompositionFaceNumber)
                      DO localFaceNodeIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                        localElementNode=COUPLED_MESH_BASIS%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                        DO derivativeIdx=1,COUPLED_MESH_DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                          derivative=COUPLED_MESH_BASIS% &
                            & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(derivativeIdx,localFaceNodeIdx,connectedFace)
                          ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                          !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                          ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(derivativeIdx,localFaceNodeIdx)
                          nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                          !Direct map between nodal lagrange parameter and corresponding coupled mesh parameter therefore PGNSI=1, so not required below
                          !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                          INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & MATRIX_COEFFICIENT
                        ENDDO !derivativeIdx
                      ENDDO !FaceNodeIdx
                    END SELECT
                  END SELECT
                ENDDO !mh
                !Scale factor adjustment
                !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix contribution to the residual for non-linear problems
                !\todo update looping of variables/components for non-zero matrix elements as done above 
                IF(INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                  & interface_matrix_idx==INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
                  !Scale factor adjustment for the Lagrange Variable (columns)
                  IF(INTERFACE_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER, &
                      & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR, &
                      & ERR,ERROR,*999)
                    mhs=0
                    !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                    !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                    DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                      !Loop over element Lagrange variable rows
                      DO ms=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        mhs=mhs+1
                        INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs) * &
                          & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
                          & INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR%SCALE_FACTORS(ms,mh)**2
                      ENDDO !ms
                    ENDDO !mh
                  ENDIF
                ELSE
                  !Scale factor adjustment for the Lagrange Variable (columns)
                  IF(INTERFACE_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER, &
                      & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR, &
                      & ERR,ERROR,*999)
                    mhs=0
                    !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                    !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                    DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                      mhc=INTERFACE_MATRIX_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                      COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% & 
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                      !Loop over element rows
                      DO ms=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        mhs=mhs+1
                        nhs=0
                        !Loop over element columns
                        DO nh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                          DO ns=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            nhs=nhs+1
                            INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs) * &
                            & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                          ENDDO !ns
                        ENDDO !nh
                      ENDDO !ms
                    ENDDO !mh
                  ENDIF
                  !Scale factor adjustment for the row dependent variable
                  IF(INTERFACE_MATRIX_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(coupledMeshElementNumber, &
                      & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                      & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                    mhs=0
                    DO mh=1,INTERFACE_MATRIX_VARIABLE%NUMBER_OF_COMPONENTS
                      !Loop over element rows
                      DO ms=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        mhs=mhs+1
                        nhs=0
                        !Loop over element columns
                        DO nh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                          DO ns=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            nhs=nhs+1
                            INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                            & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                            & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR% &
                            & SCALE_FACTORS(ms,mh)
                          ENDDO !ns
                        ENDDO !nh
                      ENDDO !ms
                    ENDDO !mh
                  ENDIF
                ENDIF
                        
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
              !Contact Mechanics starts here.                          
              CASE(INTERFACE_CONDITION_FRICTIONLESS_CONTACT_OPERATOR)
                INTERFACE_ELEMENT_MATRIX=>INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%ELEMENT_MATRIX                    
  
                INTERFACE_MATRIX_VARIABLE=>INTERFACE_EQUATIONS%INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS &
                  & (interface_matrix_idx)%VARIABLE
                INTERFACE_MATRIX_VARIABLE_TYPE=INTERFACE_MATRIX_VARIABLE%VARIABLE_TYPE
                InterfaceNumberOfGeometricComponent=INTERFACE_GEOMETRIC_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS 
                IF (.NOT. (INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                    & interface_matrix_idx==INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES)) THEN
                  COUPLED_MESH_GEOMETRIC_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                    & GEOMETRIC_FIELD 
                  coupledMeshNumberOfGeometricComponents=COUPLED_MESH_GEOMETRIC_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS 
                  INTERFACE_MATRIX_DEPENDENT_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION% &
                    & VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_FIELD   
                ENDIF    
                SELECT CASE(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD%VARIABLES(1)%COMPONENTS(1)%INTERPOLATION_TYPE)
                CASE(FIELD_NODE_BASED_INTERPOLATION) !Mesh connectivity
                  !********************************** CODE GOES HERE *********************************************************
                  !Local number of the line in contact
                  SELECT CASE(INTERFACE_DEPENDENT_BASIS%NUMBER_OF_XI)
                  CASE(1) !1D interface (line)
                    localContactNumber=ELEMENT_CONNECTIVITY%CONNECTED_LINE
                    !Global number of the line in contact
                    globalContactNumber=INTERFACE_CONDITION%INTERFACE%COUPLED_MESHES(interface_matrix_idx)%PTR%DECOMPOSITIONS% &
                      & DECOMPOSITIONS(1)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(localContactNumber)   
                  CASE(2) !2D interface (face)
                    localContactNumber=ELEMENT_CONNECTIVITY%CONNECTED_FACE
                    globalContactNumber=INTERFACE_CONDITION%INTERFACE%COUPLED_MESHES(interface_matrix_idx)%PTR%DECOMPOSITIONS% &
                      & DECOMPOSITIONS(1)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(localContactNumber)
                  END SELECT   
                  
                  !The xi direction of the surface/line normal
                  localContactXiNormal=ELEMENT_CONNECTIVITY%COUPLED_MESH_CONTACT_XI_NORMAL 
                  !Decide if the normal vector should be reversed to be outward  
                  IF(ABS(localContactXiNormal)==localContactXiNormal) THEN
                    REVERSE_NORMAL=.FALSE.
                  ELSE
                    REVERSE_NORMAL=.TRUE. 
                  ENDIF                 
                  DO gaussPointIdx=1,INTERFACE_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
                    CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                      & INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR, &
                      & ERR,ERROR,*999)
                    CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(INTERFACE_GEOMETRIC_BASIS%NUMBER_OF_XI, &
                      & INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)% &
                      & PTR,ERR,ERROR,*999)
                    RWG=INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                      & JACOBIAN*INTERFACE_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gaussPointIdx) !Multiply by the interface mesh surface J since it's integrated over the surface
                    XI=INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY, &
                      & ELEMENT_NUMBER,interface_matrix_idx,gaussPointIdx,ERR,ERROR)                   
                    !Evaluate the gap at this Gauss point, note it's the distance between two points, not multiplied by the normal
                    gap=0 
                    DO interface_matrix_gap_idx=1,INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
                      elementNumberGap=INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY &
                        & (ELEMENT_NUMBER,interface_matrix_idx)%COUPLED_MESH_ELEMENT_NUMBER
                      XIGap=INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(INTERFACE_CONDITION% &
                        & INTERFACE%MESH_CONNECTIVITY,ELEMENT_NUMBER,interface_matrix_gap_idx,gaussPointIdx,ERR,ERROR)
                      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumberGap, &
                          & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_gap_idx)% &
                          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                      CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,XIGap,INTERFACE_EQUATIONS% &
                        & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_gap_idx)%GEOMETRIC_INTERPOLATION(1)% &
                        & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)!Needed for computing metrics
                      IF (interface_matrix_gap_idx==1) THEN
                        gap=gap+INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_gap_idx)% &
                          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(:,1)
                      ELSE
                        gap=gap-INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_gap_idx)% &
                          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(:,1)
                      ENDIF
                    ENDDO
                    DO rowComponentIdx=1,coupledMeshNumberOfGeometricComponents
                      rowMeshComponentNumber=INTERFACE_MATRIX_VARIABLE%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER  
                      coupledMeshNumberOfXi=COUPLED_MESH_GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)% &
                        & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS%NUMBER_OF_XI    
                      SELECT CASE(coupledMeshNumberOfXi)
                      CASE(2)
                        !Compute the other xi direction of the line normal, i.e. xi position input for calculating interpolated point
                        XI_POSITION(1)=XI(OTHER_XI_DIRECTIONS2(ABS(localContactXiNormal)))
                        !Compute interpolation parameters for the line for the couple mesh field
                        !Needs the global contact line/face number.
                        CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,globalContactNumber, &
                          & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)!Needed for computing metrics
                      CASE(3)
                        SELECT CASE(ABS(localContactXiNormal))
                        CASE(1)
                          !Compute the other xi directions of the face normal, i.e. xi positions input for calculating interpolated point
                          XI_POSITION(1)=XI(2)
                          XI_POSITION(2)=XI(3)
                        CASE(2)
                          XI_POSITION(1)=XI(1)
                          XI_POSITION(2)=XI(3)
                        CASE(3)
                          XI_POSITION(1)=XI(1)
                          XI_POSITION(2)=XI(2)
                        END SELECT
                        !Compute interpolation parameters for the face for the couple mesh field
                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,globalContactNumber, &
                          & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                      CASE DEFAULT
                        CALL FLAG_ERROR("Xi dimension of coupled mesh should be <= 3 and >=0.",ERR,ERROR,*999)
                      END SELECT    
                      CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,XI_POSITION(1:coupledMeshNumberOfXi-1),INTERFACE_EQUATIONS% &
                        & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%GEOMETRIC_INTERPOLATION(1)% &
                        & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)!Needed for computing metrics
                      CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(coupledMeshNumberOfXi,INTERFACE_EQUATIONS% &
                        & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%GEOMETRIC_INTERPOLATION(1)% &
                        & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999) ! for compute normal         
                      INTERPOLATED_POINT_METRICS=>INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                        & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR                 
                      !Calculate outward normal of the contact surface in the coupled mesh
                      CALL FIELD_POSITION_NORMAL_TANGENTS_CALCULATE_INT_PT_METRIC(INTERPOLATED_POINT_METRICS, &
                        & REVERSE_NORMAL,POSITION(1:coupledMeshNumberOfGeometricComponents), &
                        & NORMAL(1:coupledMeshNumberOfGeometricComponents),TANGENTS(1:coupledMeshNumberOfGeometricComponents, &
                        & 1:coupledMeshNumberOfXi),ERR,ERROR,*999)   
                      COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                      !Gap function test.  
                      IF (interface_matrix_idx==1) THEN
                        penetration=gap(rowComponentIdx)*NORMAL(rowComponentIdx) 
                      ELSE
                        penetration=-gap(rowComponentIdx)*NORMAL(rowComponentIdx) 
                      ENDIF
                      !IF (penetration>0.0_DP) THEN !Only add the contribution of the Gauss point if it's a penetration
                        DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS 
                          rowIdx=rowParameterIdx+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          
!                          IF (interface_matrix_idx==1) THEN
!                            PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,rowParameterIdx,NO_PART_DERIV,XI,ERR,ERROR)* &
!                            & 1!NORMAL(rowComponentIdx)
!                          ELSE
!                            PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,rowParameterIdx,NO_PART_DERIV,XI,ERR,ERROR)* &
!                            & -1!NORMAL(rowComponentIdx)
!                          ENDIF
                          PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,rowParameterIdx,NO_PART_DERIV,XI,ERR,ERROR)* &
                            & NORMAL(rowComponentIdx)
                          colComponentIdx=rowComponentIdx !Since x and y are not coupled.
                          DO colParameterIdx=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                            colIdx=colParameterIdx+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(colComponentIdx-1)
                            PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,gaussPointIdx)
                            INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)=INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)+ &
                              & PGNSI*PGMSI*RWG!*INTERPOLATED_POINT_METRICS%JACOBIAN
                          ENDDO !colParameterIdx
                      ENDDO !rowParameterIdx  
                     ! ENDIF
                    ENDDO !rowComponentIdx       
                  ENDDO !gaussPointIdx
                  
                  !scale factor update for row variable
                  CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(coupledMeshElementNumber,INTERFACE_EQUATIONS% &
                      & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                      & INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  IF(INTERFACE_MATRIX_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    rowIdx=0
                    DO rowComponentIdx=1,coupledMeshNumberOfGeometricComponents
                      DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS 
                        rowIdx=rowIdx+1
                        colComponentIdx=rowComponentIdx !X and Y are decoupled
                        DO colParameterIdx=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                          colIdx=colParameterIdx+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(colComponentIdx-1)
                          INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)=INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx, &
                            & colIdx)*INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                            & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)% &
                            & PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)           
                        ENDDO !colParameterIdx
                      ENDDO !rowParameterIdx 
                    ENDDO !rowComponentIdx
                  ENDIF
                  !scale factor update column variable
                  CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,INTERFACE_EQUATIONS% &
                      & INTERPOLATION%INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
                      & INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  IF(INTERFACE_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                    rowIdx=0
                    DO rowComponentIdx=1,coupledMeshNumberOfGeometricComponents
                      DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS 
                        rowIdx=rowIdx+1
                        colComponentIdx=rowComponentIdx !X and Y are decoupled
                        DO colParameterIdx=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                          colIdx=colParameterIdx+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(colComponentIdx-1)
                          INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)=INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx, &
                            & colIdx)*INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(1)% &
                            & PTR%SCALE_FACTORS(colParameterIdx,colComponentIdx)              
                        ENDDO !colParameterIdx
                      ENDDO !rowParameterIdx 
                    ENDDO !rowComponentIdx
                  ENDIF
                  !Evaluate RHS element vector for frictionless contact
                  RHS_VEC=>INTERFACE_EQUATIONS%INTERFACE_MATRICES%RHS_VECTOR
                  DO colComponentIdx=1,InterfaceNumberOfGeometricComponent
                    rowComponentIdx=colComponentIdx
                    rowMeshComponentNumber=INTERFACE_MATRIX_VARIABLE%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER   
                    DO colParameterIdx=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      colIdx=colParameterIdx+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(colComponentIdx-1)
                      DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        rowIdx=rowParameterIdx+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                        rowGlobalDofNumber=INTERFACE_ELEMENT_MATRIX%ROW_DOFS(rowIdx)                      
                        dofValue=INTERFACE_MATRIX_DEPENDENT_FIELD%VARIABLES(1)%PARAMETER_SETS%PARAMETER_SETS(1)% &
                          & PTR%PARAMETERS%CMISS%DATA_DP(rowGlobalDofNumber)
                        !sign changes from + to - when flipping the matrices to RHS
                        !RHS_VEC%ELEMENT_VECTOR%VECTOR(colIdx)=RHS_VEC%ELEMENT_VECTOR%VECTOR(colIdx)- &
                        !  & INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)*dofValue
                        RHS_VEC%ELEMENT_VECTOR%VECTOR(colIdx)=0.0_DP !finite elasticity RHS=0   
                      ENDDO !rowParameterIdx
                    ENDDO !colParameterIdx
                  ENDDO !colComponentIdx        
                               
                !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ CODE FINISH HERE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                
                !Points connectivity
                CASE(FIELD_DATA_POINT_BASED_INTERPOLATION) !Points connecitivity
                  IF(INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                      & interface_matrix_idx==INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN !Penalty matrix
                    DO dataPointIdx=1,DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%NUMBER_OF_PROJECTED_DATA
                      dataPointNumber=DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%DATA_INDICES(dataPointIdx)% &
                       & GLOBAL_NUMBER
                      DO colComponentIdx=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                        colIdx=dataPointIdx+DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)% &
                          & NUMBER_OF_PROJECTED_DATA*(colComponentIdx-1)
                        weight=POINTS_CONNECTIVITY%INTERFACE%DATA_POINTS%DATA_POINTS(dataPointNumber)% &
                          & weightS(colComponentIdx)
                        dofValue=INTERFACE_PENALTY_FIELD%VARIABLES(1)%PARAMETER_SETS%PARAMETER_SETS(1)%PTR%PARAMETERS%CMISS% &
                          & DATA_DP(colIdx)
                        INTERFACE_ELEMENT_MATRIX%MATRIX(colIdx,colIdx)=-(1.0_DP/dofValue)!*weight!*weight
                      ENDDO
                    ENDDO
                  ELSE
                    DO dataPointIdx=1,DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%NUMBER_OF_PROJECTED_DATA
                      dataPointNumber=DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%DATA_INDICES(dataPointIdx)% &
                       & GLOBAL_NUMBER        
                      IF(POINTS_CONNECTIVITY%DATA_POINT_PROJECTED(dataPointNumber)) THEN !If the point has been orthogonally projected
                        !Get the global number of the coupled mesh element that this point is projected on
                        coupledMeshElementNumber=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointNumber,interface_matrix_idx)% &
                          & COUPLED_MESH_ELEMENT_NUMBER
                        !Get the index of the coupled mesh element that this data point is projected on
                        coupledMeshElementIdx=1
                        DO WHILE (coupledMeshElementNumber/=POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER, &
                          & interface_matrix_idx)%ELEMENT_NUMBERS(coupledMeshElementIdx))
                          coupledMeshElementIdx=coupledMeshElementIdx+1
                        ENDDO       
                            
                        !Local number of the line in contact
                        localContactNumber=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointNumber,interface_matrix_idx)% &
                          & COUPLED_MESH_CONTACT_NUMBER
                                               
                        !The xi direction of the surface/line normal
                        localContactXiNormal=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointNumber,interface_matrix_idx)% &
                          & COUPLED_MESH_CONTACT_XI_NORMAL 
                        !Decide if the normal vector should be reversed to be outward  
                        IF(ABS(localContactXiNormal)==localContactXiNormal) THEN
                          REVERSE_NORMAL=.FALSE.
                        ELSE
                          REVERSE_NORMAL=.TRUE. 
                        ENDIF                                          
                          
                        DO rowComponentIdx=1,coupledMeshNumberOfGeometricComponents
                          rowMeshComponentNumber=INTERFACE_MATRIX_VARIABLE%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER  
                          coupledMeshNumberOfXi=COUPLED_MESH_GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)% &
                            & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS%NUMBER_OF_XI 
                          XI(1:coupledMeshNumberOfXi)=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointNumber, &
                            & interface_matrix_idx)%XI(1:coupledMeshNumberOfXi,rowMeshComponentNumber)   
                          SELECT CASE(coupledMeshNumberOfXi)
                          CASE(2)
                            !Global number of the line in contact
                            globalContactNumber=POINTS_CONNECTIVITY%INTERFACE%COUPLED_MESHES(interface_matrix_idx)%PTR% &
                              & DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)% &
                              & ELEMENT_LINES(localContactNumber)   
                            !Compute interpolation parameters for the line for the couple mesh field
                            CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,globalContactNumber, &
                              & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                              & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)!Needed for computing metrics
                            !Compute the other xi direction of the line normal, i.e. xi position input for calculating interpolated point
                            XI_POSITION(1)=XI(OTHER_XI_DIRECTIONS2(ABS(localContactXiNormal)))
                          CASE(3)
                            !Global number of the line in contact
                            globalContactNumber=POINTS_CONNECTIVITY%INTERFACE%COUPLED_MESHES(interface_matrix_idx)%PTR% &
                              & DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)% &
                              & ELEMENT_FACES(localContactNumber)   
                            !Compute interpolation parameters for the face for the couple mesh field
                            CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,globalContactNumber, &
                              & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                              & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                            SELECT CASE(ABS(localContactXiNormal))
                            CASE(1)
                            !Compute the other xi directions of the face normal, i.e. xi positions input for calculating interpolated point
                              XI_POSITION(1)=XI(2)
                              XI_POSITION(2)=XI(3)
                            CASE(2)
                              XI_POSITION(1)=XI(1)
                              XI_POSITION(2)=XI(3)
                            CASE(3)
                              XI_POSITION(1)=XI(1)
                              XI_POSITION(2)=XI(2)
                            END SELECT
                          CASE DEFAULT
                              CALL FLAG_ERROR("Xi dimension of coupled mesh should be <= 3 and >=0.",ERR,ERROR,*999)
                          END SELECT                                          
                   
                          CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,XI_POSITION(1:coupledMeshNumberOfXi-1),INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)!Needed for computing metrics
                          ALLOCATE(INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT_METRICS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated points metrics.",ERR,ERROR,*999)
                          NULLIFY(INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR)
                          CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                            
                          CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(coupledMeshNumberOfXi,INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)                      
                          INTERPOLATED_POINT_METRICS=>INTERFACE_EQUATIONS%INTERPOLATION% &
                            & VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR                 
                          !Calculate outward normal of the contact surface in the coupled mesh
                          CALL FIELD_POSITION_NORMAL_TANGENTS_CALCULATE_INT_PT_METRIC(INTERPOLATED_POINT_METRICS, &
                            & REVERSE_NORMAL,POSITION(1:coupledMeshNumberOfGeometricComponents), &
                            & NORMAL(1:coupledMeshNumberOfGeometricComponents),TANGENTS(1:coupledMeshNumberOfGeometricComponents, &
                            & 1:coupledMeshNumberOfXi),ERR,ERROR,*999)  
                          CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)  
                          COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN &
                            & (rowMeshComponentNumber)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS   
                          DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS                    
                            rowIdx=COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*POINTS_CONNECTIVITY% &
                              & COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER,interface_matrix_idx)%NUMBER_OF_COUPLED_MESH_ELEMENTS* &
                              & (rowComponentIdx-1)+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS* &
                              & (coupledMeshElementIdx-1)+rowParameterIdx
                            PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,rowParameterIdx,NO_PART_DERIV,XI,ERR,ERROR)* &
                              & POINTS_CONNECTIVITY%NORMAL(rowComponentIdx,dataPointIdx)*MATRIX_COEFFICIENT*-1.0_dp
                            ROWBASIS=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,rowParameterIdx,NO_PART_DERIV,XI,ERR,ERROR)
                            colComponentIdx=rowComponentIdx !Since x and y are not coupled.
                            colIdx=dataPointIdx+DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)% &
                              & NUMBER_OF_PROJECTED_DATA*(colComponentIdx-1)
                            weight=POINTS_CONNECTIVITY%INTERFACE%DATA_POINTS%DATA_POINTS(dataPointNumber)% &
                              & weightS(colComponentIdx)
                            INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)=INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx, &
                              & colIdx)+PGMSI  
                          ENDDO !rowParameterIdx
                        ENDDO !component_idx     
                      ENDIF !If the data point has been orthogonally projected      
                    ENDDO !dataPointIdx               
                    !scale factor update
                    IF(INTERFACE_MATRIX_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                      rowIdx=0
                      DO rowComponentIdx=1,coupledMeshNumberOfGeometricComponents
                      rowMeshComponentNumber=INTERFACE_MATRIX_VARIABLE%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER   
                        DO coupledMeshElementIdx=1,POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER,interface_matrix_idx)% &
                          & NUMBER_OF_COUPLED_MESH_ELEMENTS
                          coupledMeshElementNumber=POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER,interface_matrix_idx)% &
                            & ELEMENT_NUMBERS(coupledMeshElementIdx)
                          COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN &
                            & (rowMeshComponentNumber)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                          CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(coupledMeshElementNumber,INTERFACE_EQUATIONS% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                          DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS 
                            rowIdx=rowIdx+1
                            colComponentIdx=rowComponentIdx !X and Y are decoupled
                            DO dataPointIdx=1,DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%NUMBER_OF_PROJECTED_DATA
                              colIdx=dataPointIdx+DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)% &
                                & NUMBER_OF_PROJECTED_DATA*(colComponentIdx-1)
                              dataPointNumber=DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%DATA_INDICES(dataPointIdx)% &
                                & GLOBAL_NUMBER
                              IF(POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointNumber,interface_matrix_idx)% &
                                  & COUPLED_MESH_ELEMENT_NUMBER==coupledMeshElementNumber) THEN
                                INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)=INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx, &
                                  & colIdx)*INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                                  & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)% &
                                  & PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                              ENDIF
                            ENDDO !dataPointIdx
                          ENDDO !rowParameterIdx 
                        ENDDO !coupledMeshElementIdx
                      ENDDO !rowComponentIdx
                    ENDIF
                    !Evaluate RHS element vector for frictionless contact
                    RHS_VEC=>INTERFACE_EQUATIONS%INTERFACE_MATRICES%RHS_VECTOR
                    DO colComponentIdx=1,InterfaceNumberOfGeometricComponent
                      rowComponentIdx=colComponentIdx
                      rowMeshComponentNumber=INTERFACE_MATRIX_VARIABLE%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER   
                      DO dataPointIdx=1,DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%NUMBER_OF_PROJECTED_DATA
                        colIdx=dataPointIdx+DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)% &
                          & NUMBER_OF_PROJECTED_DATA*(colComponentIdx-1)
                        dataPointNumber=DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%DATA_INDICES(dataPointIdx)% &
                          & GLOBAL_NUMBER
                        coupledMeshElementNumber=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointNumber,interface_matrix_idx)% &
                          & COUPLED_MESH_ELEMENT_NUMBER
                        coupledMeshElementIdx=1
                        DO WHILE (coupledMeshElementNumber/=POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER, &
                            & interface_matrix_idx)%ELEMENT_NUMBERS(coupledMeshElementIdx))
                          coupledMeshElementIdx=coupledMeshElementIdx+1
                        ENDDO  
                        !junk=0.0 
                        COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN &
                          & (rowMeshComponentNumber)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS                
                        DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                          rowIdx=COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*POINTS_CONNECTIVITY% &
                            & COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER,interface_matrix_idx)%NUMBER_OF_COUPLED_MESH_ELEMENTS* &
                            & (rowComponentIdx-1)+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS* &
                            & (coupledMeshElementIdx-1)+rowParameterIdx
                          rowGlobalDofNumber=INTERFACE_ELEMENT_MATRIX%ROW_DOFS(rowIdx)                      
                          dofValue=INTERFACE_MATRIX_DEPENDENT_FIELD%VARIABLES(1)%PARAMETER_SETS%PARAMETER_SETS(1)% &
                            & PTR%PARAMETERS%CMISS%DATA_DP(rowGlobalDofNumber)
                          !sign changes from + to - when flipping the matrices to RHS
                          RHS_VEC%ELEMENT_VECTOR%VECTOR(colIdx)=RHS_VEC%ELEMENT_VECTOR%VECTOR(colIdx)- &
                            & INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)*dofValue     
                          !junk=junk+INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)*dofValue            
                          !RHS_VEC%ELEMENT_VECTOR%VECTOR(colIdx)=0 !finite elasticity RHS=0
                        ENDDO !rowParameterIdx
                      ENDDO !dataPointIdx  
                    ENDDO !colComponentIdx        
                  ENDIF !Penalty method
                END SELECT !points connectivity
              END SELECT !Frictionless contact
            ENDIF !UPDATE_MATRIX
          ENDDO ! interface_matrix_idx
          
          ! Add penentration test for point connectivity frictionless contact (penetration has been performed earlier for mesh con case)
          IF (INTERFACE_CONDITION%OPERATOR==INTERFACE_CONDITION_FRICTIONLESS_CONTACT_OPERATOR) THEN
            IF(INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD) THEN
              numberOfCoupledMesh=INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES-1
            ELSE
              numberOfCoupledMesh=INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
            ENDIF
            DO interface_matrix_idx=1,numberOfCoupledMesh
              IF(INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%UPDATE_MATRIX) THEN
                !Pointers to the interface_matrix_idx'th coupled mesh variables (rows of interface element matrix)
                INTERFACE_MATRIX_DEPENDENT_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION% &
                  & VARIABLE_INTERPOLATION(interface_matrix_idx)%DEPENDENT_FIELD
                INTERFACE_MATRIX_VARIABLE=>INTERFACE_EQUATIONS%INTERFACE_MAPPING% & 
                  & INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE
                INTERFACE_MATRIX_VARIABLE_TYPE=INTERFACE_MATRIX_VARIABLE%VARIABLE_TYPE
                INTERFACE_ELEMENT_MATRIX=>INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%ELEMENT_MATRIX                   
                COUPLED_MESH_GEOMETRIC_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                  & GEOMETRIC_FIELD  
                INTERFACE_MATRIX_DEPENDENT_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION% &
                  & VARIABLE_INTERPOLATION(interface_matrix_idx)% DEPENDENT_FIELD     
                INTERFACE_MATRIX_VARIABLE=>INTERFACE_EQUATIONS%INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS &
                  & (interface_matrix_idx)%VARIABLE
                INTERFACE_MATRIX_VARIABLE_TYPE=INTERFACE_MATRIX_VARIABLE%VARIABLE_TYPE
                InterfaceNumberOfGeometricComponent=COUPLED_MESH_GEOMETRIC_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS          
                coupledMeshNumberOfGeometricComponents=INTERFACE_GEOMETRIC_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS    
                IF (INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD%VARIABLES(1)%COMPONENTS(1)%INTERPOLATION_TYPE== &
                    & FIELD_DATA_POINT_BASED_INTERPOLATION) THEN !Points connectivity
                  INTERFACE_GEOMETRIC_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%GEOMETRIC_FIELD      
                  POINTS_CONNECTIVITY=>INTERFACE_CONDITION%INTERFACE%POINTS_CONNECTIVITY
                  DATA_POINTS_TOPOLOGY=>POINTS_CONNECTIVITY%INTERFACE_MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR% &
                    & DOMAIN(INTERFACE_GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%DATA_POINTS 
                  DO colComponentIdx=1,InterfaceNumberOfGeometricComponent
                    rowComponentIdx=colComponentIdx
                    rowMeshComponentNumber=INTERFACE_MATRIX_VARIABLE%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER   
                    DO dataPointIdx=1,DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%NUMBER_OF_PROJECTED_DATA
                      colIdx=dataPointIdx+DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)% &
                        & NUMBER_OF_PROJECTED_DATA*(colComponentIdx-1)
                      dataPointNumber=DATA_POINTS_TOPOLOGY%ELEMENT_DATA_POINTS(ELEMENT_NUMBER)%DATA_INDICES(dataPointIdx)% &
                        & GLOBAL_NUMBER
                      coupledMeshElementNumber=POINTS_CONNECTIVITY%POINTS_CONNECTIVITY(dataPointNumber,interface_matrix_idx)% &
                        & COUPLED_MESH_ELEMENT_NUMBER
                      coupledMeshElementIdx=1
                      DO WHILE (coupledMeshElementNumber/=POINTS_CONNECTIVITY%COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER, &
                          & interface_matrix_idx)%ELEMENT_NUMBERS(coupledMeshElementIdx))
                        coupledMeshElementIdx=coupledMeshElementIdx+1
                      ENDDO    
                      DO rowParameterIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        rowIdx=COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*POINTS_CONNECTIVITY% &
                          & COUPLED_MESH_ELEMENTS(ELEMENT_NUMBER,interface_matrix_idx)%NUMBER_OF_COUPLED_MESH_ELEMENTS* &
                          & (rowComponentIdx-1)+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS* &
                          & (coupledMeshElementIdx-1)+rowParameterIdx
                        IF(-RHS_VEC%ELEMENT_VECTOR%VECTOR(colIdx)<0.000000001_DP) THEN !If the gap is a separation
                          INTERFACE_ELEMENT_MATRIX%MATRIX(rowIdx,colIdx)=0.0_DP !Only needed if reproject in every iteration. Since initial penetration may be zero
                        ENDIF
                      ENDDO !rowParameterIdx
                    ENDDO !dataPointIdx  
                  ENDDO !colComponentIdx
                  
                ENDIF !point connectivity
              ENDIF !UPDATE_MATRIX
            ENDDO ! interface_matrix_idx
            RHS_VEC%ELEMENT_VECTOR%VECTOR=0.0_DP !finite elasticity RHS=0  
          ENDIF
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Interface condition method "//TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
            & " is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
          INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
          IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element interface matrices:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",INTERFACE_MATRICES% &
              & NUMBER_OF_INTERFACE_MATRICES,ERR,ERROR,*999)
            DO interface_matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",interface_matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",INTERFACE_MATRICES% & 
                & MATRICES(interface_matrix_idx)%PTR%UPDATE_MATRIX,ERR,ERROR,*999)
              IF(INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%UPDATE_MATRIX) THEN
                ELEMENT_MATRIX=>INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%ELEMENT_MATRIX
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                  & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                  & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                  & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                  & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                  & '(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
            ENDDO !interface_matrix_idx
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE()")
#endif
       
    CALL EXITS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  FUNCTION INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(INTERFACE_MESH_CONNECTIVITY,ELEMENT_NUMBER,COUPLED_MESH_NUMBER, &
    & GAUSS_POINT,ERR,ERROR)
  
    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface meshes connectivity
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !< Index of the interface mesh element number 
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_NUMBER !< Index of the coupled mesh to which the interface gausspoint should be transformed
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT !< Index to the gauss point which needs to be transformed
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(SIZE(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY &
      & (ELEMENT_NUMBER,COUPLED_MESH_NUMBER)%XI,1))
    !Local Variables
    INTEGER(INTG) :: ms

    CALL ENTERS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM",ERR,ERROR,*999)
    
    INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM=0.0_DP
    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(ELEMENT_NUMBER,COUPLED_MESH_NUMBER)%XI)) THEN
        DO ms = 1,INTERFACE_MESH_CONNECTIVITY%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
          !INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:)= INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:) + &
          !  & INTERFACE_MESH_CONNECTIVITY%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP &
          !  & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,GAUSS_POINT) * INTERFACE_MESH_CONNECTIVITY% & 
          !  & ELEMENT_CONNECTIVITY(ELEMENT_NUMBER,COUPLED_MESH_NUMBER)%XI(:,1,ms)
          INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:)= INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:) + &
          & INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH%GENERATED_MESH%REGULAR_MESH%BASES(1)%PTR%QUADRATURE%QUADRATURE_SCHEME_MAP &
            & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,GAUSS_POINT) * INTERFACE_MESH_CONNECTIVITY% & 
            & ELEMENT_CONNECTIVITY(ELEMENT_NUMBER,COUPLED_MESH_NUMBER)%XI(:,1,ms)         
        ENDDO
      ELSE
        CALL FLAG_ERROR("Coupled Mesh xi array not allocated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Mesh Connectivity is not associated.",ERR,ERROR,*999)
    ENDIF
     
    CALL EXITS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM")
    RETURN
999 CALL ERRORS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM",ERR,ERROR)
    CALL EXITS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM")
    RETURN
    
  END FUNCTION INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM
  
  !
  !================================================================================================================================
  !

  !> Apply translation to the interface associated dependent fields, only for frictionless contact, to introduce penetration
  SUBROUTINE INTERFACE_CONDITION_TRANSLATION_INCREMENT_APPLY(INTERFACE_CONDITION,ITERATION_NUMBER, &
      & MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION 
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dependent_variable_idx,node_idx,comp_idx
    INTEGER(INTG) :: VERSION_NUMBER,DERIVATIVE_NUMBER,MESH_COMPONENT_NUMBER,NODE_USER_NUMBER
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the dependent field to translate
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_APPLY",ERR,ERROR,*999)

    VERSION_NUMBER=1;
    DERIVATIVE_NUMBER=1;
    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
    DO dependent_variable_idx=1,INTERFACE_CONDITION%DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
      !Get the field from the interface dependent field variables
      FIELD=>INTERFACE_CONDITION%DEPENDENT%FIELD_VARIABLES(dependent_variable_idx)%PTR%FIELD
      !Get the mesh component number for this field
      MESH_COMPONENT_NUMBER=FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
      NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(INTERFACE_CONDITION%BOUNDARY_CONDITIONS, &
        & INTERFACE_CONDITION%DEPENDENT%FIELD_VARIABLES(dependent_variable_idx)%PTR,BOUNDARY_CONDITIONS_VARIABLE, &
        & ERR,ERROR,*999)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
        IF(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_FIXED_INCREMENTED)>0) THEN
          DO comp_idx=1,INTERFACE_CONDITION%INTERFACE%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
            !Loop through the nodes in the mesh
            DO node_idx=1,INTERFACE_CONDITION%INTERFACE%COUPLED_MESHES(dependent_variable_idx)%PTR% &
                & TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%NODES%NUMBER_OF_NODES
              NODE_USER_NUMBER=INTERFACE_CONDITION%INTERFACE%COUPLED_MESHES(dependent_variable_idx)%PTR% &
                & TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%NODES%NODES(node_idx)%USER_NUMBER
              ! Add corresponding translation to the nodal coordinates
              CALL FIELD_PARAMETER_SET_ADD_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,VERSION_NUMBER, &
                & DERIVATIVE_NUMBER,NODE_USER_NUMBER,comp_idx,INTERFACE_CONDITION%TRANSLATIONS(comp_idx,ITERATION_NUMBER, &
                & dependent_variable_idx),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_ADD_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,VERSION_NUMBER, &
                & DERIVATIVE_NUMBER,NODE_USER_NUMBER,comp_idx,INTERFACE_CONDITION%TRANSLATIONS(comp_idx,ITERATION_NUMBER, &
                & dependent_variable_idx),ERR,ERROR,*999)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDDO
      
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_APPLY")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_APPLY",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_APPLY")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_TRANSLATION_INCREMENT_APPLY
  
  !
  !================================================================================================================================
  !

  !> Set translation to the interface associated dependent fields at a specific load increment, only for frictionless contact, to introduce penetration
  SUBROUTINE INTERFACE_CONDITION_TRANSLATION_INCREMENT_SET(INTERFACE_CONDITION,ITERATION_NUMBER, &
      & DEPENDENT_FIELD_NUMBER,TRANSLATION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION 
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER
    INTEGER(INTG), INTENT(IN) :: DEPENDENT_FIELD_NUMBER
    REAL(DP), INTENT(IN) :: TRANSLATION(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dependent_variable_idx,node_idx,comp_idx
    INTEGER(INTG) :: VERSION_NUMBER,DERIVATIVE_NUMBER,MESH_COMPONENT_NUMBER
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the dependent field to translate
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ITERATION_NUMBER<=SIZE(INTERFACE_CONDITION%TRANSLATIONS,2)) THEN
        IF(DEPENDENT_FIELD_NUMBER<=SIZE(INTERFACE_CONDITION%TRANSLATIONS,3)) THEN
          IF(SIZE(TRANSLATION,1)==SIZE(INTERFACE_CONDITION%TRANSLATIONS,1)) THEN
            INTERFACE_CONDITION%TRANSLATIONS(:,ITERATION_NUMBER,DEPENDENT_FIELD_NUMBER)=TRANSLATION
          ELSE
            CALL FLAG_ERROR("Number of spatial components in translation does not match with the that for the field.", &
              & ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Coupled mesh/ dependent field number is incorrect.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Load increment number is incorrect.",ERR,ERROR,*999)
      ENDIF
      
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_SET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_TRANSLATION_INCREMENT_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_TRANSLATION_INCREMENT_SET

  !
  !================================================================================================================================
  !

  !>Finds and returns in INTERFACE_CONDITION a pointer to the interface condition identified by USER_NUMBER in the given INTERFACE. If no interface condition with that USER_NUMBER exists INTERFACE_CONDITION is left nullified.
  SUBROUTINE INTERFACE_CONDITION_USER_NUMBER_FIND(USER_NUMBER,INTERFACE,INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<The interface to find the interface condition in.
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<On return a pointer to the interface condition with the given user number. If no interface condition with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_condition_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
        CALL FLAG_ERROR("Interface condition is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(INTERFACE_CONDITION)
        IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS)) THEN
          interface_condition_idx=1
          DO WHILE(interface_condition_idx<=INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS.AND. &
            & .NOT.ASSOCIATED(INTERFACE_CONDITION))
            IF(INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interface_condition_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
              INTERFACE_CONDITION=>INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            ELSE
              interface_condition_idx=interface_condition_idx+1
            ENDIF
          ENDDO
        ELSE
          LOCAL_ERROR="The interface conditions on interface number "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" are not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises an interface conditions and deallocates all memory.
  SUBROUTINE INTERFACE_CONDITIONS_FINALISE(INTERFACE_CONDITIONS,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_CONDITIONS_TYPE), POINTER :: INTERFACE_CONDITIONS !<A pointer to the interface conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    
    CALL ENTERS("INTERFACE_CONDITIONS_FINALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(INTERFACE_CONDITIONS)) THEN
      DO WHILE(INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS>0)
        INTERFACE_CONDITION=>INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(1)%PTR
        CALL INTERFACE_CONDITION_DESTROY(INTERFACE_CONDITION,ERR,ERROR,*999)
      ENDDO
      IF(ASSOCIATED(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)) DEALLOCATE(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)
      DEALLOCATE(INTERFACE_CONDITIONS)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITIONS_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITIONS_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITIONS_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITIONS_FINALISE

  !
  !================================================================================================================================
  !
  
  !>Initialises an interface conditions for an interface.
  SUBROUTINE INTERFACE_CONDITIONS_INITIALISE(INTERFACE,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
     
    CALL ENTERS("INTERFACE_CONDITIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS)) THEN
        LOCAL_ERROR="Interface conditions is already associated for interface number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE%INTERFACE_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface interface conditions.",ERR,ERROR,*999)
        INTERFACE%INTERFACE_CONDITIONS%INTERFACE=>INTERFACE
        INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=0
        NULLIFY(INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITIONS_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITIONS_FINALISE(INTERFACE%INTERFACE_CONDITIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITIONS_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE INTERFACE_CONDITIONS_ROUTINES
