!> \file
!> \author Chris Bradley
!> \brief This module handles all elasticity routines.
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
!> Contributor(s): Xiani (Nancy) Yan
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

!> This module handles all elasticity class routines.
MODULE RIGID_BODY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  
  PUBLIC RigidBody_ApplyTransformation
  
  PUBLIC RigidBody_EquationsSetClassTypeSet
  
  PUBLIC RigidBody_EquationsSetSetup
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Apply transformation to the geometric field based on calculated rigid body dofs 
  SUBROUTINE RigidBody_ApplyTransformation(equationsSetRigid,equationsSetNodal,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSetRigid !<A pointer to the equations set where rigid body dofs are stored
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSetNodal !<A pointer to the equations set where nodal dependent dofs are stored
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: geometricField,dependentField,rigidDofField
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    INTEGER(INTG) :: noGeoComp,localNodeNumber,userNodeNumber
    INTEGER(INTG) :: componentIdx,versionIdx,derivativeIdx,domainNodeIdx
    REAL(DP) :: theta(3),centreOfMass(3),rotationMatrix(3,3),nodalParameters(3),nodalParametersTrans(3)
    LOGICAL :: sameBases
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("RigidBody_ApplyTransformation",err,error,*999)
    
    !\todo: XY-rigid deformable contact, this subroutine needs tidy-up when LHS is ready, default to be 3 component, need to generalise
    IF(ASSOCIATED(equationsSetRigid)) THEN
      IF(ASSOCIATED(equationsSetNodal)) THEN
        ! The centred geometric field stored in rigid body equations set (equationsSetNodal) will be transformed according to 
        ! the rigid body dofs defined in the 5-10 components of equationsSetRigid, and stored in dependent field of equationsSetNodal
        ! i.e. geometricField + RigidBodyMotion = dependentField (where the object is currently)
        geometricField=>equationsSetNodal%GEOMETRY%GEOMETRIC_FIELD
        dependentField=>equationsSetNodal%DEPENDENT%DEPENDENT_FIELD
        rigidDofField=>equationsSetRigid%DEPENDENT%DEPENDENT_FIELD
        noGeoComp=geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%NUMBER_OF_COMPONENTS
        ! Get the 6 dof for rigid body position 
        DO componentIdx=1,noGeoComp
          CALL FIELD_PARAMETER_SET_GET_CONSTANT(rigidDofField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx+4, &
            & centreOfMass(componentIdx),err,error,*999)
        ENDDO !componentIdx
        DO componentIdx=1,noGeoComp
          CALL FIELD_PARAMETER_SET_GET_CONSTANT(rigidDofField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx+7, &
            & theta(componentIdx),err,error,*999)
        ENDDO !componentIdx
        ! Construct rotation matrix
        rotationMatrix(1,1)=COS(theta(2))*COS(theta(3))
        rotationMatrix(1,2)=COS(theta(2))*SIN(theta(3))
        rotationMatrix(1,3)=-SIN(theta(2))
        rotationMatrix(2,1)=SIN(theta(1))*SIN(theta(2))*COS(theta(3))-COS(theta(1))*SIN(theta(3))
        rotationMatrix(2,2)=SIN(theta(1))*SIN(theta(2))*SIN(theta(3))+COS(theta(1))*COS(theta(3))
        rotationMatrix(2,3)=SIN(theta(1))*COS(theta(2))
        rotationMatrix(3,1)=COS(theta(1))*SIN(theta(2))*COS(theta(3))+SIN(theta(1))*SIN(theta(3))
        rotationMatrix(3,2)=COS(theta(1))*SIN(theta(2))*SIN(theta(3))-SIN(theta(1))*COS(theta(3))
        rotationMatrix(3,3)=COS(theta(1))*COS(theta(2))
        !**********************************************************************************************************************
        ! Transform the field
        ! Determine if the all components have the same mesh components/ bases
        DO componentIdx=1,noGeoComp-1
          IF(geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(componentIdx)%MESH_COMPONENT_NUMBER/= &
            & geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(componentIdx+1)%MESH_COMPONENT_NUMBER) sameBases=.FALSE.
        ENDDO
        IF(sameBases) THEN
          domain=>geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(1)%DOMAIN !Use the 1st component domain since they are the same for all components
          domainNodes=>domain%TOPOLOGY%NODES
          DO domainNodeIdx=1,domainNodes%NUMBER_OF_NODES
            localNodeNumber=domainNodes%NODES(domainNodeIdx)%LOCAL_NUMBER !localNodeNumber is the same as domainNodeIdx
            userNodeNumber=domainNodes%NODES(domainNodeIdx)%USER_NUMBER
            DO derivativeIdx=1,domainNodes%NODES(domainNodeIdx)%NUMBER_OF_DERIVATIVES
              DO versionIdx=1,domainNodes%NODES(domainNodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                DO componentIdx=1,noGeoComp !Get all component for a nodal derivative
                  CALL FIELD_PARAMETER_SET_GET_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                    & derivativeIdx,userNodeNumber,componentIdx,nodalParameters(componentIdx),err,error,*999)
                ENDDO !componentIdx
                !Rotate the nodal parameters
                nodalParametersTrans(1:noGeoComp)=MATMUL(rotationMatrix(1:noGeoComp,1:noGeoComp), &
                  & nodalParameters(1:noGeoComp))
                DO componentIdx=1,noGeoComp !Update all component for a nodal derivative
                  CALL FIELD_PARAMETER_SET_UPDATE_NODE(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                    & derivativeIdx,userNodeNumber,componentIdx,nodalParametersTrans(componentIdx),err,error,*999)
                  IF(derivativeIdx==1) THEN ! Translate nodal coordinate
                    CALL FIELD_PARAMETER_SET_ADD_NODE(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                      & derivativeIdx,userNodeNumber,componentIdx,centreOfMass(componentIdx),err,error,*999)
                  ENDIF !derivativeIdx==1
                ENDDO !componentIdx
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ELSE
          CALL FLAG_ERROR("Transformation for different component bases not implemented.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set storing rigid body dependent dofs is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set storing rigid body dofs is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("RigidBody_ApplyTransformation")
    RETURN
999 CALL ERRORS("RigidBody_ApplyTransformation",err,error)
    CALL EXITS("RigidBody_ApplyTransformation")
    RETURN 1
  END SUBROUTINE RigidBody_ApplyTransformation

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the problem type and subtype for a rigid body equation set class.
  SUBROUTINE RigidBody_EquationsSetClassTypeSet(equationsSet,equationsType,equationsSubtype,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: equationsType !<The equation type
    INTEGER(INTG), INTENT(IN) :: equationsSubtype !<The equation subtype
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("RigidBody_EquationsSetClassTypeSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsType)
      CASE(EQUATIONS_SET_NO_TYPE)
        SELECT CASE(equationsType)
        CASE(EQUATIONS_SET_NO_SUBTYPE)
          equationsSet%CLASS=EQUATIONS_SET_RIGID_BODY_CLASS
          equationsSet%TYPE=equationsType
          equationsSet%SUBTYPE=equationsSubtype
        CASE DEFAULT
          localError="Equations set equation subtype "//TRIM(NUMBER_TO_VSTRING(equationsSubtype,"*",err,error))// &
            & " is not valid for an rigid body equations set class."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(equationsType,"*",err,error))// &
          & " is not valid for an rigid body equations set class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",err,error,*999)
    ENDIF
       
    CALL EXITS("RigidBody_EquationsSetClassTypeSet")
    RETURN
999 CALL ERRORS("RigidBody_EquationsSetClassTypeSet",err,error)
    CALL EXITS("RigidBody_EquationsSetClassTypeSet")
    RETURN 1
  END SUBROUTINE RigidBody_EquationsSetClassTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets up the rigid body equations set class.
  SUBROUTINE RigidBody_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)
  
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to setup a rigid body equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error  !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDimensions,dependentFieldNumberOfComponents,geometricMeshComponent,independentFieldNumberOfComponents
    INTEGER(INTG) :: geometricScalingType
    INTEGER(INTG) :: compIdx
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("RigidBody_EquationsSetSetup",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      IF(equationsSet%SUBTYPE==EQUATIONS_SET_NO_SUBTYPE) THEN
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL RigidBody_EquationEquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a rigid body equation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          !------------------------------------------------- dependent field ----------------------------------------------------
          
          !\todo: XY rigid-doformable contact, ATM dependent field stores nodal parameters, should store 6 dof of rigid body motion
          ! when LHS mapping is ready
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet%DEPENDENT% &
                & DEPENDENT_FIELD,err,error,*999)
              CALL FIELD_LABEL_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricDecomposition, &
                & err,error,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,equationsSet%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)

              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,1,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"current geometry", &
                & err,error,*999)

              CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)

              !calculate number of components with one component for each dimension and one for pressure
              dependentFieldNumberOfComponents=numberOfDimensions
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & dependentFieldNumberOfComponents,err,error,*999)

              DO compIdx=1,dependentFieldNumberOfComponents
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,compIdx, &
                  & geometricMeshComponent,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,compIdx, &
                  & geometricMeshComponent,err,error,*999)
              ENDDO !compIdx

              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO compIdx=1,dependentFieldNumberOfComponents
                  !\todo: XY rigid-doformable contact, ATM dependent field stores nodal parameters therefore nodal interpolation, 
                  ! should store 6 dof of rigid body motion and constant interpolation when LHS mapping is ready
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !compIdx
                CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType,err,error,*999)
                !Other solutions not defined yet
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            ELSE
            !Check the user specified field
              CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,1,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, & 
                & err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)

              CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)

              dependentFieldNumberOfComponents=numberOfDimensions
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, & 
                & dependentFieldNumberOfComponents,err,error,*999)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a moving mesh Laplace equation"
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          !------------------------------------------------- independent field ----------------------------------------------------
          
          !\todo: XY rigid-doformable contact, ATM independent field stores rigid body centre of mass
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created independent field
              CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet%INDEPENDENT% &
                & INDEPENDENT_FIELD,err,error,*999)
              CALL FIELD_LABEL_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,"independent field",err,error,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                & err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,geometricDecomposition, &
                & err,error,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,equationsSet%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)

              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,1,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE], &
                & err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"centre of mass", &
                & err,error,*999)

              CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)

              !calculate number of components with one component for each dimension and one for pressure
              independentFieldNumberOfComponents=numberOfDimensions
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                & independentFieldNumberOfComponents,err,error,*999)

              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO compIdx=1,independentFieldNumberOfComponents
                  !\todo: XY rigid-doformable contact, ATM dependent field stores nodal parameters therefore nodal interpolation, 
                  ! should store 6 dof of rigid body motion and constant interpolation when LHS mapping is ready
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                ENDDO !compIdx
                CALL FIELD_SCALING_TYPE_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_NO_SCALING,err,error,*999)
                !Other solutions not defined yet
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            ELSE
            !Check the user specified field
              CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,1,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, & 
                & err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)

              CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)

              independentFieldNumberOfComponents=numberOfDimensions
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, & 
                & independentFieldNumberOfComponents,err,error,*999)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a moving mesh Laplace equation"
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          IF(equationsSetSetup%ACTION_TYPE==EQUATIONS_SET_SETUP_FINISH_ACTION) THEN
            !\todo: XY rigid-doformable contact, need to finish allocate equations for rigid body equations set class once LHS mapping
            ! is ready.
            ! Calculate the centre of mass
            CALL RigidBody_CentreOfMassCalculate(equationsSet,err,error,*999)
            ! Centre the field where the nodal parameters are stored, making the centre of mass the origin, and update the first 3 
            ! components of rigid body dofs (position of COM)
            CALL RigidBody_CentreField(equationsSet,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
          & " does not equal a rigid body no subtype."
        CALL FLAG_ERROR(localError,err,error,*999)
      ENDIF !subtype match
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("RigidBody_EquationsSetSetup")
    RETURN
999 CALL ERRORS("RigidBody_EquationsSetSetup",err,error)
    CALL EXITS("RigidBody_EquationsSetSetup")
    RETURN 1
  END SUBROUTINE RigidBody_EquationsSetSetup
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a rigid body equation class.
  SUBROUTINE RigidBody_EquationEquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("RigidBody_EquationEquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SUBTYPE)        
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          equationsSet%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(solutionMethod,"*",err,error))//" is invalid."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
          & " is not valid for a rigid body equation class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("RigidBody_EquationEquationsSetSolutionMethodSet")
    RETURN
999 CALL ERRORS("RigidBody_EquationEquationsSetSolutionMethodSet",err,error)
    CALL EXITS("RigidBody_EquationEquationsSetSolutionMethodSet")
    RETURN 1
  END SUBROUTINE RigidBody_EquationEquationsSetSolutionMethodSet
  
  !
  !================================================================================================================================
  !

  !>Evaluate the centre of mass of a rigid body based on its geometric field and store the centre of mass in independent field
  SUBROUTINE RigidBody_CentreOfMassCalculate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to calculate the centre of mass for the rigid body
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: geometricField,independentField
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(BASIS_TYPE), POINTER :: geometricBasis
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: geometricQuadratureScheme
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointsMetrics(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    INTEGER(INTG) :: elementIdx,gaussIdx,componentIdx
    INTEGER(INTG) :: meshComponentNumber,localElementNumber,noGeoComp
    REAL(DP) :: gaussWeight,mass,centreOfMass(3)
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("RigidBody_CentreOfMassCalculate",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD ! The geometric field where nodal parameters are stored
      IF(ASSOCIATED(geometricField)) THEN
        IF(geometricField%FIELD_FINISHED) THEN
          independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD ! Independent field stores centre of mass
          IF(ASSOCIATED(independentField)) THEN
            IF(independentField%FIELD_FINISHED) THEN
              !\todo: XY-default mesh component number to be one for the first geometric component, same for all geometric components
              meshComponentNumber=geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(1)%MESH_COMPONENT_NUMBER
              elementsMapping=>geometricField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR%MAPPINGS%ELEMENTS
              ! Initialise interpolated points
              NULLIFY(interpolationParameters)
              NULLIFY(interpolatedPoints)
              NULLIFY(interpolatedPointsMetrics)
              CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,interpolationParameters,err,error,*999)
              CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
                & FIELD_GEOMETRIC_COMPONENTS_TYPE)
              CALL FIELD_INTERPOLATED_POINTS_METRICS_INITIALISE(interpolatedPoints,interpolatedPointsMetrics,err,error,*999)
              interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
              pointMetrics=>interpolatedPointsMetrics(FIELD_U_VARIABLE_TYPE)%PTR
              noGeoComp=geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%NUMBER_OF_COMPONENTS
              mass=0.0_DP
              centreOfMass=0.0_DP
              DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
                localElementNumber=elementsMapping%DOMAIN_LIST(elementIdx)
                geometricBasis=>geometricField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
                  & ELEMENTS(localElementNumber)%BASIS
                geometricQuadratureScheme=>geometricBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                ! Get the interpolation parameter for geometric field element
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
                  & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                DO gaussIdx=1,geometricQuadratureScheme%NUMBER_OF_GAUSS
                  gaussWeight=geometricQuadratureScheme%GAUSS_WEIGHTS(gaussIdx)
                  ! Evaluate the gauss points position in the rigid body element
                  CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                    & interpolatedPoint,ERR,ERROR,*999)
                  ! Evaluate the volume at the Gauss point, i.e. Jacobian
                  CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)
                  mass=mass+gaussWeight*pointMetrics%JACOBIAN
                  centreOfMass(1:noGeoComp)=centreOfMass(1:noGeoComp)+interpolatedPoint%VALUES(1:noGeoComp,NO_PART_DERIV)* &
                    & gaussWeight*pointMetrics%JACOBIAN
                ENDDO !gaussIdx
              ENDDO !elementIdx
              centreOfMass(1:noGeoComp)=centreOfMass(1:noGeoComp)/mass
              ! Update independent field with the calculate centre of mass
              DO componentIdx=1,noGeoComp
                CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & componentIdx,centreOfMass(componentIdx),err,error,*999)
              ENDDO !componentIdx
              CALL FIELD_INTERPOLATED_POINTS_METRICS_FINALISE(interpolatedPointsMetrics,err,error,*999)
              CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParameters,err,error,*999)
              CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPoints,err,error,*999)
            ELSE
              CALL FLAG_ERROR("Independent field is not finished.",err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Independent field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Geometric field is not finished.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Geometric field is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equation set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("RigidBody_CentreOfMassCalculate")
    RETURN
999 CALL ERRORS("RigidBody_CentreOfMassCalculate",err,error)
    CALL EXITS("RigidBody_CentreOfMassCalculate")
    RETURN 1
  END SUBROUTINE RigidBody_CentreOfMassCalculate
  
  !
  !================================================================================================================================
  !

  !>Centre the geometric field and update first 3 components of rigid body motion with the current centre of mass
  SUBROUTINE RigidBody_CentreField(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to calculate the centre of mass for the rigid body
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: geometricField,independentField,rigidDofField
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    INTEGER(INTG) :: meshComponentNumber,noGeoComp,localNodeNumber,userNodeNumber
    INTEGER(INTG) :: domainNodeIdx,derivativeIdx,versionIdx,componentIdx
    REAL(DP) :: centreOfMass(3)
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("RigidBody_CentreField",err,error,*999)
    
    !\todo: XY- subroutine needs tidy-up when LHS mapping is ready
    IF(ASSOCIATED(equationsSet)) THEN
      geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD ! The geometric field where nodal parameters are stored
      IF(ASSOCIATED(geometricField)) THEN
        IF(geometricField%FIELD_FINISHED) THEN
          independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD ! Independent field stores centre of mass
          IF(ASSOCIATED(independentField)) THEN
            IF(independentField%FIELD_FINISHED) THEN
              !\todo: XY-default mesh component number to be one for the first geometric component, same for all geometric components
              meshComponentNumber=geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(1)%MESH_COMPONENT_NUMBER
              domain=>geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(1)%DOMAIN !\todo: XY default to use the 1st component domain, need to generalise
              IF(ASSOCIATED(domain)) THEN
                !--------------------------------------------- centre geometric field-----------------------------------------------
                domainNodes=>domain%TOPOLOGY%NODES
                noGeoComp=geometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%NUMBER_OF_COMPONENTS
                DO componentIdx=1,noGeoComp
                  CALL FIELD_PARAMETER_SET_GET_CONSTANT(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx, &
                    & centreOfMass(componentIdx),err,error,*999)
                ENDDO !componentIdx
                DO domainNodeIdx=1,domainNodes%NUMBER_OF_NODES !loop through all the local nodes
                  !Get user number for the nodes to be transformed
                  localNodeNumber=domainNodes%NODES(domainNodeIdx)%LOCAL_NUMBER !localNodeNumber is the same as domainNodeIdx
                  userNodeNumber=domainNodes%NODES(domainNodeIdx)%USER_NUMBER
                  derivativeIdx=1
                  DO versionIdx=1,domainNodes%NODES(domainNodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                    DO componentIdx=1,noGeoComp !Get all component for a nodal derivative
                      CALL FIELD_PARAMETER_SET_ADD_NODE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                        & derivativeIdx,userNodeNumber,componentIdx,-centreOfMass(componentIdx),err,error,*999)
                    ENDDO !componentIdx
                  ENDDO !versionIdx
                ENDDO !domainNodeIdx
                
                !---------------------------------- update rigid body dof with Centre of Mass --------------------------------------
!                equationsSet%region%parent_region%sub_regions(2)%ptr%equations_sets%equations_sets(1)
                rigidDofField=>equationsSet%REGION%PARENT_REGION%SUB_REGIONS(1)%PTR%EQUATIONS_SETS%EQUATIONS_SETS(1)%PTR% &
                  & DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(rigidDofField)) THEN
                  DO componentIdx=1,noGeoComp
                    CALL FIELD_PARAMETER_SET_UPDATE_CONSTANT(rigidDofField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & componentIdx+4,centreOfMass(componentIdx),err,error,*999)
                  ENDDO !componentIdx
                ELSE
                  CALL FLAG_ERROR("Dependent field storing rigid body dofs is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Independent field is not finished.",err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Independent field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Geometric field is not finished.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Geometric field is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equation set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("RigidBody_CentreField")
    RETURN
999 CALL ERRORS("RigidBody_CentreField",err,error)
    CALL EXITS("RigidBody_CentreField")
    RETURN 1
  END SUBROUTINE RigidBody_CentreField

END MODULE RIGID_BODY_ROUTINES

