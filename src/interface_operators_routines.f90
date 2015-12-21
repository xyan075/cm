!> \file
!> \author Chris Bradley
!> \brief This module contains all interface conditions operators routines.
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
!> Contributor(s): Xiani (Nancy) Yan, Thiranja Prasad Babarenda Gamage
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
MODULE INTERFACE_OPERATORS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_EQUATIONS_ROUTINES
  USE INTERFACE_MAPPING_ROUTINES
  USE INTERFACE_MATRICES_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE MATHS
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FieldContinuity_FiniteElementCalculate
  
  PUBLIC FrictionlessContact_FiniteElementCalculate
  
  PUBLIC FrictionlessContact_contactMetricsCalculate
  
  PUBLIC InterfaceContactMetrics_Initialise, InterfaceContactMetrics_Finalise
  
  PUBLIC InterfaceContactMetrics_IterationAddGeoTermSet
  
  PUBLIC InterfaceContactMetrics_RigidBodySet
  
  PUBLIC SolidFluidOperator_FiniteElementCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for field continuity operator 
  SUBROUTINE FieldContinuity_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: GaussPoint, rowComponentIdx, rowIdx, rowParameterIdx, colComponentIdx, colIdx, colParameterIdx
    INTEGER(INTG) :: rowMeshComponentNumber,derivativeIdx,derivative,localElementNode,interfaceNode,interfaceDerivative
    INTEGER(INTG) :: coupledMeshElementNumber,coupledMeshIdx,coupledMeshVariableType,lagrangeVariableType
    INTEGER(INTG) :: connectedLine,decompositionLineNumber,localLineNodeIdx,connectedFace,decompositionFaceNumber,localFaceNodeIdx
    REAL(DP) :: XI(3),rwg,PGMSI,PGNSI,matrixCoefficient
    TYPE(BASIS_TYPE), POINTER :: interfaceDependentBasis,coupledMeshBasis,interfaceGeometricBasis, &
      & interfacePenaltyBasis,interfaceConnectivityBasis
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: interfaceQuadratureScheme
    TYPE(FIELD_TYPE), POINTER :: coupledMeshDependentField,interfaceDependentField,interfaceGeometricField, &
      & interfacePenaltyField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: interfaceMatrixVariable,lagrangeVariable
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: interfaceElementMatrix
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE), POINTER :: interfaceInterpolation
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: elementConnectivity
    TYPE(DOMAIN_LINE_TYPE), POINTER :: coupledMeshDomainLine
    TYPE(DOMAIN_FACE_TYPE), POINTER :: coupledMeshDomainFace
    TYPE(VARYING_STRING) :: localError

    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData !<A pointer to the decomposition data point topology
    TYPE(BASIS_TYPE), POINTER :: coupledMeshDependentBasis
    TYPE(FIELD_TYPE), POINTER :: coupledMeshGeometricField
    INTEGER(INTG) :: meshComponentNumber,numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi, &
      & numberOfMatrixCoupledElements
    INTEGER(INTG) :: dataPointIdx,localElementNumber,matrixElementIdx
    INTEGER(INTG) :: matrixCoefficients(2),interfaceelementnumber

    CALL ENTERS("FieldContinuity_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FLAG_error("Interface condition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE_EQUATIONS)) CALL FLAG_error("Interface equations is not associated." &
      & ,err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE)) CALL FLAG_error("Interface is not associated.",err,error,*999)

    interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS

    SELECT CASE(interfaceCondition%METHOD)

    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FLAG_error("Not implemented.",err,error,*999)

    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)

      SELECT CASE(interfaceCondition%integrationType)

      CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
        !Pointers to interface variables (columns of interface element matrix)
        interfaceInterpolation=>interfaceEquations%INTERPOLATION%INTERFACE_INTERPOLATION
        interfaceGeometricField=>interfaceInterpolation%GEOMETRIC_FIELD
        interfaceDependentField=>interfaceInterpolation%DEPENDENT_FIELD
        interfaceGeometricBasis=>interfaceGeometricField%DECOMPOSITION%DOMAIN(interfaceGeometricField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        interfaceDependentBasis=>interfaceDependentField%DECOMPOSITION%DOMAIN(interfaceDependentField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        SELECT CASE(interfaceCondition%METHOD)
        CASE(INTERFACE_CONDITION_PENALTY_METHOD)
          interfacePenaltyField=>interfaceInterpolation%PENALTY_FIELD
          interfacePenaltyBasis=>interfacePenaltyField%DECOMPOSITION%DOMAIN(interfacePenaltyField% &
            & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,interfaceInterpolation% &
            & PENALTY_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        ENDSELECT
        !Integrate using the interface quadrature scheme
        interfaceQuadratureScheme=>interfaceGeometricBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        lagrangeVariable=>interfaceEquations%INTERFACE_MAPPING%LAGRANGE_VARIABLE
        lagrangeVariableType=lagrangeVariable%VARIABLE_TYPE
        !Get element interpolation parameters from the first geometric interpolation set (to get Jacobian for interface surface integral)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,interfaceInterpolation% &
          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
        matrixCoefficient=1.0_DP
        DO coupledMeshIdx=1,interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
            IF(coupledMeshIdx>1) THEN
              matrixCoefficient=-1.0_DP
            ENDIF 
            !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
            coupledMeshDependentField=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_FIELD
            elementConnectivity=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(elementNumber,coupledMeshIdx)
            coupledMeshElementNumber=elementConnectivity%COUPLED_MESH_ELEMENT_NUMBER
            interfaceMatrixVariable=> &
              & interfaceEquations%INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(coupledMeshIdx)%VARIABLE
            coupledMeshVariableType=interfaceMatrixVariable%VARIABLE_TYPE
            interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
            interfaceConnectivityBasis=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%BASIS

            !coupledMeshDependentInterpolation=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
            !  & DEPENDENT_INTERPOLATION

            !Loop over gauss points
            DO GaussPoint=1,interfaceQuadratureScheme%NUMBER_OF_GAUSS
              !CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint, &
              !  & coupledMeshDependentInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR, &
              !  & err,error,*999)
              !CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
              !  & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)

              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              rwg=interfaceInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                & JACOBIAN*interfaceQuadratureScheme%GAUSS_WEIGHTS(GaussPoint)
              IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                  & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint,interfaceInterpolation% &
                  & PENALTY_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  !Loop over the Lagrange variable matrix rows
                  DO rowParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(rowParameterIdx,NO_PART_DERIV,GaussPoint)
                    rowIdx=rowIdx+1
                    interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,rowIdx)- &
                      & (1.0_DP/interfaceInterpolation%PENALTY_INTERPOLATION(1)% &
                      & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1))*PGNSI**2.0_DP*rwg
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ELSE
                !\todo defaults to first mesh component, generalise
                !XI=INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM( &
                !  & elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
                XI(1:interfaceDependentBasis%NUMBER_OF_XI)=INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM( &
                  & elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
                !XI=interfaceCondition%interface%pointsConnectivity%pointsConnectivity(GaussPoint,coupledMeshIdx)%xi
                ! Loop over number of Lagrange variable components as not all components in the dependent field variable may be coupled
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% & 
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS

                  SELECT CASE(interfaceDependentBasis%NUMBER_OF_XI)

                  CASE(1) !1D interface (line)
                    connectedLine = elementConnectivity%CONNECTED_LINE
                    decompositionLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(connectedLine)
                    coupledMeshDomainLine=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                      & LINES%LINES(decompositionLineNumber)
                    DO localLineNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_LINE(connectedLine)
                      localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                      DO derivativeIdx=1,coupledMeshDomainLine%BASIS%NUMBER_OF_DERIVATIVES(localLineNodeIdx)
                        derivative=coupledMeshBasis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                        derivative=coupledMeshDomainLine%DERIVATIVES_IN_LINE(1,derivativeIdx,localLineNodeIdx)
                        rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                        PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV,XI,err,error)
                        rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                        DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                          DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                            PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                            colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                            interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                              & PGNSI*PGMSI*rwg*matrixCoefficient
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivativeIdx
                    ENDDO !localLineNodeIdx

                  CASE(2) !2D interface (face)

                    SELECT CASE(coupledMeshBasis%NUMBER_OF_XI)

                    CASE(2) !Coupled Mesh has 2 xi directions
                      DO localElementNode=1,coupledMeshBasis%NUMBER_OF_NODES
                        DO derivative=1,coupledMeshBasis%NUMBER_OF_DERIVATIVES(localElementNode)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivative
                      ENDDO !localElementNode

                    CASE(3) !Coupled Mesh has 3 xi directions
                      connectedFace = elementConnectivity%CONNECTED_FACE
                      decompositionFaceNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(connectedFace)
                      coupledMeshDomainFace=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                        & FACES%FACES(decompositionFaceNumber)
                      DO localFaceNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                        localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                        DO derivativeIdx=1,coupledMeshDomainFace%BASIS%NUMBER_OF_DERIVATIVES(localFaceNodeIdx)
                          derivative=coupledMeshBasis% &
                            & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(derivativeIdx,localFaceNodeIdx,connectedFace)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivativeIdx
                      ENDDO !FaceNodeIdx

                    END SELECT !coupledMeshBasis%NUMBER_OF_XI

                  END SELECT !interfaceDependentBasis%NUMBER_OF_XI

                ENDDO !rowComponentIdx
              ENDIF
            ENDDO !GaussPoint

            !Scale factor adjustment
            !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix contribution to the residual for non-linear problems
            !\todo update looping of variables/components for non-zero matrix elements as done above 
            IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
              & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(elementNumber, &
                  & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  !Loop over element Lagrange variable rows
                  DO rowParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    interfaceElementMatrix%MATRIX(rowIdx,rowIdx)=interfaceElementMatrix%MATRIX(rowIdx,rowIdx) * &
                      & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)% &
                      & INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)**2
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ELSE
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(elementNumber, &
                  & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% & 
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx) * &
                        & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)% &
                        & INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR%SCALE_FACTORS(colParameterIdx,colComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
              !Scale factor adjustment for the row dependent variable
              IF(coupledMeshDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(coupledMeshElementNumber, &
                  & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                  & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR,err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,interfaceMatrixVariable%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)* &
                        & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                        & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR% &
                        & SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! coupledMeshIdx

      CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)

        matrixCoefficients(1)=1; !\todo: Change to interface mapping matrix coefficients
        matrixCoefficients(2)=-1;
        interfaceElementNumber = elementNumber!todo simplify
        interface=>interfaceCondition%INTERFACE
        pointsConnectivity=>interface%pointsConnectivity
        numberOfInterfaceMeshXi=pointsConnectivity%interfaceMesh%NUMBER_OF_DIMENSIONS
        IF(ASSOCIATED(pointsConnectivity)) THEN
          decompositionElementData=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD%DECOMPOSITION%TOPOLOGY%dataPoints% &
            & elementDataPoint(interfaceElementNumber)
         
          !Calculate PGSMI, update interface matrices with PGSMI, and update scale factors
          DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
            IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
              numberOfMatrixCoupledElements=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                & numberOfCoupledElements
              numberOfCoupledMeshXi=interface%COUPLED_MESHES(coupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
              coupledMeshGeometricField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                & GEOMETRY%GEOMETRIC_FIELD
              coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                & DEPENDENT%DEPENDENT_FIELD

              numberOfCoupledMeshGeoComp=coupledMeshGeometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%NUMBER_OF_COMPONENTS
              interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
              !mesh component number is the same for all geometric components in elasticity problems
              meshComponentNumber=coupledMeshDependentField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(1)% &
                & MESH_COMPONENT_NUMBER
              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                  & coupledMeshElementNumber

                !Calculate the element index (non-conforming element) for this interface matrix
                matrixElementIdx=1
                DO WHILE ((localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                    & elementNumbers(matrixElementIdx)).AND.(matrixElementIdx/= &
                    & pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)%numberOfCoupledElements))
                  matrixElementIdx=matrixElementIdx+1
                ENDDO   
                xi(1:numberOfCoupledMeshXi)=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                  & xi(1:numberOfCoupledMeshXi)
                !Calculate PGSMI for each data point component
                coupledMeshDependentBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
                DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                  DO rowParameterIdx=1,coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    PGMSI=BASIS_EVALUATE_XI(coupledMeshDependentBasis,rowParameterIdx,NO_PART_DERIV, &
                      & xi(1:numberOfCoupledMeshXi),ERR,ERROR)*matrixCoefficients(coupledMeshIdx)
                    rowIdx=rowParameterIdx+coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                    colIdx=dataPointIdx+decompositionElementData%numberOfProjectedData*(rowComponentIdx-1)
                    interfaceElementMatrix%MATRIX(rowIdx,colIdx)=PGMSI !Update interface element matrix with contact point contribution
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDDO !dataPointIdx

              !scale factor update
              IF(coupledMeshDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                  localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                    & coupledMeshElementNumber
                  CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(localElementNumber,interfaceEquations% &
                    & INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
                    & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  !Calculate the element index (non-conforming element) for this interface matrix
                  matrixElementIdx=1
                  DO WHILE ((localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                      & elementNumbers(matrixElementIdx)).AND.(matrixElementIdx/= &
                      & pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)%numberOfCoupledElements))
                    matrixElementIdx=matrixElementIdx+1
                  ENDDO
                  coupledMeshDependentBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
                  DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                    DO rowParameterIdx=1,coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                      rowIdx=rowParameterIdx+coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                      colIdx=dataPointIdx+decompositionElementData%numberOfProjectedData*(rowComponentIdx-1)
                      interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)* &
                        & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                        & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                        & PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                    ENDDO !rowParameterIdx
                  ENDDO !rowComponentIdx
                ENDDO !dataPointIdx
              ENDIF !.NOT. FIELD_NO_SCALING

            ENDIF !UPDATE_MATRIX
          ENDDO !coupledMeshIdx
        ELSE
          CALL FLAG_ERROR("Interface points connectivity is not associated.",err,error,*999)
        ENDIF

      CASE DEFAULT
        localError="Interface condition integration type "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType, &
          & "*",err,error))// " is not valid."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT !interfaceCondition%integrationType

    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FLAG_error("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Interface condition method "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%METHOD,"*",err,error))// &
        & " is not valid."
      CALL FLAG_error(localError,err,error,*999)
    END SELECT

    CALL EXITS("FieldContinuity_FiniteElementCalculate")
    RETURN
999 CALL ERRORS("FieldContinuity_FiniteElementCalculate",err,error)
    CALL EXITS("FieldContinuity_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE FieldContinuity_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for frictionless contact operator
  SUBROUTINE FrictionlessContact_FiniteElementCalculate(interfaceCondition,interfaceElementNumber,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: interfaceElementNumber !<The interface element number to calcualte the interface element matrix for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData !<A pointer to the decomposition data point topology
    TYPE(FIELD_TYPE), POINTER :: coupledMeshDependentField,penaltyField
    TYPE(INTERFACE_PENALTY_TYPE), POINTER :: interfacePenalty
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointsMetrics(:)
    TYPE(BASIS_TYPE), POINTER :: coupledMeshDependentBasis
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: interfaceElementMatrix
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: penaltyMatrix
    INTEGER(INTG) :: meshComponentNumber,numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi, &
      & numberOfMatrixCoupledElements,localDof
    INTEGER(INTG) :: dataPointIdx,coupledMeshIdx,xiIdx,localElementNumber,localFaceLineNumber,matrixElementIdx,rowComponentIdx, &
      & rowParameterIdx,rowIdx,colIdx,componentIdx,globalDataPointNumber
    INTEGER(INTG) :: matrixCoefficients(2)
    REAL(DP) :: PGMSI,contactStiffness
    REAL(DP) :: positionPoint(3),normalPoint(3),tangentsPoint(3,3),xi(3)
    LOGICAL :: reverseNormal
    REAL(DP), ALLOCATABLE :: gaps(:),gapsComponents(:,:),normals(:,:)
    LOGICAL, ALLOCATABLE :: orthogonallyProjected(:)
    
    
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("FrictionlessContact_FiniteElementCalculate",err,error,*999)
    
    IF(ASSOCIATED(interfaceCondition)) THEN
      interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
      IF(ASSOCIATED(interfaceEquations)) THEN
        interface=>interfaceCondition%INTERFACE
        IF(ASSOCIATED(interface)) THEN
          SELECT CASE(interfaceCondition%METHOD)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",err,error,*999)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            SELECT CASE(interfaceCondition%integrationType)
            CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
              CALL FLAG_ERROR("Mesh connectivity is not implemented for frictionless contact.",err,error,*999)
            CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
              matrixCoefficients(1)=1; !\todo: Change to interface mapping matrix coefficients
              matrixCoefficients(2)=-1;
              pointsConnectivity=>interface%pointsConnectivity
              numberOfInterfaceMeshXi=pointsConnectivity%interfaceMesh%NUMBER_OF_DIMENSIONS
              IF(ASSOCIATED(pointsConnectivity)) THEN
                decompositionElementData=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD%DECOMPOSITION%TOPOLOGY%dataPoints% &
                  & elementDataPoint(interfaceElementNumber)
                !###################################################################################################################
                
                !Test if datapoints were orthogonally projected.  
                !\todo: Allow the user to choose to only include orthogonally projected points or not (check is commented when populating element matrix below).  
                ALLOCATE(orthogonallyProjected(decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate orthogonal projected logicals.",err,error,*999)
                orthogonallyProjected=.TRUE. !Initialise orthogonal projected logicals
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                    & DEPENDENT%DEPENDENT_FIELD
                  !mesh component number is the same for all geometric components in elasticity problems
                  DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                    globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
                    DO xiIdx=1,SIZE(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%reducedXi,1)
                      IF(ABS(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%reducedXi(xiIdx)) &
                          & < ZERO_TOLERANCE) THEN
                        orthogonallyProjected(dataPointIdx)=.FALSE.
                      ENDIF
                    ENDDO !xiIdx
                  ENDDO !dataPointIdx
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Allocate memory for local allocatable variables
                ALLOCATE(gaps(decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate gaps.",err,error,*999)
                gaps=0.0_DP !Initialise gap functions
                ALLOCATE(gapsComponents(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component gaps.",err,error,*999)
                gapsComponents=0.0_DP !Initialise gap functions
                ALLOCATE(normals(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate normals.",err,error,*999)
                normals=0.0_DP !Initialise gap functions

                !Calculate Gap for each data point 
                !\todo: This is only required if only penetration is to penalized (ie seperation of meshes allowed.)
                ! If a no seperation condition is also required then calculation of the gap is not required.
                ! Need to allow user to choose which type of problem to solve.
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  coupledMeshDependentField=>interfaceCondition%DEPENDENT%FIELD_VARIABLES(coupledMeshIdx)%PTR%FIELD
                  numberOfCoupledMeshGeoComp=coupledMeshDependentField%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)% &
                    & PTR%NUMBER_OF_COMPONENTS
                  NULLIFY(interpolatedPoints)
                  NULLIFY(interpolationParameters)
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(coupledMeshDependentField,interpolationParameters,err,error, &
                    & *999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                  CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
                    & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                  interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
                  DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                    globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
                    !Only interpolate if orthogonally projected
                    !\todo: Allow the user to choose to only include orthogonally projected points or not (currenlty commented out).  
                    !IF(orthogonallyProjected(dataPointIdx)) THEN
                      localElementNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                        & coupledMeshElementNumber
                      localFaceLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)% &
                        & ELEMENT_FACES(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                        & elementLineFaceNumber)
                      SELECT CASE(numberOfInterfaceMeshXi) !Use face/line interpolation parameters for normal calculation
                      CASE(1)
                        CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      CASE(2)
                        SELECT CASE(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                            & elementLineFaceNumber)
                        CASE(1,3,5)
                          reverseNormal=.FALSE.
                        CASE(2,4,6)
                          reverseNormal=.TRUE.
                        END SELECT
                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      END SELECT
                      ! Determine the gap. 
                      ! \todo: Note that FIELD_INTERPOLATE_XI(FIRST_PART_DERIV by default calculates NO_PART_DERIV too
                      ! and is used because the FIRST_PART_DERIV is also need for the surface normal calculation. However the
                      ! normal is only calculated for one of the coupled bodies so unnecessary computation. Need to generalize
                      ! FIELD_INTERPOLATE_XI to allow the user to specify which PART_DERIV to calculate.  
                      CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,pointsConnectivity%pointsConnectivity(globalDataPointNumber, &
                        & coupledMeshIdx)%reducedXi(:),interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
                      gapsComponents(1:numberOfCoupledMeshGeoComp,dataPointIdx)=gapsComponents(1:numberOfCoupledMeshGeoComp, &
                        & dataPointIdx)+interpolatedPoint%VALUES(1:numberOfCoupledMeshGeoComp,NO_PART_DERIV)* &
                        & matrixCoefficients(coupledMeshIdx) !Calculate 3 components gap function for each contact point
                      !Calculate surface normal (use 2nd coupled mesh surface normal)
                      !\todo: Allow the user to choose which surface normal to calculate or alternatively allow for a weighted average of the two.  
                      IF (coupledMeshIdx==2) THEN
                        CALL FIELD_INTERPOLATED_POINTS_METRICS_INITIALISE(interpolatedPoints,interpolatedPointsMetrics, &
                          & err,error,*999)
                        CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(numberOfCoupledMeshGeoComp,interpolatedPointsMetrics &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                        CALL FIELD_POSITION_NORMAL_TANGENTS_CALCULATE_INT_PT_METRIC(interpolatedPointsMetrics &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,reverseNormal,positionPoint,normalPoint,tangentsPoint,err,error,*999)
                        normals(1:numberOfCoupledMeshGeoComp,dataPointIdx)=normalPoint(1:numberOfCoupledMeshGeoComp)
                        CALL FIELD_INTERPOLATED_POINTS_METRICS_FINALISE(interpolatedPointsMetrics,err,error,*999)
                      ENDIF !coupledMeshIdx==1
                    !ENDIF !orthogonallyProjected(dataPointIdx)
                  ENDDO !dataPointIdx
                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParameters,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPoints,err,error,*999)
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Calcualte 1 component gap
                DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                  gaps(dataPointIdx)=DOT_PRODUCT(gapsComponents(1:numberOfCoupledMeshGeoComp,dataPointIdx), &
                    & normals(1:numberOfCoupledMeshGeoComp,dataPointIdx))
                ENDDO !dataPointIdx
                
                !###################################################################################################################
                
                !Calculate PGSMI, update interface matrices with PGSMI, and update scale factors
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
                    numberOfMatrixCoupledElements=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                      & numberOfCoupledElements
                    numberOfCoupledMeshXi=interface%COUPLED_MESHES(coupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
                    numberOfCoupledMeshGeoComp=interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR%GEOMETRY% &
                      & GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%NUMBER_OF_COMPONENTS
                    coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                      & DEPENDENT%DEPENDENT_FIELD
                    interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
                    DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                      globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
                      !\todo: Allow the user to choose gap tolerance or default to zero tolerance (currently commented out).  
                      !IF(gaps(dataPointIdx)>1.0E-10) THEN !Only add contact point contribution if the gap is a penetration
                        localElementNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                          & coupledMeshElementNumber
                        !Calculate the element index (non-conforming element) for this interface matrix
                        matrixElementIdx=1
                        DO WHILE (localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                            & elementNumbers(matrixElementIdx))
                          matrixElementIdx=matrixElementIdx+1
                        ENDDO
                        CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(localElementNumber,interfaceEquations% &
                          & INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
                          & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                        xi(1:numberOfCoupledMeshXi)=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                          & xi(1:numberOfCoupledMeshXi)                  
                        DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                          meshComponentNumber=coupledMeshDependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                          !Calculate PGSMI for each data point component
                          coupledMeshDependentBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR% &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
                          !\todo: Loop over the number of coupled mesh dependent basis element parameters on the contact face to save a bit of computation.
                          DO rowParameterIdx=1,coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                            PGMSI=BASIS_EVALUATE_XI(coupledMeshDependentBasis,rowParameterIdx,NO_PART_DERIV, &
                              & xi(1:numberOfCoupledMeshXi),ERR,ERROR)*normals(rowComponentIdx,dataPointIdx)* &
                              & matrixCoefficients(coupledMeshIdx)
                            rowIdx=coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*numberOfMatrixCoupledElements* &
                              & (rowComponentIdx-1)+coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS* &
                              & (matrixElementIdx-1)+rowParameterIdx
                            colIdx=dataPointIdx
                            !Update interface element matrix with contact point contribution
                            !\todo: Seperate multiplication of scale factors if required.  
                            interfaceElementMatrix%MATRIX(rowIdx,colIdx)=PGMSI*interfaceEquations%INTERPOLATION% &
                              & VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
                              & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                          ENDDO !rowParameterIdx
                        ENDDO !rowComponentIdx
                      !ENDIF !gaps(dataPointIdx)>ZERO_TOLERANCE
                    ENDDO !dataPointIdx

                    IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%FIRST_ASSEMBLY) &
                      & interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%FIRST_ASSEMBLY=.FALSE.
                  ENDIF !UPDATE_MATRIX
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Deallocate memory
                IF(ALLOCATED(orthogonallyProjected)) DEALLOCATE(orthogonallyProjected)
                IF(ALLOCATED(gapsComponents)) DEALLOCATE(gapsComponents)
                IF(ALLOCATED(gaps)) DEALLOCATE(gaps)
                
                !###################################################################################################################
                
                !Calculate penalty matrix if required
                IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD) THEN
                  interfacePenalty=>interfaceCondition%PENALTY
                  IF(ASSOCIATED(interfacePenalty)) THEN
                    penaltyField=>interfacePenalty%PENALTY_FIELD
                    IF(ASSOCIATED(penaltyField)) THEN
                      penaltyMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(interfaceEquations% &
                        & INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES)%PTR
                      IF(ASSOCIATED(penaltyMatrix)) THEN
                        IF(penaltyMatrix%FIRST_ASSEMBLY .AND. penaltyMatrix%UPDATE_MATRIX) THEN
                          DO componentIdx=1,penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%NUMBER_OF_COMPONENTS
                            SELECT CASE(penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                              & COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_CONSTANT(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                & componentIdx,contactStiffness,err,error,*999)
                              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                                penaltyMatrix%ELEMENT_MATRIX%MATRIX(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                              ENDDO !dataPointIdx
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              localDof=penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS(componentIdx)% &
                                & PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(interfaceElementNumber)
                              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                & localDof,contactStiffness,err,error,*999)
                              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                                penaltyMatrix%ELEMENT_MATRIX%MATRIX(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                              ENDDO !dataPointIdx
                            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                                localDof=penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS(componentIdx)% &
                                  & PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointIdx)
                                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & localDof,contactStiffness,err,error,*999)
                                penaltyMatrix%ELEMENT_MATRIX%MATRIX(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                              ENDDO !dataPointIdx
                            CASE DEFAULT
                              localError="The interpolation type for component number "// &
                                & TRIM(NUMBER_TO_VSTRING(componentIdx,"*",err,error))// &
                                & " of variable type "//TRIM(NUMBER_TO_VSTRING(FIELD_U_VARIABLE_TYPE,"*",err,error))// &
                                & " of field number "//TRIM(NUMBER_TO_VSTRING(penaltyField%USER_NUMBER,"*",err,error))//" is "// &
                                & TRIM(NUMBER_TO_VSTRING(penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS &
                                & (componentIdx)%INTERPOLATION_TYPE,"*", err,error))// " which is invalid for penalty field."
                              CALL FLAG_ERROR(localError,err,error,*999)
                            END SELECT
                          ENDDO !componentIdx
                        ENDIF              
                      ELSE
                        CALL FLAG_ERROR("Interface penalty matrix is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Interface penalty field is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface penalty is not associated.",err,error,*999)
                  ENDIF
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface points connectivity is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="Interface condition integration type "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType, &
                & "*",err,error))// " is not valid."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT !interfaceCondition%integrationType
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="Interface condition method "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%METHOD,"*",err,error))// &
              & " is not valid."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT !interfaceCondition%METHOD
        ELSE
          CALL FLAG_ERROR("Interface is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("FrictionlessContact_FiniteElementCalculate")
    RETURN
999 CALL ERRORS("FrictionlessContact_FiniteElementCalculate",err,error)
    CALL EXITS("FrictionlessContact_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE FrictionlessContact_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !
  
  !>Finalise contact point metrics 
  SUBROUTINE InterfaceContactMetrics_ContactPointFinalise(contactMetrics,err,error,*)

    !Argument variables
    TYPE(InterfaceContactPointMetricsType) :: contactMetrics !<A pointer to the interface contact metrics
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL ENTERS("InterfaceContactMetrics_ContactPointFinalise",err,error,*999)
    
    contactMetrics%signedGapNormal=0.0_DP
    contactMetrics%contactStiffness=0.0_DP
    contactMetrics%contactForce=0.0_DP
    IF(ALLOCATED(contactMetrics%normal)) DEALLOCATE(contactMetrics%normal)
    IF(ALLOCATED(contactMetrics%tangents)) DEALLOCATE(contactMetrics%tangents)
    IF(ALLOCATED(contactMetrics%tangentDerivatives)) DEALLOCATE(contactMetrics%tangentDerivatives)
    IF(ALLOCATED(contactMetrics%covariantMetricTensor)) DEALLOCATE(contactMetrics%covariantMetricTensor)
    IF(ALLOCATED(contactMetrics%contravariantMetricTensor)) DEALLOCATE(contactMetrics%contravariantMetricTensor)
    IF(ALLOCATED(contactMetrics%inverseA)) DEALLOCATE(contactMetrics%inverseA)
    
    CALL EXITS("InterfaceContactMetrics_ContactPointFinalise")
    RETURN
999 CALL ERRORS("InterfaceContactMetrics_ContactPointFinalise",err,error)
    CALL EXITS("InterfaceContactMetrics_ContactPointFinalise")
    RETURN 1
    
  END SUBROUTINE InterfaceContactMetrics_ContactPointFinalise
  
  !
  !================================================================================================================================
  !

  !>Initilise contact point metrics individually
  SUBROUTINE InterfaceContactMetrics_ContactPointInitialise(contactPointMetrics,numberOfGeometricComp,numberOfDimensions, &
       & err,error,*)

    !Argument variables
    TYPE(InterfaceContactPointMetricsType) :: contactPointMetrics !<A pointer to the individual contact point metrics
    INTEGER(INTG), INTENT(IN) :: numberOfGeometricComp, numberOfDimensions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
    
    CALL ENTERS("InterfaceContactMetrics_ContactPointInitialise",err,error,*998)
    
    ! Allocate memory for each interface contact metric
    ALLOCATE(contactPointMetrics%normal(numberOfGeometricComp),STAT=err)
    IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics - normal.",err,error,*998)
    ALLOCATE(contactPointMetrics%tangents(numberOfDimensions,numberOfGeometricComp),STAT=err)
    IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics - tangents.",err,error,*999)
    contactPointMetrics%contactStiffness=0.0_DP
    contactPointMetrics%signedGapNormal=0.0_DP
    contactPointMetrics%contactForce=0.0_DP
    contactPointMetrics%Jacobian=0.0_DP
    ALLOCATE(contactPointMetrics%tangentDerivatives(numberOfDimensions,numberOfDimensions, &
      & numberOfGeometricComp),STAT=err)
    IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics - tangent derivatives.",err,error,*999)
    ALLOCATE(contactPointMetrics%covariantMetricTensor(numberOfDimensions,numberOfDimensions), &
      & STAT=err)
    IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics - covariants.",err,error,*999)
    ALLOCATE(contactPointMetrics%contravariantMetricTensor(numberOfDimensions,numberOfDimensions), &
      & STAT=err)
    IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics - contravariants.",err,error,*999)
    ALLOCATE(contactPointMetrics%inverseA(numberOfDimensions,numberOfDimensions),STAT=err)
    IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics - inverseA.",err,error,*999)

    
    CALL EXITS("InterfaceContactMetrics_ContactPointInitialise")
    RETURN
999 CALL InterfaceContactMetrics_ContactPointFinalise(contactPointMetrics,dummyErr,dummyError,*998) 
998 CALL ERRORS("InterfaceContactMetrics_ContactPointInitialise",err,error)
    CALL EXITS("InterfaceContactMetrics_ContactPointInitialise")
    RETURN 1
    
  END SUBROUTINE InterfaceContactMetrics_ContactPointInitialise
  
  !
  !================================================================================================================================
  !
  
  !>Finalise contact points metrics 
  SUBROUTINE InterfaceContactMetrics_Finalise(contactMetrics,err,error,*)

    !Argument variables
    TYPE(InterfaceContactMetricsType), POINTER :: contactMetrics !<A pointer to the interface contact metrics
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: contactPtIdx
    
    IF(ASSOCIATED(contactMetrics)) THEN
      contactMetrics%numberOfContactPts=0
      contactMetrics%addGeometricTerm=.FALSE.
      contactMetrics%iterationGeometricTerm=0
      IF(ALLOCATED(contactMetrics%orthogonallyProjected)) DEALLOCATE(contactMetrics%orthogonallyProjected)
      IF(ALLOCATED(contactMetrics%inContact)) DEALLOCATE(contactMetrics%inContact)
      IF(ALLOCATED(contactMetrics%residualOriginal)) DEALLOCATE(contactMetrics%residualOriginal)
      IF(ALLOCATED(contactMetrics%residualPerturbed)) DEALLOCATE(contactMetrics%residualPerturbed)
      IF(ALLOCATED(contactMetrics%contactPointMetrics)) THEN
        DO contactPtIdx=1,SIZE(contactMetrics%contactPointMetrics,1)
          CALL InterfaceContactMetrics_ContactPointFinalise(contactMetrics%contactPointMetrics(contactPtIdx),err,error,*999)
        ENDDO
        DEALLOCATE(contactMetrics%contactPointMetrics)
      ENDIF
      DEALLOCATE(contactMetrics)
    ENDIF
    
    CALL EXITS("InterfaceContactMetrics_Finalise")
    RETURN
999 CALL ERRORS("InterfaceContactMetrics_Finalise",err,error)
    CALL EXITS("InterfaceContactMetrics_Finalise")
    RETURN 1
    
  END SUBROUTINE InterfaceContactMetrics_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initilise contact point metrics 
  SUBROUTINE InterfaceContactMetrics_Initialise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(INTERFACE_GEOMETRY_TYPE), POINTER :: geometry
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(InterfaceContactMetricsType), POINTER :: contactMetrics 
    INTEGER(INTG) :: numberOfGeometricComp, numberOfDimensions, pointIdx
    
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
    
    CALL ENTERS("InterfaceContactMetrics_Initialise",err,error,*998)
    
    IF(ASSOCIATED(interfaceCondition)) THEN
      interface=>interfaceCondition%INTERFACE
      IF(ASSOCIATED(interface)) THEN
        pointsConnectivity=>interface%pointsConnectivity
        IF(ASSOCIATED(pointsConnectivity)) THEN
          IF(ASSOCIATED(interfaceCondition%interfaceContactMetrics)) THEN
            CALL FLAG_ERROR("Contact metrics are already associated.",err,error,*998)
          ELSE
            geometry=>interfaceCondition%GEOMETRY
            IF(ASSOCIATED(geometry)) THEN
              geometricField=>geometry%GEOMETRIC_FIELD
              IF(ASSOCIATED(geometricField)) THEN
                IF(ALLOCATED(geometricField%VARIABLES)) THEN
                  ! Initialise interface contact metrices information
                  ALLOCATE(interfaceCondition%interfaceContactMetrics,STAT=err)
                  IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics.",err,error,*998)
                  contactMetrics=>interfaceCondition%interfaceContactMetrics
                  contactMetrics%numberOfContactPts=SIZE(pointsConnectivity%pointsConnectivity,1)
                  ! Get number of geometric components and number of mesh dimensions
                  numberOfGeometricComp=geometricField%VARIABLES(1)%NUMBER_OF_COMPONENTS
                  numberOfDimensions=geometricField%VARIABLES(1)%DIMENSION
                  ALLOCATE(contactMetrics%contactPointMetrics(contactMetrics%numberOfContactPts),STAT=err)
                  IF(err/=0) CALL FLAG_ERROR("Could not allocate interface contact metrics.",err,error,*999)
                  ALLOCATE(contactMetrics%orthogonallyProjected(contactMetrics%numberOfContactPts),STAT=err)
                  IF(err/=0) CALL FLAG_ERROR("Could not allocate orthogonally projected logical.",err,error,*999)
                  ALLOCATE(contactMetrics%inContact(contactMetrics%numberOfContactPts),STAT=err)
                  IF(err/=0) CALL FLAG_ERROR("Could not allocate in contact logical.",err,error,*999)
                  DO pointIdx=1,contactMetrics%numberOfContactPts
                    CALL InterfaceContactMetrics_ContactPointInitialise(contactMetrics%contactPointMetrics(pointIdx), &
                      & numberOfGeometricComp, numberOfDimensions,err,error,*999)
                  ENDDO !pointIdx
                  contactMetrics%iterationGeometricTerm=0
                  contactMetrics%addGeometricTerm=.FALSE.
                ELSE
                  CALL FLAG_ERROR("Interface geometric field variables are not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface geometric field is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface geometry is not associated.",err,error,*999)
            ENDIF
          ENDIF !contact matrices not associated
        ELSE
          CALL FLAG_ERROR("Interface points connectivity is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("InterfaceContactMetrics_Initialise")
    RETURN
999 CALL InterfaceContactMetrics_Finalise(interfaceCondition%interfaceContactMetrics,dummyErr,dummyError,*998) 
998 CALL ERRORS("InterfaceContactMetrics_Initialise",err,error)
    CALL EXITS("InterfaceContactMetrics_Initialise")
    RETURN 1
    
  END SUBROUTINE InterfaceContactMetrics_Initialise
  
  !
  !================================================================================================================================
  !

  !>Initilise contact point metrics 
  SUBROUTINE InterfaceContactMetrics_IterationAddGeoTermSet(interfaceCondition,iterationNumber,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The iteration where geometric term is added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceContactMetricsType), POINTER :: interfaceContactMetrics
    
    CALL ENTERS("InterfaceContactMetrics_IterationAddGeoTermSet",err,error,*999)
    
    IF(ASSOCIATED(interfaceCondition)) THEN
      IF(interfaceCondition%INTERFACE_CONDITION_FINISHED) THEN
        interfaceContactMetrics=>interfaceCondition%interfaceContactMetrics
        IF(ASSOCIATED(interfaceContactMetrics)) THEN
          interfaceContactMetrics%iterationGeometricTerm=iterationNumber
        ELSE
          CALL FLAG_ERROR("Interface contact metrics is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition is not finished.",err,error,*999)
      ENDIF
      
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("InterfaceContactMetrics_IterationAddGeoTermSet")
    RETURN
999 CALL ERRORS("InterfaceContactMetrics_IterationAddGeoTermSet",err,error)
    CALL EXITS("InterfaceContactMetrics_IterationAddGeoTermSet")
    RETURN 1
    
  END SUBROUTINE InterfaceContactMetrics_IterationAddGeoTermSet
  
  !
  !================================================================================================================================
  !

  !> Set rigid body forces
  SUBROUTINE InterfaceContactMetrics_RigidBodySet(interfaceCondition,forces,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    REAL(DP), INTENT(IN) :: forces(:) !<The external forces
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceContactMetricsType), POINTER :: interfaceContactMetrics
    
    CALL ENTERS("InterfaceContactMetrics_RigidBodySet",err,error,*999)
    
    IF(ASSOCIATED(interfaceCondition)) THEN
      IF(interfaceCondition%INTERFACE_CONDITION_FINISHED) THEN
        interfaceContactMetrics=>interfaceCondition%interfaceContactMetrics
        IF(ASSOCIATED(interfaceContactMetrics)) THEN
          interfaceContactMetrics%rigidBody%forces(:)=forces(:)
        ELSE
          CALL FLAG_ERROR("Interface contact metrics is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition is not finished.",err,error,*999)
      ENDIF
      
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("InterfaceContactMetrics_RigidBodySet")
    RETURN
999 CALL ERRORS("InterfaceContactMetrics_RigidBodySet",err,error)
    CALL EXITS("InterfaceContactMetrics_RigidBodySet")
    RETURN 1
    
  END SUBROUTINE InterfaceContactMetrics_RigidBodySet
  
  !
  !================================================================================================================================
  !

  !>Calculates the metrics of contact points for linearisation 
  SUBROUTINE FrictionlessContact_ContactMetricsCalculate(interfaceCondition,iterationNumber,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The current iteration number in non-linear solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(DATA_POINTS_TYPE), POINTER :: dataPoints
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(InterfaceContactMetricsType), POINTER :: contactMetrics 
    TYPE(InterfaceContactPointMetricsType), POINTER :: contactPointMetrics
    TYPE(FIELD_TYPE), POINTER :: projectedDependentField,penaltyField,slaveDependentField,slaveGeometricField,LagrangeField
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParametersMaster(:),interpolationParametersSlave(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPointsMaster(:),interpolatedPointsSlave(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointsMetricsMaster(:),interpolatedPointsMetricsSlave(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPointMaster,interpolatedPointSlave
    
    INTEGER(INTG) :: projectedMeshIdx,noGeoComp,noXi,localElementNumber,localFaceLineNumber,normalStiffnessComp,penaltyPtDof, &
      & ContPtElementNum,slaveMeshIdx
    INTEGER(INTG) :: contactPtIdx,xiIdx,i,j,componentIdx,countPts
    REAL(DP) :: gapsComponents(3),junkPosition(3),tangents(3,2), A(2,2), detA, centreOfMass(3)
    LOGICAL :: reverseNormal
    
    
    CALL ENTERS("FrictionlessContact_FiniteElementCalculate",err,error,*999)
    
    IF(ASSOCIATED(interfaceCondition)) THEN
      interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
      IF(ASSOCIATED(interfaceEquations)) THEN
        interface=>interfaceCondition%INTERFACE
        IF(ASSOCIATED(interface)) THEN
          dataPoints=>interface%DATA_POINTS
          IF(ASSOCIATED(dataPoints)) THEN
            pointsConnectivity=>interface%pointsConnectivity
            IF(ASSOCIATED(pointsConnectivity)) THEN
              contactMetrics=>interfaceCondition%interfaceContactMetrics
              IF(ASSOCIATED(contactMetrics)) THEN
                penaltyField=>interfaceCondition%PENALTY%PENALTY_FIELD
                IF(ASSOCIATED(penaltyField)) THEN
                  !Determine if geometric term is going to be added in this iteration
                  contactMetrics%addGeometricTerm=.FALSE.
                  IF((iterationNumber>=contactMetrics%iterationGeometricTerm) .AND. (contactMetrics%iterationGeometricTerm>0)) THEN
                    IF(interfaceCondition%operator==INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR) THEN
!                      contactMetrics%addGeometricTerm=.TRUE.
!                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"******************* Add Geom Term ********************",err,error,*999)
                    ENDIF 
                  ENDIF

                  !#################################################################################################################
                  
                  ! Get master dependent field and initalise interpolation points, parameters and metrics
                  projectedMeshIdx=2; ! The mesh where contact points are projected to, i.e. master mesh
                  
                  ! Get the dependent field for the master body
                  projectedDependentField=>interfaceCondition%DEPENDENT%FIELD_VARIABLES(projectedMeshIdx)%PTR%FIELD !master
                  noGeoComp=projectedDependentField%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)% &
                    & PTR%NUMBER_OF_COMPONENTS
                  noXi=2
                  NULLIFY(interpolationParametersMaster)
                  NULLIFY(interpolatedPointsMaster)
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(projectedDependentField,interpolationParametersMaster,err,error, &
                    & *999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                  CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParametersMaster,interpolatedPointsMaster,err,error,*999, &
                    & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                  CALL FIELD_INTERPOLATED_POINTS_METRICS_INITIALISE(interpolatedPointsMaster,interpolatedPointsMetricsMaster, &
                    & err,error,*999)
                  interpolatedPointMaster=>interpolatedPointsMaster(FIELD_U_VARIABLE_TYPE)%PTR
                  
                  slaveMeshIdx=1; 
                  
                  ! Get the current centre of mass for the master
                  slaveDependentField=>interfaceCondition%DEPENDENT%FIELD_VARIABLES(slaveMeshIdx)%PTR%FIELD
                  DO componentIdx=1,noGeoComp
                    CALL FIELD_PARAMETER_SET_GET_CONSTANT(slaveDependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & componentIdx+4,centreOfMass(componentIdx),err,error,*999)
                  ENDDO !componentIdx
                  

                  ! Gete geometric field for the slave
                  slaveGeometricField=>interfaceCondition%DEPENDENT%FIELD_VARIABLES(slaveMeshIdx)%PTR%FIELD%GEOMETRIC_FIELD !slave
                  
                  NULLIFY(interpolationParametersSlave)
                  NULLIFY(interpolatedPointsSlave)
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(slaveGeometricField,interpolationParametersSlave,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParametersSlave,interpolatedPointsSlave,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_METRICS_INITIALISE(interpolatedPointsSlave,interpolatedPointsMetricsSlave, &
                    & err,error,*999)
                  interpolatedPointSlave=>interpolatedPointsSlave(FIELD_U_VARIABLE_TYPE)%PTR
                  
!                  centreOfMass=0.0_DP
!                  countPts=0
                  !#################################################################################################################
                  
                  !Initialise contact logicals
                  IF(.NOT.dataPoints%DATA_PROJECTIONS(3)%PTR%perturbation) THEN
!                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"********************Reset******************",ERR,ERROR,*999)
                  contactMetrics%orthogonallyProjected=.TRUE. !Initialise orthogonal projected logicals
                  contactMetrics%inContact=.FALSE. !Initialise in contact logicals
                  ENDIF
                    
                  DO contactPtIdx=1,contactMetrics%numberOfContactPts !contactPtIdx is a global contact point index
                    IF(dataPoints%DATA_PROJECTIONS(3)%PTR%projectData(contactPtIdx)) THEN
                      ! Get the metric structure for this contact point
                      contactPointMetrics=>contactMetrics%contactPointMetrics(contactPtIdx)
                      ! Get the contact stiffness for this data point, so that it only need to access the data structure once
                      DO normalStiffnessComp=1,3 ! Xiani hardcode for penalising flexion of head
                        SELECT CASE(penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & COMPONENTS(normalStiffnessComp)%INTERPOLATION_TYPE)
                        CASE(FIELD_CONSTANT_INTERPOLATION)
                          CALL FIELD_PARAMETER_SET_GET_CONSTANT(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & normalStiffnessComp,contactPointMetrics%contactStiffness(normalStiffnessComp),err,error,*999)
                        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                          !\todo mesh 1 element no is assumed to be where contact point embed
    !                      ContPtElementNum=pointsConnectivity%pointsConnectivity(contactPtIdx,1)%coupledMeshElementNumber
    !                      WRITE(*,'(1X,''elem: '',I4)') ContPtElementNum
                          ContPtElementNum=pointsConnectivity%interfaceMesh%topology(1)%ptr%datapoints%datapoints(contactPtIdx)% &
                            & elementnumber
                          penaltyPtDof=penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS(normalStiffnessComp)% &
                            & PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ContPtElementNum)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & penaltyPtDof,contactPointMetrics%contactStiffness(normalStiffnessComp),err,error,*999)
                        CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                          ! \todo: 1 is for the frictionless contact, normal contact stiffness 
                          penaltyPtDof=penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS(normalStiffnessComp)% &
                            & PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(contactPtIdx)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                            & penaltyPtDof,contactPointMetrics%contactStiffness(normalStiffnessComp),err,error,*999)
                        CASE DEFAULT
                          CALL FLAG_ERROR("Interface penalty field can only be constant, element based or data point based.", &
                            & err,error,*999)
                        END SELECT
                      ENDDO
                      normalStiffnessComp=1
                      
                      ! Determine if a contact point has been orthogonally projected for this Newton step
                      
                      DO xiIdx=1,SIZE(pointsConnectivity%pointsConnectivity(contactPtIdx,projectedMeshIdx)%reducedXi,1)
                        IF((pointsConnectivity%pointsConnectivity(contactPtIdx,projectedMeshIdx)%reducedXi(xiIdx)==0.0_DP) .OR.  &
                            & (pointsConnectivity%pointsConnectivity(contactPtIdx,projectedMeshIdx)%reducedXi(xiIdx)==1.0_DP))THEN
                          contactMetrics%orthogonallyProjected(contactPtIdx)=.FALSE.
                        ENDIF
                      ENDDO !xiIdx
                    
                      IF(contactMetrics%orthogonallyProjected(contactPtIdx)) THEN !Only calculate metrics if orthogonally projected
                      
                        !#############################################################################################################
                        
                        ! Get the local element of the master where this point is projected on
                        localElementNumber=pointsConnectivity%pointsConnectivity(contactPtIdx,projectedMeshIdx)% &
                          & coupledMeshElementNumber
                        ! Get the local face number of the master
                        localFaceLineNumber=projectedDependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)% &
                          & ELEMENT_FACES(pointsConnectivity%pointsConnectivity(contactPtIdx,projectedMeshIdx)% &
                          & elementLineFaceNumber)
                        ! Get the interpolation parameter for evaluation of projected contact points on the master
                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParametersMaster(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                        ! Evaluate the projected data points on the master surface
                        CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,pointsConnectivity%pointsConnectivity(contactPtIdx, &
                          & projectedMeshIdx)%reducedXi(:),interpolatedPointMaster,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                        ! Calculate gap vector
                        gapsComponents(1:noGeoComp)=-dataPoints%DATA_POINTS(contactPtIdx)%position(1:noGeoComp)+ &
                          & interpolatedPointMaster%VALUES(1:noGeoComp,NO_PART_DERIV)
                        
                        !#############################################################################################################
                        
                        ! Get the local element of the slave where this point is defined
                        localElementNumber=pointsConnectivity%pointsConnectivity(contactPtIdx,slaveMeshIdx)%coupledMeshElementNumber
                        ! Get the local face number of the slave
                        localFaceLineNumber=slaveGeometricField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)% &
                          & ELEMENT_FACES(pointsConnectivity%pointsConnectivity(contactPtIdx,slaveMeshIdx)%elementLineFaceNumber)
                        ! Get the interpolation parameter for evaluation of contact points on the slave
                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParametersSlave(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                        ! Evaluate contact data point on the slave surface
                        CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,pointsConnectivity%pointsConnectivity(contactPtIdx, &
                          & slaveMeshIdx)%reducedXi(:),interpolatedPointSlave,err,error,*999)
                        CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(2,interpolatedPointsMetricsSlave &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999) !2=COORDINATE_JACOBIAN_AREA_TYPE
                        contactPointMetrics%Jacobian=interpolatedPointsMetricsSlave(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN
                        
                        !\todo: XY -rigid deformable contact, temporarily store contact points evaluated on rigid body (relative to centre
                        ! of mass) in the 3-components Lagrange field
                        LagrangeField=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD
                        IF(LagrangeField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%NUMBER_OF_COMPONENTS==noGeoComp) THEN
                          ! Store the new relative position of the contact points on the master body
                          DO componentIdx=1,noGeoComp
                            !contactPtIdx is the same as global number
                            CALL Field_ParameterSetUpdateDataPoint(LagrangeField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & contactPtIdx,componentIdx,interpolatedPointMaster%VALUES(componentIdx,NO_PART_DERIV)- &
                              & centreOfMass(componentIdx),ERR,ERROR,*999) 
                          ENDDO !componentIdx
                        ENDIF
                          
                        !#############################################################################################################
                        
                        ! Get the local element of the slave where this point is defined
                        localElementNumber=pointsConnectivity%pointsConnectivity(contactPtIdx,slaveMeshIdx)%coupledMeshElementNumber
                        ! Get the local face number of the slave
                        localFaceLineNumber=slaveGeometricField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)% &
                          & ELEMENT_FACES(pointsConnectivity%pointsConnectivity(contactPtIdx,slaveMeshIdx)%elementLineFaceNumber)
                        ! Get the interpolation parameter for evaluation of contact points on the slave
                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParametersSlave(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                        ! Evaluate contact data point on the slave surface
                        CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,pointsConnectivity%pointsConnectivity(contactPtIdx, &
                          & slaveMeshIdx)%reducedXi(:),interpolatedPointSlave,err,error,*999)
                        CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(2,interpolatedPointsMetricsSlave &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999) !2=COORDINATE_JACOBIAN_AREA_TYPE
                        contactPointMetrics%Jacobian=interpolatedPointsMetricsSlave(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN
                            
                        !#############################################################################################################
                        
                        ! Calculate normal and tangent vectors defined on the master surface, assumed 3D mesh in contact
                        ! Determine the sign of normal inward/outward
                        SELECT CASE(pointsConnectivity%pointsConnectivity(contactPtIdx,projectedMeshIdx)%elementLineFaceNumber)
                        CASE(1,3,5)
                          reverseNormal=.FALSE.
                        CASE(2,4,6)
                          reverseNormal=.TRUE.
                        END SELECT
                        CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(2,interpolatedPointsMetricsMaster &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999) !2=COORDINATE_JACOBIAN_AREA_TYPE
                        CALL FIELD_POSITION_NORMAL_TANGENTS_CALCULATE_INT_PT_METRIC(interpolatedPointsMetricsMaster &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,reverseNormal,junkPosition, &
                          & contactPointMetrics%normal,tangents,err,error,*999,.TRUE.)
                        ! Calculate signed gap
                        contactPointMetrics%signedGapNormal=DOT_PRODUCT(gapsComponents,contactPointMetrics%normal)
                        
                        ! \todo: 1 is for the frictionless contact, normal contact stiffness 
                        contactPointMetrics%contactForce=contactPointMetrics%signedGapNormal* &
                          & contactPointMetrics%contactStiffness(normalStiffnessComp)
                        IF(contactPointMetrics%signedGapNormal>ZERO_TOLERANCE) contactMetrics%inContact(contactPtIdx)=.TRUE.
  !                      IF(contactMetrics%inContact(contactPtIdx)) THEN
  !                        countPts=countPts+1
  !                        DO componentIdx=1,3
  !                          centreOfMass(componentIdx)=centreOfMass(componentIdx)+ &
  !                            & dataPoints%DATA_POINTS(contactPtIdx)%position(componentIdx)
  !                        ENDDO !componentIdx
  !                      ENDIF
                        !#############################################################################################################
                        
                        ! These terms are only required if geometric term is added to the contact stiffness matrix
                        IF((contactMetrics%addGeometricTerm) .AND. (contactMetrics%inContact(contactPtIdx))) THEN
                          ! Store the second derivative information
                          contactPointMetrics%tangentDerivatives(1,1,:)=interpolatedPointMaster%VALUES(1:noGeoComp,PART_DERIV_S1_S1)
                          contactPointMetrics%tangentDerivatives(1,2,:)=interpolatedPointMaster%VALUES(1:noGeoComp,PART_DERIV_S1_S2)
                          contactPointMetrics%tangentDerivatives(2,1,:)=contactPointMetrics%tangentDerivatives(1,2,:)
                          contactPointMetrics%tangentDerivatives(2,2,:)=interpolatedPointMaster%VALUES(1:noGeoComp,PART_DERIV_S2_S2)
                          
                          ! Re-populate tangent vectors into the appropriate format
                          contactPointMetrics%tangents(1,1:noGeoComp)=tangents(1:noGeoComp,1)
                          contactPointMetrics%tangents(2,1:noGeoComp)=tangents(1:noGeoComp,2)
                          ! Store the covariant and contravariant information
                          contactPointMetrics%covariantMetricTensor=interpolatedPointsMetricsMaster(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & GL(1:noXi,1:noXi)
                          contactPointMetrics%contravariantMetricTensor=interpolatedPointsMetricsMaster(FIELD_U_VARIABLE_TYPE)% &
                            & PTR%GU(1:noXi,1:noXi)
                          
                          ! Calculate inverse A
                          ! Calculate A first
                          A=0.0_DP !Initialise to be 0.0
                          DO i=1,noXi
                            DO j=1,noXi
                              A(i,j)=contactPointMetrics%covariantMetricTensor(i,j)+ DOT_PRODUCT(contactPointMetrics%normal, &
                                & contactPointMetrics%tangentDerivatives(i,j,:))*contactPointMetrics%signedGapNormal
                            ENDDO
                          ENDDO
                          CALL INVERT(A,contactPointMetrics%inverseA,detA,ERR,ERROR,*999)
                        ENDIF !add geometric term
                        !#############################################################################################################
                      ENDIF !orthogonally projected
                    ENDIF ! data projected
                  ENDDO !contactPtIdx
                  CALL FIELD_INTERPOLATED_POINTS_METRICS_FINALISE(interpolatedPointsMetricsMaster,err,error,*999)
                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParametersMaster,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPointsMaster,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_METRICS_FINALISE(interpolatedPointsMetricsSlave,err,error,*999)
                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParametersSlave,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPointsSlave,err,error,*999)
                  
                  
                  !#############################################################################################################
                  
!                  DO componentIdx=1,3
!                    centreOfMass(componentIdx)=centreOfMass(componentIdx)/countPts
!                    DO contactPtIdx=1,contactMetrics%numberOfContactPts !contactPtIdx is a global contact point index
!                      IF(contactMetrics%inContact(contactPtIdx)) THEN
!                        CALL Field_ParameterSetUpdateDataPoint(LagrangeField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!                          & contactPtIdx,componentIdx,dataPoints%DATA_POINTS(contactPtIdx)%position(componentIdx)- &
!                          & centreOfMass(componentIdx),ERR,ERROR,*999) 
!                        
!                        
!                      ENDIF
!                    ENDDO !contactPtIdx
!                  ENDDO !componentIdx
                  !#############################################################################################################
                ELSE
                  CALL FLAG_ERROR("Interface penalty field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface contact metrices is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface points connectivity is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface data points is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("FrictionlessContact_contactMetricsCalculate")
    RETURN
999 CALL ERRORS("FrictionlessContact_contactMetricsCalculate",err,error)
    CALL EXITS("FrictionlessContact_contactMetricsCalculate")
    RETURN 1
    
  END SUBROUTINE FrictionlessContact_contactMetricsCalculate

  !
  !================================================================================================================================
  !
  
  !>Calculates the element stiffness matrices for the given element number for solid fluid operator 
  !Note First interface matrix must be solid equations set's interface matrix
  SUBROUTINE SolidFluidOperator_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)
  
    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: GaussPoint, rowComponentIdx, rowIdx, rowParameterIdx, colComponentIdx, colIdx, colParameterIdx
    INTEGER(INTG) :: rowMeshComponentNumber,derivativeIdx,derivative,localElementNode,interfaceNode,interfaceDerivative
    INTEGER(INTG) :: coupledMeshElementNumber,coupledMeshIdx,coupledMeshVariableType,lagrangeVariableType
    INTEGER(INTG) :: connectedLine,decompositionLineNumber,localLineNodeIdx,connectedFace,decompositionFaceNumber,localFaceNodeIdx
    REAL(DP) :: XI(3),rwg,PGMSI,PGNSI,matrixCoefficient
    TYPE(BASIS_TYPE), POINTER :: interfaceDependentBasis,coupledMeshBasis,interfaceGeometricBasis, &
      & interfacePenaltyBasis,interfaceConnectivityBasis
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: interfaceQuadratureScheme
    TYPE(FIELD_TYPE), POINTER :: coupledMeshDependentField,interfaceDependentField,interfaceGeometricField, &
      & interfacePenaltyField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: interfaceMatrixVariable,lagrangeVariable
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: interfaceElementMatrix
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE), POINTER :: interfaceInterpolation
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: elementConnectivity
    TYPE(DOMAIN_LINE_TYPE), POINTER :: coupledMeshDomainLine
    TYPE(DOMAIN_FACE_TYPE), POINTER :: coupledMeshDomainFace
    TYPE(VARYING_STRING) :: localError

    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData !<A pointer to the decomposition data point topology
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointsMetrics(:)
    TYPE(BASIS_TYPE), POINTER :: coupledMeshDependentBasis
    TYPE(FIELD_TYPE), POINTER :: coupledMeshGeometricField
    INTEGER(INTG) :: meshComponentNumber,numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi, &
      & numberOfMatrixCoupledElements
    INTEGER(INTG) :: dataPointIdx,localElementNumber,localFaceLineNumber,matrixElementIdx
    INTEGER(INTG) :: matrixCoefficients(2),interfaceelementnumber

    CALL ENTERS("SolidFluidOperator_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FLAG_ERROR("Interface condition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE_EQUATIONS)) CALL FLAG_ERROR("Interface equations is not associated." &
      & ,err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE)) CALL FLAG_ERROR("Interface is not associated.",err,error,*999)

    interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS

    !===============================================================================================================================
    !Select Interface method
    SELECT CASE(interfaceCondition%METHOD)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FLAG_ERROR("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
      !=============================================================================================================================
      !Select Integration type
      SELECT CASE(interfaceCondition%integrationType)
      CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
        !Pointers to interface variables (columns of interface element matrix)
        interfaceInterpolation=>interfaceEquations%INTERPOLATION%INTERFACE_INTERPOLATION
        interfaceGeometricField=>interfaceInterpolation%GEOMETRIC_FIELD
        interfaceDependentField=>interfaceInterpolation%DEPENDENT_FIELD
        interfaceGeometricBasis=>interfaceGeometricField%DECOMPOSITION%DOMAIN(interfaceGeometricField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        interfaceDependentBasis=>interfaceDependentField%DECOMPOSITION%DOMAIN(interfaceDependentField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        !Integrate using the interface quadrature scheme
        interfaceQuadratureScheme=>interfaceGeometricBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        lagrangeVariable=>interfaceEquations%INTERFACE_MAPPING%LAGRANGE_VARIABLE
        !lagrangeVariableNumberOfComponents=>interfaceEquations%INTERFACE_MAPPING%LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
        lagrangeVariableType=lagrangeVariable%VARIABLE_TYPE
        !Get element interpolation parameters from the first geometric interpolation set (to get Jacobian for interface surface integral)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,interfaceInterpolation% &
          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
        matrixCoefficient=1.0_DP
        !Loop over interface matrices (1st solid, 2nd fluid)
        DO coupledMeshIdx=1,interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
            IF(coupledMeshIdx>1) THEN
              matrixCoefficient=-1.0_DP
            ENDIF 
            !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
            coupledMeshDependentField=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_FIELD
            elementConnectivity=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(elementNumber,coupledMeshIdx)
            coupledMeshElementNumber=elementConnectivity%COUPLED_MESH_ELEMENT_NUMBER
            interfaceMatrixVariable=> &
              & interfaceEquations%INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(coupledMeshIdx)%VARIABLE
            coupledMeshVariableType=interfaceMatrixVariable%VARIABLE_TYPE
            interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
            interfaceConnectivityBasis=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%BASIS

            !coupledMeshDependentInterpolation=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
            !  & DEPENDENT_INTERPOLATION
            
            !=======================================================================================================================
            !Loop over gauss points
            DO GaussPoint=1,interfaceQuadratureScheme%NUMBER_OF_GAUSS
              !CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint, &
              !  & coupledMeshDependentInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR, &
              !  & err,error,*999)
              !CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
              !  & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              !=====================================================================================================================
              !Interpolates field at given gauss point, includes first partial derivatives
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              !Calculates the interpolated point metrics and the associated interpolated point
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              !=====================================================================================================================
              ! R W G = GAUSSWEIGTHS * JACOBIAN
              rwg=interfaceInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                & JACOBIAN*interfaceQuadratureScheme%GAUSS_WEIGHTS(GaussPoint)
              IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                  & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              ELSE
                !===================================================================================================================
                !\todo defaults to first mesh component, generalise
                !TODO Originally XI=...
                XI(1:SIZE(elementConnectivity%XI,1))=INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM( &
                  & elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
                ! Loop over number of Lagrange variable components as not all components in the dependent field variable may be coupled
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS

                  SELECT CASE(interfaceDependentBasis%NUMBER_OF_XI)

                  CASE(1) !1D interface (line)
                    connectedLine=elementConnectivity%CONNECTED_LINE
                    decompositionLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(connectedLine)
                    coupledMeshDomainLine=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                      & LINES%LINES(decompositionLineNumber)
                    DO localLineNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_LINE(connectedLine)
                      localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                      DO derivativeIdx=1,coupledMeshDomainLine%BASIS%NUMBER_OF_DERIVATIVES(localLineNodeIdx)
                      !???????????????????
                        derivative=coupledMeshBasis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                        derivative=coupledMeshDomainLine%DERIVATIVES_IN_LINE(1,derivativeIdx,localLineNodeIdx)
                      !???????????????????
                        rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                        !===========================================================================================================
                        ! P G M S I - this represents the D E P E N D E N T _ F I E L D S (solid, fluid)
                        !Evaluates the appropriate partial derivative index at position XI for the basis
                        PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV,XI,err,error)
                        rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                        DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                          DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                            !=======================================================================================================
                            ! P G N S I - this represents the L A M B D A
                            PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                            colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                            interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                              & PGNSI*PGMSI*rwg*matrixCoefficient
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivativeIdx
                    ENDDO !localLineNodeIdx

                  CASE(2) !2D interface (face)

                    SELECT CASE(coupledMeshBasis%NUMBER_OF_XI)

                    CASE(2) !Coupled Mesh has 2 xi directions
                      DO localElementNode=1,coupledMeshBasis%NUMBER_OF_NODES
                        DO derivative=1,coupledMeshBasis%NUMBER_OF_DERIVATIVES(localElementNode)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          !=========================================================================================================
                          ! P G M S I
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              !=====================================================================================================
                              ! P G N S I
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivative
                      ENDDO !localElementNode

                    CASE(3) !Coupled Mesh has 3 xi directions
                      connectedFace = elementConnectivity%CONNECTED_FACE
                      decompositionFaceNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(connectedFace)
                      coupledMeshDomainFace=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                        & FACES%FACES(decompositionFaceNumber)
                      DO localFaceNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                        localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                        DO derivativeIdx=1,coupledMeshDomainFace%BASIS%NUMBER_OF_DERIVATIVES(localFaceNodeIdx)
                          derivative=coupledMeshBasis% &
                            & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(derivativeIdx,localFaceNodeIdx,connectedFace)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          !=========================================================================================================
                          ! P G M S I
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              !=====================================================================================================
                              ! P G N S I
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivativeIdx
                      ENDDO !FaceNodeIdx

                    END SELECT !coupledMeshBasis%NUMBER_OF_XI

                  END SELECT !interfaceDependentBasis%NUMBER_OF_XI

                ENDDO !rowComponentIdx
              ENDIF
            ENDDO !GaussPoint

            !Scale factor adjustment
            !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix contribution to the residual for non-linear problems
            !\todo update looping of variables/components for non-zero matrix elements as done above 
            IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
              & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
              CALL FLAG_ERROR("Not implemented.",err,error,*999)
            ELSE
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(elementNumber, &
                  & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% & 
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx) * &
                        & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)% &
                        & INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR%SCALE_FACTORS(colParameterIdx,colComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
              !Scale factor adjustment for the row dependent variable
              IF(coupledMeshDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(coupledMeshElementNumber, &
                  & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                  & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR,err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,interfaceMatrixVariable%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)* &
                        & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                        & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR% &
                        & SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! coupledMeshIdx
        
      CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Interface condition integration type "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType, &
          & "*",err,error))// " is not valid."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT !interfaceCondition%integrationType

    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FLAG_ERROR("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Interface condition method "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%METHOD,"*",err,error))// &
        & " is not valid."
      CALL FLAG_ERROR(localError,err,error,*999)
    END SELECT

    CALL EXITS("SolidFluidOperator_FiniteElementCalculate")
    RETURN
999 CALL ERRORS("SolidFluidOperator_FiniteElementCalculate",err,error)
    CALL EXITS("SolidFluidOperator_FiniteElementCalculate")
    RETURN 1
  
  END SUBROUTINE SolidFluidOperator_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  FUNCTION INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
  
    !Argument variables
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: elementConnectivity !<A pointer to the element connectivity
    TYPE(BASIS_TYPE), POINTER :: interfaceConnectivityBasis !<A pointer to the interface mesh connectivity basis
    INTEGER(INTG), INTENT(IN) :: GaussPoint !< Index to the gauss point which needs to be transformed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(SIZE(elementConnectivity%XI,1))
    !Local Variables
    INTEGER(INTG) :: rowParameterIdx

    CALL ENTERS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM",err,error,*999)
    
    INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM=0.0_DP
    DO rowParameterIdx = 1,interfaceConnectivityBasis%NUMBER_OF_ELEMENT_PARAMETERS
      INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:)= INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:) + &
        & interfaceConnectivityBasis%QUADRATURE%QUADRATURE_SCHEME_MAP &
        & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_BASIS_FNS(rowParameterIdx,NO_PART_DERIV,GaussPoint) * &
        & elementConnectivity%XI(:,1,rowParameterIdx)
    ENDDO
     
    CALL EXITS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM")
    RETURN
999 CALL ERRORS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM",err,error)
    CALL EXITS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM")
    RETURN
    
  END FUNCTION INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM

END MODULE INTERFACE_OPERATORS_ROUTINES
