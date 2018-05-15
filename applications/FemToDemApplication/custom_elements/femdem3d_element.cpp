#include "includes/define.h"
#include <string>
#include "includes/constitutive_law.h"
#include "custom_constitutive/zarate_law.hpp"
#include "femdem3d_element.hpp"
#include "includes/element.h"
#include "includes/node.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "solid_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{
	//***********************DEFAULT CONSTRUCTOR******************************************
	//************************************************************************************

	FemDem3DElement::FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: SmallDisplacementElement(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}
	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	FemDem3DElement::FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		:SmallDisplacementElement(NewId, pGeometry, pProperties)
	{
		//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	FemDem3DElement::FemDem3DElement(FemDem3DElement const& rOther)
		:SmallDisplacementElement(rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
		//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
	}

	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	FemDem3DElement&  FemDem3DElement::operator=(FemDem3DElement const& rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

		SmallDisplacementElement::operator=(rOther);
		return *this;
	}

	//*********************************OPERATIONS*****************************************
	//************************************************************************************

	Element::Pointer FemDem3DElement::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
	{
		//NEEDED TO CREATE AN ELEMENT   
		return Element::Pointer(new FemDem3DElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
	}


	//************************************CLONE*******************************************
	//************************************************************************************

	Element::Pointer FemDem3DElement::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
	{

		//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
		//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

		FemDem3DElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

		return Element::Pointer(new FemDem3DElement(NewElement));
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	FemDem3DElement::~FemDem3DElement()
	{
	}

	void FemDem3DElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		if (this->GetIteration() == 0)
		{
			this->ComputeEdgeNeighbours(rCurrentProcessInfo); 
			this->CalculateLchar(); 
			
			this->IterationPlus();
		}
	}

	void FemDem3DElement::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
	{
		std::vector <std::vector<Element*>> EdgeNeighboursContainer;

		Geometry< Node < 3 > >& NodesCurrentElement = this->GetGeometry();

		Node<3>& pNode0 = NodesCurrentElement[0];
		Node<3>& pNode1 = NodesCurrentElement[1];
		Node<3>& pNode2 = NodesCurrentElement[2];
		Node<3>& pNode3 = NodesCurrentElement[3];
		
		// Neighbour elements of each node of the current element 
		WeakPointerVector< Element >& NeighNode0 = pNode0.GetValue(NEIGHBOUR_ELEMENTS);
		WeakPointerVector< Element >& NeighNode1 = pNode1.GetValue(NEIGHBOUR_ELEMENTS);
		WeakPointerVector< Element >& NeighNode2 = pNode2.GetValue(NEIGHBOUR_ELEMENTS);
		WeakPointerVector< Element >& NeighNode3 = pNode3.GetValue(NEIGHBOUR_ELEMENTS);

		// Nodal neighbours container
		std::vector<WeakPointerVector< Element >> NodalNeighbours;
		NodalNeighbours.push_back(NeighNode0);
		NodalNeighbours.push_back(NeighNode1);
		NodalNeighbours.push_back(NeighNode2);
		NodalNeighbours.push_back(NeighNode3);

		// Aux indexes
		Matrix NodeIndexes = ZeroMatrix(6, 2);
		this->SetNodeIndexes(NodeIndexes);

		// Loop over EDGES to assign the elements that share that edge -> Fill mEdgeNeighboursContainer
		for (int edge = 0 ; edge < 6; edge++)
		{
			int NodeIndex1 = NodeIndexes(edge, 0);
			int NodeIndex2 = NodeIndexes(edge, 1);

			// Neigh elements of local node 1 and 2  // 
			WeakPointerVector< Element >& NeighOfNode1 = NodalNeighbours[NodeIndex1];
			WeakPointerVector< Element >& NeighOfNode2 = NodalNeighbours[NodeIndex2];

			int NodeId1 = NodesCurrentElement[NodeIndex1].Id();
			int NodeId2 = NodesCurrentElement[NodeIndex2].Id();

			std::vector<Element*> EdgeSharedElementsNode1;
			// Loop over neigh elements of the node 1
			for (int neigh_elem = 0; neigh_elem < NeighOfNode1.size();neigh_elem++)
			{
				// Nodes of the neigh element
				Geometry< Node < 3 > >& NodesNeighElem = NeighOfNode1[neigh_elem].GetGeometry();

				// Loop over the nodes of the neigh element
				for (int neigh_elem_node = 0; neigh_elem_node < 4 ; neigh_elem_node++)
				{
					int NeighElementNodeId = NodesNeighElem[neigh_elem_node].Id();

					if (NeighElementNodeId == NodeId2 & this->Id() != NeighOfNode1[neigh_elem].Id())
					{
						EdgeSharedElementsNode1.push_back(&NeighOfNode1[neigh_elem]); // ( [] returns an Element object!!)
					} 
				}
			}

			std::vector<Element*> EdgeSharedElementsNode2;
			// Loop over neigh elements of the node 2
			for (int neigh_elem = 0; neigh_elem < NeighOfNode2.size();neigh_elem++)
			{
				// Nodes of the neigh element
				Geometry< Node < 3 > >& NodesNeighElem = NeighOfNode2[neigh_elem].GetGeometry();

				// Loop over the nodes of the neigh element
				for (int neigh_elem_node = 0; neigh_elem_node < 4 ;neigh_elem_node++)
				{
					int NeighElementNodeId = NodesNeighElem[neigh_elem_node].Id();

					if (NeighElementNodeId == NodeId1 & this->Id() != NeighOfNode2[neigh_elem].Id())
					{
						EdgeSharedElementsNode2.push_back(&NeighOfNode2[neigh_elem]);
					}
				}
			}

			// Let's create the vector of neighbour elements for this edge
			std::vector<Element*> EdgeSharedElements = EdgeSharedElementsNode1;

			// Add the neigh elements from the node 2
			for (int i = 0; i < EdgeSharedElementsNode2.size(); i++)
			{
				int aux = 0;

				for (int j = 0; j < EdgeSharedElements.size(); j++)
				{
					if (EdgeSharedElementsNode2[i]->Id() == EdgeSharedElements[j]->Id()) aux++;
				}

				if (aux == 0) EdgeSharedElements.push_back(EdgeSharedElementsNode2[i]);
			}

			EdgeNeighboursContainer.push_back(EdgeSharedElements);

		} // End loop edges

		// Storages the information inside the element
		this->SaveEdgeNeighboursContainer(EdgeNeighboursContainer);

	} // End finding edge neighbour elements

	void FemDem3DElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		double CurrentfSigma = 0.0, damage_element = 0.0;

		//Loop over edges
		for (int cont = 0; cont < 6; cont++)
		{
			this->Set_Convergeddamages(this->Get_NonConvergeddamages(cont), cont);
			this->SetConverged_f_sigmas(this->Get_NonConvergedf_sigma(cont), cont);
			CurrentfSigma = this->GetConverged_f_sigmas(cont);
			if (CurrentfSigma > this->Get_threshold(cont)) { this->Set_threshold(CurrentfSigma, cont); }

		} // End Loop over edges

		damage_element = this->Get_NonConvergeddamage();
		this->Set_Convergeddamage(damage_element);

		if (damage_element > 0.0) 
		{
			this->SetValue(IS_DAMAGED, 1);
		}
		
		if (damage_element >= 0.98)
		{
			this->Set(ACTIVE, false);  
			double old_threshold = this->GetValue(STRESS_THRESHOLD);
			this->SetValue(INITIAL_THRESHOLD, old_threshold);
		}

		this->ResetNonConvergedVars();
		this->SetValue(DAMAGE_ELEMENT, damage_element);


		// Reset the nodal force flag for the next time step
		Geometry< Node < 3 > >& NodesElement = this->GetGeometry();
		for (int i = 0; i < 3; i++)
		{
			#pragma omp critical
			{
				NodesElement[i].SetValue(NODAL_FORCE_APPLIED, false);
			}
		}
	}

	void FemDem3DElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
	{
		//*****************************
		KRATOS_TRY

		bool is_active = true;
		if (this->IsDefined(ACTIVE))
		{
			is_active = this->Is(ACTIVE);
		}

		// Inactive elements can have negative determinant of the Jacobian
		if (is_active == true)
		{
			//1.-Initialize sizes for the system components:
			const unsigned int number_of_nodes = GetGeometry().size();
			const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

			Vector StrainVector(voigt_size);
			noalias(StrainVector) = ZeroVector(voigt_size);
			Vector StressVector(voigt_size);
			noalias(StressVector) = ZeroVector(voigt_size);
			Matrix ConstitutiveMatrix(voigt_size, voigt_size);
			noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
			Matrix B(voigt_size, dimension*number_of_nodes);
			noalias(B) = ZeroMatrix(voigt_size, dimension*number_of_nodes);
			Matrix DN_DX(number_of_nodes, dimension);
			noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);


			//deffault values for the infinitessimal theory
			double detF = 1;
			Matrix F(dimension, dimension);
			noalias(F) = identity_matrix<double>(dimension);

			//3.-Calculate elemental system:

			//reading integration points
			const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

			//get the shape functions [N] (for the order of the default integration method)
			const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

			//get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
			const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

			//calculate delta position (here coincides with the current displacement)
			Matrix DeltaPosition(number_of_nodes, dimension);
			noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
			DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

			//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
			GeometryType::JacobiansType J;
			J.resize(1, false);
			J[0].resize(dimension, dimension, false);
			noalias(J[0]) = ZeroMatrix(dimension, dimension);
			J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

			// Loop Over Integration Points
			for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
			{
				Matrix InvJ(dimension, dimension);
				noalias(InvJ) = ZeroMatrix(dimension, dimension);
				double detJ = 0;
				MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

				if (detJ < 0)
				{
					this->Set(ACTIVE, false); // element alone inside a crack
					detJ = fabs(detJ);
				}

				if (detJ < 0) KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

				//compute cartesian derivatives for this integration point  [dN/dx_n]
				noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

				//set shape functions for this integration point
				Vector N = row(Ncontainer, PointNumber);

				//b.-compute infinitessimal strain
				this->CalculateInfinitesimalStrain(StrainVector, DN_DX);
				this->SetValue(STRAIN_VECTOR, StrainVector);

				ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

				//set constitutive law variables: (it passes only references to this local variables)
				Values.SetStrainVector(StrainVector);
				Values.SetStressVector(StressVector);
				Values.SetConstitutiveMatrix(ConstitutiveMatrix);
				Values.SetShapeFunctionsDerivatives(DN_DX);
				Values.SetShapeFunctionsValues(N);
				//values to be set:
				Values.SetDeterminantF(detF);
				Values.SetDeformationGradientF(F);

				//set constitutive law flags:
				Flags &ConstitutiveLawOptions = Values.GetOptions();

				//compute stress and constitutive matrix
				ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
				ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

				//CALL THE CONSTITUTIVE LAW (for this integration point)
				//(after calling the constitutive law StressVector and ConstitutiveMatrix are set and can be used)
				mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);  
				this->SetValue(STRESS_VECTOR, Values.GetStressVector());

				this->CalculateDeformationMatrix(B, DN_DX);
				this->SetBMatrix(B);

			}
		}
		
		KRATOS_CATCH("")
	}

	void FemDem3DElement::CalculateLocalSystem(
		MatrixType& rLeftHandSideMatrix, 
		VectorType& rRightHandSideVector, 
		ProcessInfo& rCurrentProcessInfo
	)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
		unsigned int system_size = number_of_nodes * dimension;
		if (rLeftHandSideMatrix.size1() != system_size) rLeftHandSideMatrix.resize(system_size, system_size, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size); 
																 
		if (rRightHandSideVector.size() != system_size) rRightHandSideVector.resize(system_size, false);
		noalias(rRightHandSideVector) = ZeroVector(system_size); 

		Matrix DeltaPosition(number_of_nodes, dimension);
		noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0].resize(dimension, dimension, false);
		noalias(J[0]) = ZeroMatrix(dimension, dimension);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);
		
		for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
		{
			const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
			Vector N = row(Ncontainer, PointNumber);

			double detJ = 0;
			Matrix InvJ(dimension, dimension);
			noalias(InvJ) = ZeroMatrix(dimension, dimension);
			MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

			double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;
			const Matrix& B = this->GetBMatrix();
			Vector IntegratedStressVector = ZeroVector(voigt_size);
			Vector DamagesOnEdges = ZeroVector(6);
			
			// Loop over edges of the element
			for (int edge = 0; edge < 6; edge++)
			{
				std::vector<Element*> EdgeNeighbours = this->GetEdgeNeighbourElements(edge);

				Vector AverageStressVector, AverageStrainVector, IntegratedStressVectorOnEdge;
				this->CalculateAverageStressOnEdge(AverageStressVector, EdgeNeighbours);
				this->CalculateAverageStrainOnEdge(AverageStrainVector, EdgeNeighbours);

				double DamageEdge = 0.0; 
				double Lchar = this->Get_l_char(edge);
				this->IntegrateStressDamageMechanics(IntegratedStressVectorOnEdge, DamageEdge,
					AverageStrainVector, AverageStressVector, edge, Lchar );
				
				this->Set_NonConvergeddamages(DamageEdge, edge);
				DamagesOnEdges[edge] = DamageEdge;

			} // End loop over edges

			double damage_element = this->CalculateElementalDamage(DamagesOnEdges);
			if (damage_element >= 0.999) { damage_element = 0.999; }
			this->Set_NonConvergeddamage(damage_element);
			
			const Vector& StressVector = this->GetValue(STRESS_VECTOR);
			IntegratedStressVector = (1 - damage_element)*StressVector;
			this->SetIntegratedStressVector(IntegratedStressVector);

			Matrix ConstitutiveMatrix = ZeroMatrix(voigt_size, voigt_size);
			double E  = this->GetProperties()[YOUNG_MODULUS];
			double nu = this->GetProperties()[POISSON_RATIO];
			this->CalculateConstitutiveMatrix(ConstitutiveMatrix, E, nu);

			noalias(rLeftHandSideMatrix) += prod(trans(B), IntegrationWeight *(1 - damage_element)* Matrix(prod(ConstitutiveMatrix, B))); // LHS

			Vector VolumeForce = ZeroVector(dimension);
			VolumeForce = this->CalculateVolumeForce(VolumeForce, N);

			// RHS
			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				int index = dimension * i;
				for (unsigned int j = 0; j < dimension; j++)
				{
					rRightHandSideVector[index + j] += IntegrationWeight * N[i] * VolumeForce[j];
				}
			}

			//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
			noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(B), IntegratedStressVector);

			// Add nodal DEM forces
			Vector NodalRHS = ZeroVector(system_size);
			this->AddDEMContactForces(NodalRHS);
		
			// Add nodal contact forces from the DEM
			noalias(rRightHandSideVector) += NodalRHS;

		}
		KRATOS_CATCH("")
		//*****************************
	}

	void FemDem3DElement::AddDEMContactForces(Vector& rNodalRHS)
	{
		#pragma omp critical
		{
			Geometry< Node < 3 > >& NodesElement = this->GetGeometry();

			// Loop Over nodes to apply the DEM contact forces to the FEM
			for (int i = 0; i < 4; i++)
			{
				bool IsDEM = NodesElement[i].GetValue(IS_DEM);
				bool NodalForceApplied = NodesElement[i].GetValue(NODAL_FORCE_APPLIED);

				if (IsDEM == true && NodalForceApplied == false)
				{
					double ForceX = NodesElement[i].GetValue(NODAL_FORCE_X);
					double ForceY = NodesElement[i].GetValue(NODAL_FORCE_Y);
					double ForceZ = NodesElement[i].GetValue(NODAL_FORCE_Z);

					rNodalRHS[3 * i]     += ForceX;
					rNodalRHS[3 * i + 1] += ForceY;
					rNodalRHS[3 * i + 2] += ForceZ;

					NodesElement[i].SetValue(NODAL_FORCE_APPLIED, true);
				}
			}
		}
	}

	void FemDem3DElement::CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

		for (unsigned int i = 0; i < number_of_nodes; i++)
		{
			unsigned int index = 3 * i;

			rB( 0, index + 0 ) = rDN_DX( i, 0 );
			rB( 1, index + 1 ) = rDN_DX( i, 1 );
			rB( 2, index + 2 ) = rDN_DX( i, 2 );

			rB( 3, index + 0 ) = rDN_DX( i, 1 );
			rB( 3, index + 1 ) = rDN_DX( i, 0 );

			rB( 4, index + 1 ) = rDN_DX( i, 2 );
			rB( 4, index + 2 ) = rDN_DX( i, 1 );

			rB( 5, index + 0 ) = rDN_DX( i, 2 );
			rB( 5, index + 2 ) = rDN_DX( i, 0 );
		}	
	}

	void FemDem3DElement::CalculateConstitutiveMatrix(
		Matrix& rConstitutiveMatrix, 
		const double rYoungModulus,
		const double rPoissonCoefficient
	)
	{
		rConstitutiveMatrix.clear();

		// 3D linear elastic constitutive matrix
		rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
		rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );
		rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 );

		rConstitutiveMatrix ( 3 , 3 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-2.0*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
		rConstitutiveMatrix ( 4 , 4 ) = rConstitutiveMatrix ( 3 , 3 );
		rConstitutiveMatrix ( 5 , 5 ) = rConstitutiveMatrix ( 3 , 3 );

		rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
		rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

		rConstitutiveMatrix ( 0 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
		rConstitutiveMatrix ( 2 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

		rConstitutiveMatrix ( 1 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
		rConstitutiveMatrix ( 2 , 1 ) = rConstitutiveMatrix ( 0 , 1 );
	}

	void FemDem3DElement::CalculateDN_DX(Matrix& rDN_DX, int PointNumber)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		//reading integration points
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		//get the shape functions [N] (for the order of the default integration method)
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

		//get the shape functions parent coordinates derivative [dN/d�] (for the order of the default integration method)
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
		//calculate delta position (here coincides with the current displacement)
		Matrix DeltaPosition = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

		//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0] = ZeroMatrix(1, 1);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);
		//a.-compute element kinematics

		//calculating the inverse of the jacobian for this integration point[d�/dx_n]
		Matrix InvJ = ZeroMatrix(dimension, dimension);
		double detJ = 0;
		MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

		if (detJ < 0)
			KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

			//compute cartesian derivatives for this integration point  [dN/dx_n]
			rDN_DX = prod(DN_De[PointNumber], InvJ);
	}

	void FemDem3DElement::CalculateInfinitesimalStrain(Vector& rStrainVector, const Matrix& rDN_DX)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		Matrix H = zero_matrix<double>(dimension); //[dU/dx_n]
		
		for (unsigned int i = 0; i < number_of_nodes; i++)
		{

			array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            H ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            H ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            H ( 0 , 2 ) += Displacement[0]*rDN_DX ( i , 2 );
            H ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            H ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
            H ( 1 , 2 ) += Displacement[1]*rDN_DX ( i , 2 );
            H ( 2 , 0 ) += Displacement[2]*rDN_DX ( i , 0 );
            H ( 2 , 1 ) += Displacement[2]*rDN_DX ( i , 1 );
            H ( 2 , 2 ) += Displacement[2]*rDN_DX ( i , 2 );
		}

		//Infinitesimal Strain Calculation
		if (rStrainVector.size() != 6) rStrainVector.resize(6, false);

        rStrainVector[0] = H( 0, 0 );

        rStrainVector[1] = H( 1, 1 );

        rStrainVector[2] = H( 2, 2 );

        rStrainVector[3] = ( H( 0, 1 ) + H( 1, 0 ) ); // xy

        rStrainVector[4] = ( H( 1, 2 ) + H( 2, 1 ) ); // yz

        rStrainVector[5] = ( H( 0, 2 ) + H( 2, 0 ) ); // xz
		
		KRATOS_CATCH("")
	}

	void FemDem3DElement::CalculateStressVector(
		Vector& rStressVector, 
		const Matrix& rConstitutiveMAtrix, 
		const Vector& rInfinitesimalStrainVector
	)
	{
		noalias(rStressVector) = prod(rConstitutiveMAtrix, rInfinitesimalStrainVector);
	}

	void FemDem3DElement::CalculatePrincipalStresses(Vector& rPrincipalStressVector, const Vector& StressVector)
	{
		rPrincipalStressVector.resize(3);
		double I1, I2, I3, phi, Num, Denom, II1;
		I1 = this->Calculate_I1_Invariant(StressVector);
		I2 = this->Calculate_I2_Invariant(StressVector);
		I3 = this->Calculate_I3_Invariant(StressVector);
		II1 = I1*I1;

		Num = (2.0*II1 - 9.0*I2)*I1 + 27.0*I3;
		Denom = (II1 - 3.0*I2);

		if (Denom != 0.0)
		{
			phi = Num / (2.0*Denom*sqrt(Denom));

			if (abs(phi) > 1.0)
			{
				if (phi > 0.0) phi = 1.0;
				else phi = -1.0;
			}

			double acosphi = acos(phi);
			phi = acosphi / 3.0;

			double aux1 = 0.666666666666667*sqrt(II1 - 3.0*I2);
			double aux2 = I1 / 3.0;

			rPrincipalStressVector[0] = aux2 + aux1*cos(phi);
			rPrincipalStressVector[1] = aux2 + aux1*cos(phi - 2.09439510239);
			rPrincipalStressVector[2] = aux2 + aux1*cos(phi - 4.18879020478);
		}
		else 
		{
			rPrincipalStressVector = ZeroVector(3);
		}
		
	}

	void FemDem3DElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
	}

	void FemDem3DElement::AverageVector(Vector& rAverageVector, const Vector& v, const Vector& w)
	{
		int n = v.size();
		int m = w.size();
		if (n != m) KRATOS_ERROR << "The dimension of the vectors are different or null";
		rAverageVector.resize(n);

		for (int cont = 0;cont < n;cont++)
		{
			rAverageVector[cont] = (v[cont] + w[cont])*0.5;
		}
	}

	void FemDem3DElement::CalculateAverageStressOnEdge(
		Vector& rAverageVector, 
		const std::vector<Element*> VectorOfElems
	)
	{
		Vector CurrentElementStress = this->GetValue(STRESS_VECTOR);
		rAverageVector = CurrentElementStress;
		int counter = 0;

		for (int elem = 0; elem < VectorOfElems.size(); elem++)
		{
			// Only take into account the active elements
			bool is_active = true;
			if (VectorOfElems[elem]->IsDefined(ACTIVE))
			{
				is_active = VectorOfElems[elem]->Is(ACTIVE);
			}

			if (is_active == true)
			{
				rAverageVector += VectorOfElems[elem]->GetValue(STRESS_VECTOR);
				counter++;
			}
		}
		rAverageVector /= (counter + 1);
	}

	void FemDem3DElement::CalculateAverageStrainOnEdge(
		Vector& rAverageVector, 
		const std::vector<Element*> VectorOfElems
	)
	{
		Vector CurrentElementStress = this->GetValue(STRAIN_VECTOR);
		rAverageVector = CurrentElementStress;
		int counter = 0;

		for (int elem = 0; elem < VectorOfElems.size(); elem++)
		{
			// Only take into account the active elements
			bool is_active = true;
			if (VectorOfElems[elem]->IsDefined(ACTIVE))
			{
				is_active = VectorOfElems[elem]->Is(ACTIVE);
			}

			if (is_active == true)
			{
				rAverageVector += VectorOfElems[elem]->GetValue(STRAIN_VECTOR);
				counter++;
			}
		}
		rAverageVector /= (counter + 1);
	}

	// Double values
	void FemDem3DElement::GetValueOnIntegrationPoints(
		const Variable<double>& rVariable, 
		std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		if (rVariable == DAMAGE_ELEMENT || rVariable == IS_DAMAGED || rVariable == STRESS_THRESHOLD)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}	
	}

	// Vector Values
	void FemDem3DElement::GetValueOnIntegrationPoints(
		const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		if (rVariable == STRAIN_VECTOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_VECTOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_VECTOR_INTEGRATED)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}

	// Tensor variables
	void FemDem3DElement::GetValueOnIntegrationPoints( 
		const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo 
	)
	{
		if (rVariable == STRAIN_TENSOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_TENSOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_TENSOR_INTEGRATED)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}

	// DOUBLE VARIABLES
	void FemDem3DElement::CalculateOnIntegrationPoints(
		const Variable<double>& rVariable, 
		std::vector<double>& rOutput, 
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		if (rVariable == DAMAGE_ELEMENT)
		{
			rOutput.resize(1);
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
				//rOutput[PointNumber] = double(this->Get_Convergeddamage());
				rOutput[PointNumber] = double(this->GetValue(DAMAGE_ELEMENT));
			}
		}

		if (rVariable == IS_DAMAGED)
		{
			rOutput.resize(1);
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
				rOutput[PointNumber] = double(this->GetValue(IS_DAMAGED));
			}
		}

		if (rVariable == STRESS_THRESHOLD)
		{
			rOutput.resize(1);
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
				rOutput[PointNumber] = double(this->GetValue(STRESS_THRESHOLD));
			}
		}
	}

	// 	VECTOR VARIABLES
	void FemDem3DElement::CalculateOnIntegrationPoints(
		const Variable<Vector>& rVariable, 
		std::vector<Vector>& rOutput, 
		const ProcessInfo& rCurrentProcessInfo
	)
	{
		KRATOS_TRY

		if (rVariable == STRESS_VECTOR)
		{
			rOutput[0] = this->GetValue(STRESS_VECTOR);
		}

		if (rVariable == STRAIN_VECTOR)
		{
			rOutput[0] = this->GetValue(STRAIN_VECTOR);
		}

		if (rVariable == STRESS_VECTOR_INTEGRATED)
		{
			rOutput[0] = this->GetIntegratedStressVector();
		}
		KRATOS_CATCH("")
	}

	// 	TENSOR VARIABLES
	void FemDem3DElement::CalculateOnIntegrationPoints(
		const Variable<Matrix >& rVariable, 
		std::vector< Matrix >& rOutput, 
		const ProcessInfo& rCurrentProcessInfo
	)
	{
    	const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    	const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

        if ( rOutput[0].size2() != dimension )
            rOutput[0].resize( dimension, dimension, false );

		if (rVariable == STRESS_TENSOR)
		{
			rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(STRESS_VECTOR));
		}

		if (rVariable == STRAIN_TENSOR)
		{
			rOutput[0] =  MathUtils<double>::StrainVectorToTensor(this->GetValue(STRAIN_VECTOR));
		}

		if (rVariable == STRESS_TENSOR_INTEGRATED)
		{
			rOutput[0] =  MathUtils<double>::StressVectorToTensor(this->GetIntegratedStressVector());
		}

	}

	// Fills the array of characteristic lengths of the element
	void FemDem3DElement::CalculateLchar()
	{
		Geometry< Node < 3 > >& NodesElem = this->GetGeometry(); 
		Matrix Indexes;
		this->SetNodeIndexes(Indexes);
		
		for (int edge = 0; edge < 6; edge++)
		{
			double X1 = NodesElem[Indexes(edge,0)].X();
			double X2 = NodesElem[Indexes(edge,1)].X();
			double Y1 = NodesElem[Indexes(edge,0)].Y();
			double Y2 = NodesElem[Indexes(edge,1)].Y();
			double Z1 = NodesElem[Indexes(edge,0)].Z();
			double Z2 = NodesElem[Indexes(edge,1)].Z();

			double lchar = sqrt((X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) + (Z1-Z2)*(Z1-Z2));
			this->Set_l_char(lchar, edge);
		}
	}

	void FemDem3DElement::Get2MaxValues(Vector& MaxValues, double a, double b, double c)
	{
		MaxValues.resize(2);
		Vector V;
		V.resize(3);
		V[0] = a;
		V[1] = b;
		V[2] = c;
		int n = 3, imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}
		MaxValues[0] = V[2];
		MaxValues[1] = V[1];
	}

	void FemDem3DElement::Get2MinValues(Vector& MaxValues, double a, double b, double c)
	{
		MaxValues.resize(2);
		Vector V;
		V.resize(3);
		V[0] = a;
		V[1] = b;
		V[2] = c;
		int n = 3, imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}
		MaxValues[0] = V[1];
		MaxValues[1] = V[0];
	}
	double FemDem3DElement::Calculate_I1_Invariant(Vector StressVector) 
	{
		return  StressVector[0] + StressVector[1] + StressVector[2];
	}

	double FemDem3DElement::Calculate_I2_Invariant(const Vector StressVector)
	{
		return (StressVector[0] + StressVector[2])*StressVector[1] + StressVector[0]*StressVector[2] +
			- StressVector[3]*StressVector[3] - StressVector[4]*StressVector[4] - StressVector[5]*StressVector[5];
	}

	double FemDem3DElement::Calculate_I3_Invariant(const Vector StressVector)
	{
		return (StressVector[1]*StressVector[2] - StressVector[4]*StressVector[4])*StressVector[0] -
			StressVector[1]*StressVector[5]*StressVector[5] - StressVector[2]*StressVector[3]*StressVector[3] +
			2.0*StressVector[3]*StressVector[4]*StressVector[5];
	}

	void FemDem3DElement::CalculateDeviatorVector(Vector& rDeviator, const Vector StressVector, const double I1)
	{
		rDeviator.resize(6);
		rDeviator = StressVector;
		double Pmean = I1 / 3;

		rDeviator[0] -= Pmean;
		rDeviator[1] -= Pmean;
		rDeviator[2] -= Pmean;
	}

	double FemDem3DElement::Calculate_J2_Invariant(const Vector Deviator)
	{
		return  0.5*(Deviator[0]*Deviator[0] + Deviator[1]*Deviator[1] + Deviator[2]*Deviator[2]) +
			(Deviator[3]*Deviator[3] + Deviator[4]*Deviator[4] + Deviator[5]*Deviator[5]);
	}
	double FemDem3DElement::Calculate_J3_Invariant(const Vector Deviator)
	{
		return	 Deviator[0]*(Deviator[1]*Deviator[2] - Deviator[4]*Deviator[4]) +
			Deviator[3]*(-Deviator[3]*Deviator[2] + Deviator[5]*Deviator[4])     +
			Deviator[5]*(Deviator[3]*Deviator[4] - Deviator[5]*Deviator[1]);
	}
	double FemDem3DElement::Calculate_Theta_Angle(double J2, double J3)
	{
		double sint3;
		sint3 = (-3.0*sqrt(3)*J3) / (2 * J2*sqrt(J2));
		if (sint3 < -0.95) { sint3 = -1;}
		if (sint3 > 0.95)  { sint3 = 1; }
		return asin(sint3) / 3;
	}

	void FemDem3DElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			bool ComputeLumpedMassMatrix = false;
		if (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX))
			if (rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == true)
				ComputeLumpedMassMatrix = true;

		if (ComputeLumpedMassMatrix == false) {

			//create local system components
			LocalSystemComponents LocalSystem;

			//calculation flags   
			LocalSystem.CalculationFlags.Set(SmallDisplacementElement::COMPUTE_LHS_MATRIX);

			VectorType RightHandSideVector = Vector();

			//Initialize sizes for the system components:
			this->InitializeSystemMatrices(rMassMatrix, RightHandSideVector, LocalSystem.CalculationFlags);

			//Set Variables to Local system components
			LocalSystem.SetLeftHandSideMatrix(rMassMatrix);
			LocalSystem.SetRightHandSideVector(RightHandSideVector);

			//Calculate elemental system
			CalculateDynamicSystem(LocalSystem, rCurrentProcessInfo);
		}
		else {
			//lumped
			unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			const unsigned int number_of_nodes = GetGeometry().PointsNumber();
			unsigned int MatSize = dimension * number_of_nodes;

			if (rMassMatrix.size1() != MatSize)
				rMassMatrix.resize(MatSize, MatSize, false);

			noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize);

			double TotalMass = 0;
			TotalMass = this->CalculateTotalMass(TotalMass, rCurrentProcessInfo);

			Vector LumpFact(number_of_nodes);
			noalias(LumpFact) = ZeroVector(number_of_nodes);

			LumpFact = GetGeometry().LumpingFactors(LumpFact);

			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				double temp = LumpFact[i] * TotalMass;

				for (unsigned int j = 0; j < dimension; j++)
				{
					unsigned int index = i * dimension + j;
					rMassMatrix(index, index) = temp;
				}
			}

		}

		KRATOS_CATCH("")
	}
	Vector& FemDem3DElement::CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		if (rVolumeForce.size() != dimension)
			rVolumeForce.resize(dimension, false);

		noalias(rVolumeForce) = ZeroVector(dimension);

		for (unsigned int j = 0; j < number_of_nodes; j++)
		{
			if (GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION)) { // it must be checked once at the begining only
				array_1d<double, 3 >& VolumeAcceleration = GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
				for (unsigned int i = 0; i < dimension; i++)
					rVolumeForce[i] += rN[j] * VolumeAcceleration[i];
			}
		}

		rVolumeForce *= GetProperties()[DENSITY];

		return rVolumeForce;

		KRATOS_CATCH("")
	}

	double FemDem3DElement::GetMaxValue(Vector Strain)
	{
		Vector V;
		int n = Strain.size();
		V.resize(n);

		for (int cont = 0;cont < n;cont++)
		{
			V[cont] = Strain[cont];
		}

		int imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}
		return V[n - 1];
	}

	double FemDem3DElement::GetMaxAbsValue(Vector Strain)
	{
		Vector V;
		int n = Strain.size();
		V.resize(n);

		for (int cont = 0;cont < n;cont++)
		{
			V[cont] = abs(Strain[cont]);
		}
		int imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}
		return V[n - 1];
	}

	double FemDem3DElement::GetMinAbsValue(Vector Strain)
	{
		Vector V;
		V.resize(3);
		V[0] = abs(Strain[0]);
		V[1] = abs(Strain[1]);
		V[2] = abs(Strain[2]);
		int n = 3, imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}
		return V[0];
	}


	// ******* DAMAGE MECHANICS YIELD SURFACES AND EXPONENTIAL SOFTENING ********
	void FemDem3DElement::IntegrateStressDamageMechanics(
		Vector& rIntegratedStress, 
		double& rdamage,
		const Vector& StrainVector, 
		const Vector& StressVector, 
		int cont, 
		double l_char)
	{
		std::string yield_surface = this->GetProperties()[YIELD_SURFACE];

		if ( yield_surface     == "ModifiedMohrCoulomb") { this->ModifiedMohrCoulombCriterion(rIntegratedStress, rdamage, StressVector, cont, l_char); }
		else if (yield_surface == "SimoJu") { this->SimoJuCriterion(rIntegratedStress, rdamage, StrainVector, StressVector, cont, l_char); }
		else if (yield_surface == "Rankine") { this->RankineCriterion(rIntegratedStress, rdamage, StressVector, cont, l_char); }
		else if (yield_surface == "DruckerPrager") { this->DruckerPragerCriterion(rIntegratedStress, rdamage, StressVector, cont, l_char); }
		else if (yield_surface == "RankineFragile") { this->RankineFragileLaw(rIntegratedStress, rdamage, StressVector, cont, l_char); }
		else { KRATOS_ERROR << " Yield Surface not defined "; }
	}

	void FemDem3DElement::ModifiedMohrCoulombCriterion(
		Vector& rIntegratedStress, 
		double& rdamage, 
		const Vector& StressVector,
		int cont, 
		double l_char
	)
	{
		rIntegratedStress.resize(6);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159265359 / 180.0; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];

		// Check input variables 
		if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32 deg " << std::endl; }
		if (sigma_c < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
		if (sigma_t < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
		if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa "; }

		double K1, K2, K3, Rmorh, R, alpha_r, c_max, theta, c_threshold;
		R = abs(sigma_c / sigma_t);
		Rmorh = pow(tan((3.14159265359 / 4) + friction_angle / 2), 2);
		alpha_r = R / Rmorh;
		c_max = abs(sigma_c);
		double sinphi = sin(friction_angle);

		double I1, J2, J3;
		I1 = this->Calculate_I1_Invariant(StressVector);
		Vector Deviator = ZeroVector(6);
		this->CalculateDeviatorVector(Deviator, StressVector, I1);
		J2 = this->Calculate_J2_Invariant(Deviator);
		J3 = this->Calculate_J3_Invariant(Deviator);
		K1 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r)*sinphi;
		K2 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r) / sinphi;
		K3 = 0.5*(1 + alpha_r)*sinphi - 0.5*(1 - alpha_r);

		double n = sigma_c / sigma_t;
		double A = 1.00 / (n*n*Gt*E / (l_char * pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f = 0.0, F = 0.0; /// F = f-c = 0 classical definition of yield surface

		// Check Modified Mohr-Coulomb criterion
		if (I1 == 0.0) { f = 0.0; }
		else
		{
			theta = Calculate_Theta_Angle(J2, J3);
			f = (2.00*tan(3.14159265359*0.25 + friction_angle*0.5) / cos(friction_angle))*((I1*K3 / 3) + sqrt(J2)*(K1*cos(theta) - K2*sin(theta)*sinphi / sqrt(3)));
		}

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);
	
		F = f - c_threshold;

		if (F <= 0.0)  // Elastic region --> Damage is constant 
		{
			rdamage = this->Get_Convergeddamages(cont);
		}
		else
		{
			rdamage = 1 - (c_max / f)*exp(A*(1 - f / c_max)); // Exponential softening law
			if (rdamage > 0.99) { rdamage = 0.99; }
		}

		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - rdamage);
	}

	void FemDem3DElement::RankineCriterion(
		Vector& rIntegratedStress, 
		double& damage, 
		const Vector& StressVector, 
		int cont, 
		double l_char
	) 
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStresses(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0, c_max = 0.0, c_threshold = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159265359 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];
		c_max = abs(sigma_t);

		double A = 1.00 / (Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f, F; /// F = f-c = 0 classical definition of yield surface
		f = GetMaxValue(PrincipalStressVector);

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);

		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));  // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	void FemDem3DElement::DruckerPragerCriterion(
		Vector& rIntegratedStress, 
		double& damage, 
		const Vector& StressVector, 
		int cont, 
		double l_char
	)
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStresses(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159265359 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];

		// Check input variables 
		if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32� " << std::endl; }
		if (sigma_c < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
		if (sigma_t < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
		if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa "; }

		double  c_max, c_threshold;
		c_max = abs(sigma_t*(3 + sin(friction_angle)) / (3 * sin(friction_angle) - 3));

		double I1, J2; 
		I1 = Calculate_I1_Invariant(StressVector);
		Vector Deviator = ZeroVector(6);
		this->CalculateDeviatorVector(Deviator, StressVector, I1);
		J2 = Calculate_J2_Invariant(Deviator);

		double A = 1.00 / (Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f, F; /// F = f-c = 0 classical definition of yield surface

		double CFL = 0.0, TEN0 = 0.0;

		// Check DruckerPrager criterion
		if (I1 == 0.0) { f = 0; }
		else
		{
			CFL = -sqrt(3)*(3 - sin(friction_angle)) / (3 * sin(friction_angle) - 3);
			TEN0 = 6 * I1*sin(friction_angle) / (sqrt(3)*(3 - sin(friction_angle))) + sqrt(J2);
			f = abs(CFL*TEN0);
		}

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);
		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	void FemDem3DElement::SimoJuCriterion(
		Vector& rIntegratedStress, 
		double& damage,  
		const Vector& StrainVector,  
		const Vector& StressVector, 
		int cont, 
		double l_char
	)
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStresses(PrincipalStressVector, StressVector);

		double sigma_t = 0.0, E = 0.0, Gt = 0.0, sigma_c = 0.0, n = 0;
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];

		double c_max, c_threshold;
		n = abs(sigma_c / sigma_t);
		c_max = abs(sigma_c) / sqrt(E);

		double SumA = 0.0, SumB = 0.0, SumC = 0.0, ere0 = 0.0, ere1 = 0.0;
		for (int cont = 0;cont < 2;cont++)
		{
			SumA += abs(PrincipalStressVector[cont]);
			SumB += 0.5*(PrincipalStressVector[cont]  + abs(PrincipalStressVector[cont]));
			SumC += 0.5*(-PrincipalStressVector[cont] + abs(PrincipalStressVector[cont]));
		}
		ere0 = SumB / SumA;
		ere1 = SumC / SumA;

		double f = 0, F = 0; /// F = f-c = 0 classical definition of yield surface

		// Check SimoJu criterion
		if (StrainVector[0] == 0 && StrainVector[1] == 0) { f = 0; }
		else
		{
			double auxf = 0.0;
			for (int cont = 0; cont < 5; cont++)
			{
				auxf += StrainVector[cont] * StressVector[cont];  // E*S
			}
			f = sqrt(auxf);
			f *= (ere0*n + ere1);
		}

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);
		F = f - c_threshold;

		double A = 1.00 / (Gt*n*n*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	void FemDem3DElement::RankineFragileLaw(
		Vector& rIntegratedStress, 
		double& damage, 
		const Vector& StressVector, 
		int cont, 
		double l_char
	) 
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStresses(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0, c_max = 0.0, c_threshold = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159265359 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];
		c_max = abs(sigma_t);

		double A = 1.00 / (Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f, F; /// F = f-c = 0 classical definition of yield surface
		f = GetMaxValue(PrincipalStressVector);

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);

		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
		}
		else
		{
			damage = 0.98;  // Fragile  law
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	// Computes the damage of the element considering different fracture modes
	double FemDem3DElement::CalculateElementalDamage(const Vector& EdgeDamages)
	{
		// 7 modes of fracture of the tetrahedron
		Vector DamageModeFracture = ZeroVector(7);
		DamageModeFracture[0] = 0.3333333334*(EdgeDamages[0] + EdgeDamages[1] + EdgeDamages[2]);
		DamageModeFracture[1] = 0.3333333334*(EdgeDamages[0] + EdgeDamages[3] + EdgeDamages[4]);
		DamageModeFracture[2] = 0.3333333334*(EdgeDamages[1] + EdgeDamages[3] + EdgeDamages[5]);
		DamageModeFracture[3] = 0.25*(EdgeDamages[1] + EdgeDamages[2] + EdgeDamages[3] + EdgeDamages[4]);
		DamageModeFracture[4] = 0.25*(EdgeDamages[0] + EdgeDamages[1] + EdgeDamages[4] + EdgeDamages[5]);
		DamageModeFracture[5] = 0.3333333334*(EdgeDamages[2] + EdgeDamages[4] + EdgeDamages[5]);
		DamageModeFracture[6] = 0.25*(EdgeDamages[0] + EdgeDamages[2] + EdgeDamages[3] + EdgeDamages[5]);

		return this->GetMaxValue(DamageModeFracture);
	}

} // Element