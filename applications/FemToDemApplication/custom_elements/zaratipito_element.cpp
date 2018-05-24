//
//   Project Name:        KratosFemToDemApplication  $
//   Created by:          $Author: Alejandro Cornejo $
//   Last modified by:    $Co-Author: The Same Freak $
//   Date:                $Date:           June 2017 $
//   Revision:            $Revision:             n   $
//
//

#include "includes/define.h"
#include <string>
#include "includes/constitutive_law.h"
#include "custom_constitutive/zarate_law.hpp"
#include "zaratipito_element.hpp"
#include "includes/element.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"

namespace Kratos
{
	//***********************DEFAULT CONSTRUCTOR******************************************
	//************************************************************************************

	ZaratipitoElement::ZaratipitoElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: SmallDisplacementElement(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}
	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	ZaratipitoElement::ZaratipitoElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		:SmallDisplacementElement(NewId, pGeometry, pProperties)
	{
		//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	ZaratipitoElement::ZaratipitoElement(ZaratipitoElement const& rOther)
		:SmallDisplacementElement(rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
		//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
	}

	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	ZaratipitoElement&  ZaratipitoElement::operator=(ZaratipitoElement const& rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

		SmallDisplacementElement::operator=(rOther);
		return *this;
	}

	//*********************************OPERATIONS*****************************************
	//************************************************************************************

	Element::Pointer ZaratipitoElement::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
	{
		//NEEDED TO CREATE AN ELEMENT   
		return Element::Pointer(new ZaratipitoElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
	}


	//************************************CLONE*******************************************
	//************************************************************************************

	Element::Pointer ZaratipitoElement::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
	{

		//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
		//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

		ZaratipitoElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

		return Element::Pointer(new ZaratipitoElement(NewElement));
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	ZaratipitoElement::~ZaratipitoElement()
	{
	}

	//**************************** METHODS ALEJANDRO CORNEJO *****************************
	//************************************************************************************

	void ZaratipitoElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
	{ 
		
	}
	//************************************************************************************
	//************************************************************************************

	void ZaratipitoElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		double CurrentfSigma = 0;

		// Loop over edges
		for (int cont = 0;cont < 3;cont++)  
		{
			this->Set_Convergeddamages(this->Get_NonConvergeddamages(cont), cont);
			this->SetConverged_f_sigma(this->Get_NonConvergedf_sigma(cont), cont);
			CurrentfSigma   = this->GetConverged_f_sigma(cont);
			if (CurrentfSigma > this->Get_threshold(cont)) { this->Set_threshold(CurrentfSigma, cont);}
			
		} // Loop over edges

		double damage_element = this->Get_NonConvergeddamage();
		this->Set_Convergeddamage(damage_element);
		this->SetValue(DAMAGE_ELEMENT, damage_element);

		if (damage_element >= 0.98) 
		{ 
			this->Set(TO_ERASE, true);
			std::cout << "ELIMINADO EL ELEMENTO  " << this->Id() << std::endl;
		}

		// Calculation of the Integrated Stress Vector
		Vector IntegratedStressVector = ZeroVector(3);
		Vector StressVector = this->GetValue(STRESS_VECTOR);
		IntegratedStressVector = CalculateIntegratedStressVector(StressVector, damage_element); 
		this->SetValue(STRESS_VECTOR_INTEGRATED, IntegratedStressVector);  

		// Set to Zero the NonConvergedVars
		this->ResetNonConvergedVars();
	}

	//************************************************************************************
	//************************************************************************************
	void ZaratipitoElement::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) / 2;

		// Initialize Parameters
		Vector StrainVector       = ZeroVector(voigt_size);
		Vector StressVector       = ZeroVector(voigt_size);
		Matrix ConstitutiveMatrix = ZeroMatrix(voigt_size, voigt_size);
		Matrix DN_DX              = ZeroMatrix(number_of_nodes, dimension);

		double E  = this-> GetProperties()[YOUNG_MODULUS];
		double nu = this-> GetProperties()[POISSON_RATIO];

		this->CalculateConstitutiveMatrix(ConstitutiveMatrix, E, nu);  // Compute C0
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
		{
			this->CalculateDN_DX(DN_DX, PointNumber);									    // Compute DN_DX
			this->CalculateInfinitesimalStrain(StrainVector, DN_DX);                        // Compute Strain E 
			this->CalculateStressVector(StressVector, ConstitutiveMatrix, StrainVector);    // Compute Stress  S=C0*E
			this->SetValue(STRESS_VECTOR, StressVector);                                    // Storage of the Stress Vector
		}
		std::string yield_criterion = this->GetProperties()[YIELD_SURFACE];
		if (yield_criterion == "SimoJu" || this->GetProperties()[TANGENT_CONSTITUTIVE_TENSOR] == 1)
			{ 
				this->SetValue(STRAIN_VECTOR, StrainVector);                                // Storage of the Strain Vector
			}
	} // InitializeNonLinearIteration
	 //************************************************************************************
	 //************************************************************************************
	Vector ZaratipitoElement::CalculatePrincipalStress(const Vector StressVector)
	{
		Vector res;
		res.resize(2);

		double sigma1 = 0.5*(StressVector[0] + StressVector[1]) + sqrt(pow(0.5*(StressVector[0] - StressVector[1]), 2) + pow(StressVector[2], 2));
		double sigma2 = 0.5*(StressVector[0] + StressVector[1]) - sqrt(pow(0.5*(StressVector[0] - StressVector[1]), 2) + pow(StressVector[2], 2));
		
		res[0] = sigma1;
		res[1] = sigma2;
		return res;
	}
	//************************************************************************************
	//************************************************************************************
	void ZaratipitoElement::CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,const double &rYoungModulus,const double &rPoissonCoefficient)
	{
		rConstitutiveMatrix.clear();

		if (this->GetProperties()[THICKNESS] == 1) 
		{
			// Plane strain constitutive matrix
			rConstitutiveMatrix(0, 0) = (rYoungModulus*(1.0 - rPoissonCoefficient) / ((1.0 + rPoissonCoefficient)*(1.0 - 2.0*rPoissonCoefficient)));
			rConstitutiveMatrix(1, 1) = rConstitutiveMatrix(0, 0);
			rConstitutiveMatrix(2, 2) = rConstitutiveMatrix(0, 0)*(1.0 - 2.0*rPoissonCoefficient) / (2.0*(1.0 - rPoissonCoefficient));
			rConstitutiveMatrix(0, 1) = rConstitutiveMatrix(0, 0)*rPoissonCoefficient / (1.0 - rPoissonCoefficient);
			rConstitutiveMatrix(1, 0) = rConstitutiveMatrix(0, 1);
		}
		else
		{
			// Plane stress constitutive matrix
			rConstitutiveMatrix(0, 0) = (rYoungModulus) / (1.0 - rPoissonCoefficient*rPoissonCoefficient);
			rConstitutiveMatrix(1, 1) = rConstitutiveMatrix(0, 0);
			rConstitutiveMatrix(2, 2) = rConstitutiveMatrix(0, 0)*(1.0 - rPoissonCoefficient)*0.5;
			rConstitutiveMatrix(0, 1) = rConstitutiveMatrix(0, 0)*rPoissonCoefficient;
			rConstitutiveMatrix(1, 0) = rConstitutiveMatrix(0, 1);
		}
	}
	//************************************************************************************
	//************************************************************************************
	void ZaratipitoElement::CalculateInfinitesimalStrain(Vector& rStrainVector, const Matrix& rDN_DX)
	{
		KRATOS_TRY

			const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		Matrix H = zero_matrix<double>(dimension); //[dU/dx_n]

		if (dimension == 2)
		{

			for (unsigned int i = 0; i < number_of_nodes; i++)
			{

				array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

				H(0, 0) += Displacement[0] * rDN_DX(i, 0);
				H(0, 1) += Displacement[0] * rDN_DX(i, 1);
				H(1, 0) += Displacement[1] * rDN_DX(i, 0);
				H(1, 1) += Displacement[1] * rDN_DX(i, 1);
			}


			//Infinitesimal Strain Calculation
			if (rStrainVector.size() != 3) rStrainVector.resize(3, false);

			rStrainVector[0] = H(0, 0);
			rStrainVector[1] = H(1, 1);
			rStrainVector[2] = (H(0, 1) + H(1, 0)); // xy

		}
		else if (dimension == 3)
		{

			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

				H(0, 0) += Displacement[0] * rDN_DX(i, 0);
				H(0, 1) += Displacement[0] * rDN_DX(i, 1);
				H(0, 2) += Displacement[0] * rDN_DX(i, 2);

				H(1, 0) += Displacement[1] * rDN_DX(i, 0);
				H(1, 1) += Displacement[1] * rDN_DX(i, 1);
				H(1, 2) += Displacement[1] * rDN_DX(i, 2);

				H(2, 0) += Displacement[2] * rDN_DX(i, 0);
				H(2, 1) += Displacement[2] * rDN_DX(i, 1);
				H(2, 2) += Displacement[2] * rDN_DX(i, 2);
			}


			//Infinitesimal Strain Calculation
			if (rStrainVector.size() != 6) rStrainVector.resize(6, false);

			rStrainVector[0] = H(0, 0);
			rStrainVector[1] = H(1, 1);
			rStrainVector[2] = H(2, 2);
			rStrainVector[3] = (H(0, 1) + H(1, 0)); // xy
			rStrainVector[4] = (H(1, 2) + H(2, 1)); // yz
			rStrainVector[5] = (H(0, 2) + H(2, 0)); // zx
		}
		else
		{

			KRATOS_THROW_ERROR(std::invalid_argument, "something is wrong with the dimension", "")

		}


		KRATOS_CATCH("")

	}
	//************************************************************************************
	//************************************************************************************
	void ZaratipitoElement::CalculateStressVector(Vector& rStressVector, Matrix& ConstitutiveMAtrix, Vector& InfinitesimalStrainVector)
	{

		noalias(rStressVector) = prod(ConstitutiveMAtrix, InfinitesimalStrainVector);
	}
	//************************************************************************************
	//************************************************************************************
	void ZaratipitoElement::CalculateDN_DX(Matrix& rDN_DX, int PointNumber)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		//reading integration points
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		//get the shape functions [N] (for the order of the default integration method)
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

		//get the shape functions parent coordinates derivative [dN/d£] (for the order of the default integration method)
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
		//calculate delta position (here coincides with the current displacement)
		Matrix DeltaPosition = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);
		//KRATOS_WATCH(DeltaPosition)
		//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0] = ZeroMatrix(1, 1);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);
			//a.-compute element kinematics

			//calculating the inverse of the jacobian for this integration point[d£/dx_n]
			Matrix InvJ = ZeroMatrix(dimension, dimension);
			double detJ = 0;
			MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

			if (detJ < 0)
				KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

				//compute cartesian derivatives for this integration point  [dN/dx_n]
				rDN_DX = prod(DN_De[PointNumber], InvJ);
	}
	//************************************************************************************
	//************************************************************************************
	void ZaratipitoElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{

	} // FinalizeNonLinearIteration

	//************************************************************************************
    //************************************************************************************
	void ZaratipitoElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{

		KRATOS_TRY
		GeneralVariables Variables;
		this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

		//1.-Initialize sizes for the system components:
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		//resizing as needed the LHS
		unsigned int system_size = number_of_nodes * dimension;

		if (rLeftHandSideMatrix.size1() != system_size)
			rLeftHandSideMatrix.resize(system_size, system_size, false);

		noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size); //resetting LHS

		//resizing as needed the RHS
		if (rRightHandSideVector.size() != system_size)
			rRightHandSideVector.resize(system_size, false);

		rRightHandSideVector = ZeroVector(system_size); //resetting RHS

		//2.- Initialize local variables
		unsigned int voigt_size    = dimension * (dimension + 1) * 0.5;

		Vector StrainVector        = ZeroVector(voigt_size);
		Vector StressVector        = ZeroVector(voigt_size);
		Matrix ConstitutiveMatrix  = ZeroMatrix(voigt_size, voigt_size);
		Matrix B                   = ZeroMatrix(voigt_size, dimension*number_of_nodes);
		Matrix DN_DX               = ZeroMatrix(number_of_nodes, dimension);

		//default values for the infinitessimal theory
		double detF = 1;
		Matrix F = identity_matrix<double>(dimension);
		//3.-Calculate elemental system:

		//reading integration points
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		//get the shape functions [N] (for the order of the default integration method)
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

		//get the shape functions parent coodinates derivative [dN/d£] (for the order of the default integration method)
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

		//calculate delta position (here coincides with the current displacement)
		Matrix DeltaPosition = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);
		
		//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0] = ZeroMatrix(1, 1);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);
		
		for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
		{
			//a.-compute element kinematics
			//calculating the inverse of the jacobian for this integration point[d£/dx_n]
			Matrix InvJ = ZeroMatrix(dimension, dimension);
			double detJ = 0;
			MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);
			if (detJ<0)
				KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

			//compute cartesian derivatives for this integration point  [dN/dx_n]
			DN_DX = prod(DN_De[PointNumber], InvJ);
			//set shape functions for this integration point
			Vector N = row(Ncontainer, PointNumber);
			//b.-compute infinitessimal strain
			this->CalculateInfinitesimalStrain(StrainVector, DN_DX);

			//create constitutive law parameters structure:
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
			// see consitutive_law.h
			mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

			//calculating weights for integration on the "reference configuration"
			double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;

			if (dimension == 2) IntegrationWeight *= GetProperties()[THICKNESS];

			//compute the deformation matrix B
			this->CalculateDeformationMatrix(B, DN_DX);

			// Find Neighbour Elements
			WeakPointerVector< Element >& elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
			if (elem_neigb.size() == 0) { KRATOS_THROW_ERROR(std::invalid_argument, " Neighbour Elements not calculated --> size = ", elem_neigb.size()) }

			// Compute damage on each edge of the element
			double damage[3] = { 0,0,0 };

			// Select one yield criterion
			std::string yield_criterion = this->GetProperties()[YIELD_SURFACE];

			Matrix TangentTensor, AuxTensor;
			AuxTensor = ZeroMatrix(3, 3);
			// Loop over edges
			for (int cont = 0; cont < 3; cont++)
			{
				if (yield_criterion == "ModifiedMohrCoulomb")
				{
					damage[cont] = ModifiedMohrCoulombCriterionCalculateDamageOnEdge(*this, elem_neigb[cont], cont);  // Smoothed stress field

					if (this->GetProperties()[TANGENT_CONSTITUTIVE_TENSOR] == 1 && damage[cont] != 0.0 && StrainVector[0] != 0.0)
					{
						TangentTensor.resize(3, 3);
						CalculateTangentTensorModifiedMohrCoulombCriterion(*this, elem_neigb[cont], cont, TangentTensor, damage[cont]);
						//AuxTensor += TangentTensor; // Adds all the Tangent --> to improve
						/*KRATOS_WATCH(this->Id())
						KRATOS_WATCH(damage[cont])*/
						AuxTensor += TangentTensor;
					}

				}
				else if (yield_criterion == "SimoJu")
				{
					damage[cont] = SimoJuCriterionCalculateDamageOnEdge(*this, elem_neigb[cont], ConstitutiveMatrix, cont);  // Smoothed stress field
				}
				else if (yield_criterion == "DruckerPrager")
				{
					damage[cont] = DruckerPragerCriterionCalculateDamageOnEdge(*this, elem_neigb[cont], cont);
				}
				else if (yield_criterion == "Rankine")
				{
					damage[cont] = RankineCriterionCalculateDamageOnEdge(*this, elem_neigb[cont], cont);
				}
				else { KRATOS_ERROR << " ERROR: Yield surface not defined in the .mdpa -> Fill with ModifiedMohrCoulomb / Rankine / SimoJu or DruckerPrager "; }
			}

			// Damage of the element
			Vector TwoMaxDamages;
			TwoMaxDamages.resize(2);
			TwoMaxDamages = this->Get2MaxValues(damage[0], damage[1], damage[2]); 
			double damage_element = (TwoMaxDamages[0] + TwoMaxDamages[1])*0.5;
			if (damage_element >= 0.999) { damage_element = 0.999; }
			this->Set_NonConvergeddamage(damage_element);	
			
			if (this->GetProperties()[TANGENT_CONSTITUTIVE_TENSOR] == 0 || damage_element == 0.0 || StrainVector[0] == 0.0)
			{
				//compute and add stiffness matrix (LHS = rLeftHandSideMatrix = B*(1 - d)*C0*B
				noalias(rLeftHandSideMatrix) += prod(trans(B), IntegrationWeight * Matrix(prod((1 - damage_element)*ConstitutiveMatrix, B)));
			}
			else 
			{ 
				TangentTensor = AuxTensor / 3.0;
				noalias(rLeftHandSideMatrix) += prod(trans(B), IntegrationWeight * Matrix(prod(TangentTensor, B))); 
				KRATOS_WATCH(TangentTensor)
				KRATOS_WATCH((1 - damage_element)*ConstitutiveMatrix)
				KRATOS_WATCH(damage_element)
				KRATOS_WATCH(this->Id())
			//	KRATOS_WATCH()
				std::cout << "" << std::endl;
			}
			//compute and add external forces 
			Vector VolumeForce = ZeroVector(dimension);
			VolumeForce = this->CalculateVolumeForce(VolumeForce, N);

			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				int index = dimension * i;
				for (unsigned int j = 0; j < dimension; j++)
				{
					rRightHandSideVector[index + j] += IntegrationWeight * N[i] * VolumeForce[j];
				}
			}

			//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
			noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(B), StressVector)*(1 - damage_element);
		}
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	Vector ZaratipitoElement::AverageVector(const Vector a, const Vector b)
	{
		Vector res;
		res.resize(3);

		for (int i = 0; i < a.size(); i++)
			{
				if (a.size() != b.size()) { KRATOS_THROW_ERROR(std::invalid_argument, "Dimension of the vectors mismatch ", 0) }
				else {res[i] = 0.5*(a[i] + b[i]);}
			}

		return res;
	}
	//************************************************************************************
	//************************************************************************************
	void ZaratipitoElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
		std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == DAMAGE_ELEMENT) {
			rValues.resize(1);
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
				rValues[PointNumber] = double(this->GetValue(DAMAGE_ELEMENT));
			}
		}
	}
	void ZaratipitoElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == STRESS_VECTOR) 
			{
					rValues[0] = this->GetValue(STRESS_VECTOR);
			}
		if (rVariable == STRESS_VECTOR_INTEGRATED)
		{
			rValues[0] = this->GetValue(STRESS_VECTOR_INTEGRATED);
		}

		if (rVariable == SMOOTHED_STRESS_VECTOR)
		{
			rValues[0] = this->GetValue(SMOOTHED_STRESS_VECTOR);
		}
	}
	//************************************************************************************
	//************************************************************************************
	Vector ZaratipitoElement::Get2MaxValues(double a, double b, double c)
	{
		Vector res, V;
		res.resize(2);
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

		res[0] = V[2];
		res[1] = V[1];

		return res; // Returns the 2 max values over 3 numbers
	}
	//************************************************************************************
	//************************************************************************************
	double ZaratipitoElement::ModifiedMohrCoulombCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, int cont)  /// Elem1 should be the current one in CalculateLocalSystem
	{
		//Modified Mohr-Coulomb criterion pg.150 S.Oller
		Vector StressVector1 = Elem1.GetValue(STRESS_VECTOR);
		Vector StressVector2 = Elem2.GetValue(STRESS_VECTOR);

		Vector AverageStressVector, PrincipalStressVector;
		AverageStressVector.resize(3);
		PrincipalStressVector.resize(2);

		AverageStressVector   = AverageVector(StressVector1, StressVector2);
		PrincipalStressVector = CalculatePrincipalStress(AverageStressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
		sigma_c        = Elem1.GetProperties()[YIELD_STRESS_C];
		sigma_t        = Elem1.GetProperties()[YIELD_STRESS_T];
		friction_angle = Elem1.GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
	    E              = Elem1.GetProperties()[YOUNG_MODULUS];
		Gt             = Elem1.GetProperties()[FRAC_ENERGY_T];

		// Check input variables 
		if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32º " << std::endl; }
		if (sigma_c < 1e-24){ KRATOS_ERROR  << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
		if (sigma_t < 1e-24) {KRATOS_ERROR  << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
		if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa ";}
		
		double damage = 0.0; 
		
		double K1, K2, K3, Rmorh, R, alpha_r, c_max, theta, c_threshold;
		R       = abs(sigma_c / sigma_t);
		Rmorh   = pow(tan((3.14159 / 4) + friction_angle / 2), 2);
		alpha_r = R / Rmorh;
		c_max   = abs(sigma_c);

		double I1, J2, J3;
		I1 = Calculate_I1_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
		J2 = Calculate_J2_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
		J3 = Calculate_J3_Invariant(PrincipalStressVector[0], PrincipalStressVector[1], I1);
		K1 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r)*sin(friction_angle);
		K2 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r) / sin(friction_angle);
		K3 = 0.5*(1 + alpha_r)*sin(friction_angle) - 0.5*(1 - alpha_r);
		
		double f, F; /// F = f-c = 0 classical definition of yield surface
		
		// Check Modified Mohr-Coulomb criterion
		if (PrincipalStressVector[0] == 0 && PrincipalStressVector[1] == 0) { f = 0; }
		else 
		{
			theta = Calculate_Theta_Angle(J2, J3);
			f = (2.00*tan(3.14159*0.25+friction_angle*0.5) / cos(friction_angle))*((I1*K3 / 3) + sqrt(J2)*(K1*cos(theta) - K2*sin(theta)*sin(friction_angle) / sqrt(3)));
		}

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);
		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamages(cont); 
			this->Set_NonConvergeddamages(damage, cont);
		}  
		else         // F > 0 Increase Damage
		{
			if (this->GetIteration() < 3) // Computes the l_char on each side only once
			{
				Geometry< Node < 3 > >& NodesElem1 = Elem1.GetGeometry();  // 3 nodes of the Element 1
				Geometry< Node < 3 > >& NodesElem2 = Elem2.GetGeometry();  // "                    " 2
				Vector Xcoord, Ycoord;
				Xcoord.resize(3);
				Ycoord.resize(3);

				// Let's find the two shared nodes between the 2 elements 
				int aux = 0;
				double l_char = 0;
				for (int cont = 0; cont < 3; cont++)
				{
					for (int cont2 = 0; cont2 < 3; cont2++)
					{
						if (NodesElem1[cont].Id() == NodesElem2[cont2].Id())
						{
							Xcoord[aux] = NodesElem1[cont].X();
							Ycoord[aux] = NodesElem1[cont].Y();
							aux++;                              // aux > 3 if the two elements are the same one (in fact aux == 9)
						}
					}
				} // End finding nodes

				if (aux < 2) { std::cout << " Something wrong with the elements " << std::endl; }        // Must have at least 2 shared nodes
				double length = 0;

				// Computation of the l_char
				if (aux < 3) {                                                                           // It is not an edge element --> The 2 elements are not equal
					length = pow((pow(Xcoord[0] - Xcoord[1], 2) + pow(Ycoord[0] - Ycoord[1], 2)), 0.5);  // Length of the edge between 2 elements
					l_char = length;                                                                     // Currently the characteristic length is the edge length (can be modified)
				}
				else {  // Edge Element
					double ElementArea = Elem1.GetGeometry().Area();
					l_char = sqrt(4 * ElementArea / sqrt(3));   // Cervera's Formula
				} // l_char computed

				this->Set_l_char(l_char, cont);  // Storages the l_char of this side
				this->IterationPlus();
			}

			double l_char = this->Get_l_char(cont);
			double n = 0;
			double A = 0;
			n = sigma_c / sigma_t;                                  
			A = 1.00 / (n*n*Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);  	// Exponential softening parameter
			if (A < 0){ KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			this->Set_NonConvergeddamages(damage, cont);  
		}
		return damage; 
	} // Damage computed on edge
	//************************************************************************************
	//************************************************************************************
	double ZaratipitoElement::SimoJuCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, Matrix& ConstitutiveMatrix, int cont)
		{
			const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			unsigned int voigt_size = dimension * (dimension + 1) * 0.5;
			//Simo-Ju criterion pg.164 S.Oller
			// Initializing Strain Vectors
			Vector StrainVector1 = Elem1.GetValue(STRAIN_VECTOR);
			Vector StrainVector2 = Elem2.GetValue(STRAIN_VECTOR);
			Vector StressVector1 = Elem1.GetValue(STRESS_VECTOR);
			Vector StressVector2 = Elem2.GetValue(STRESS_VECTOR);

			Vector AverageStrainVector;
			AverageStrainVector.resize(3);
			AverageStrainVector = AverageVector(StrainVector1, StrainVector2);

			Vector AverageStressVector, PrincipalStressVector;
			AverageStressVector.resize(3);
			PrincipalStressVector.resize(2);

			AverageStressVector = AverageVector(StressVector1, StressVector2);
			PrincipalStressVector = CalculatePrincipalStress(AverageStressVector);

			double sigma_t = 0.0, E = 0.0, Gt = 0.0, sigma_c = 0.0, n = 0;
			sigma_t = Elem1.GetProperties()[YIELD_STRESS_T];
			sigma_c = Elem1.GetProperties()[YIELD_STRESS_C];
			E  = Elem1.GetProperties()[YOUNG_MODULUS];
			Gt = Elem1.GetProperties()[FRAC_ENERGY_T];

			// Check input variables 
			// todo
			double damage = 0.0;

			double c_max, c_threshold;
			n = abs(sigma_c / sigma_t);
			c_max = abs(sigma_c) / sqrt(E);

			double SumA = 0.0, SumB = 0.0, SumC = 0.0, ere0 = 0.0, ere1 = 0.0;
			for (int cont = 0;cont < 1;cont++)
				{
					SumA += abs(PrincipalStressVector[cont]);
					SumB += 0.5*(PrincipalStressVector[cont]  + abs(PrincipalStressVector[cont]));
					SumC += 0.5*(-PrincipalStressVector[cont] + abs(PrincipalStressVector[cont]));
				}
			ere0 = SumB / SumA;
			ere1 = SumC / SumA;

			double f = 0, F = 0; /// F = f-c = 0 classical definition of yield surface

			// Check SimoJu criterion
			if (AverageStrainVector[0] == 0 && AverageStrainVector[1] == 0) { f = 0; }
			else
			{
				double auxf = 0.0;
				//Vector auxf = ZeroVector(voigt_size);
				//auxf = prod(trans(AverageStrainVector), ConstitutiveMatrix);
				//f = sqrt(auxf[0] * AverageStrainVector[0] + auxf[1] * AverageStrainVector[1] + auxf[2] * AverageStrainVector[2]);

				for (int cont = 0;cont < 2;cont++)
				{
					auxf += AverageStrainVector[cont] * AverageStressVector[cont];
				}
				f = sqrt(auxf);
				f *= (ere0*n + ere1);
			}

			if (this->Get_threshold(cont) == 0.0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
			c_threshold = this->Get_threshold(cont);
			this->Set_NonConvergedf_sigma(f, cont);
			F = f - c_threshold;

			if (F <= 0)  // Elastic region --> Damage is constant 
			{
				damage = this->Get_Convergeddamages(cont);
				this->Set_NonConvergeddamages(damage, cont);
			}
			else         // F > 0 Increase Damage
			{
				if (this->GetIteration() < 3) // Computes the l_char on each side only once
				{
					Geometry< Node < 3 > >& NodesElem1 = Elem1.GetGeometry();  // 3 nodes of the Element 1
					Geometry< Node < 3 > >& NodesElem2 = Elem2.GetGeometry();  // "                    " 2
					Vector Xcoord, Ycoord;
					Xcoord.resize(3);
					Ycoord.resize(3);

					// Let's find the two shared nodes between the 2 elements 
					int aux = 0;
					double l_char = 0;
					for (int cont = 0; cont < 3; cont++)
					{
						for (int cont2 = 0; cont2 < 3; cont2++)
						{
							if (NodesElem1[cont].Id() == NodesElem2[cont2].Id())
							{
								Xcoord[aux] = NodesElem1[cont].X();
								Ycoord[aux] = NodesElem1[cont].Y();
								aux++;                              // aux > 3 if the two elements are the same one (in fact aux == 9)
							}
						}
					} // End finding nodes

					if (aux < 2) { std::cout << " Something wrong with the elements " << std::endl; }        // Must have at least 2 shared nodes
					double length = 0;

					// Computation of the l_char
					if (aux < 3) {                                                                           // It is not an edge element --> The 2 elements are not equal
						length = pow((pow(Xcoord[0] - Xcoord[1], 2) + pow(Ycoord[0] - Ycoord[1], 2)), 0.5);  // Length of the edge between 2 elements
						l_char = length;                                                                     // Currently the characteristic length is the edge length (can be modified)
					}
					else {  // Edge Element
						double ElementArea = Elem1.GetGeometry().Area();
						l_char = sqrt(4 * ElementArea / sqrt(3));   // Cervera's Formula
					} // l_char computed

					this->Set_l_char(l_char, cont);  // Storages the l_char of this side
					this->IterationPlus();
				}

				double l_char = this->Get_l_char(cont);
				double A = 0;
				A = 1.00 / (Gt*n*n / (l_char *pow(sigma_c, 2)) - 0.5);  	// Exponential softening parameter
				if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }
				damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
				this->Set_NonConvergeddamages(damage, cont);
			}
			return damage;
		}
	//************************************************************************************
	//************************************************************************************
	double ZaratipitoElement::DruckerPragerCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, int cont)
		{
			// Initializing Stress Vectors
			Vector StressVector1 = Elem1.GetValue(STRESS_VECTOR);
			Vector StressVector2 = Elem2.GetValue(STRESS_VECTOR);

			Vector AverageStressVector, PrincipalStressVector;
			AverageStressVector.resize(3);
			PrincipalStressVector.resize(2);

			AverageStressVector = AverageVector(StressVector1, StressVector2);
			PrincipalStressVector = CalculatePrincipalStress(AverageStressVector);

			double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
			sigma_c = Elem1.GetProperties()[YIELD_STRESS_C];
			sigma_t = Elem1.GetProperties()[YIELD_STRESS_T];
			friction_angle = Elem1.GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
			E = Elem1.GetProperties()[YOUNG_MODULUS];
			Gt = Elem1.GetProperties()[FRAC_ENERGY_T];

			// Check input variables 
			if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32º " << std::endl; }
			if (sigma_c < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
			if (sigma_t < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
			if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa "; }

			double damage = 0.0;

			double  c_max,  c_threshold;
			c_max = abs(sigma_t*(3 + sin(friction_angle)) / (3 * sin(friction_angle) - 3));

			double I1, J2;
			I1 = Calculate_I1_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
			J2 = Calculate_J2_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);

			double f, F; /// F = f-c = 0 classical definition of yield surface

			double CFL = 0.0, TEN0 = 0.0;

			// Check DruckerPrager criterion
			if (PrincipalStressVector[0] == 0 && PrincipalStressVector[1] == 0) { f = 0; }
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
				damage = this->Get_Convergeddamages(cont);
				this->Set_NonConvergeddamages(damage, cont);
			}
			else         // F > 0 Increase Damage
			{
				if (this->GetIteration() < 3) // Computes the l_char on each side only once
				{
					Geometry< Node < 3 > >& NodesElem1 = Elem1.GetGeometry();  // 3 nodes of the Element 1
					Geometry< Node < 3 > >& NodesElem2 = Elem2.GetGeometry();  // "                    " 2
					Vector Xcoord, Ycoord;
					Xcoord.resize(3);
					Ycoord.resize(3);

					// Let's find the two shared nodes between the 2 elements 
					int aux = 0;
					double l_char = 0;
					for (int cont = 0; cont < 3; cont++)
					{
						for (int cont2 = 0; cont2 < 3; cont2++)
						{
							if (NodesElem1[cont].Id() == NodesElem2[cont2].Id())
							{
								Xcoord[aux] = NodesElem1[cont].X();
								Ycoord[aux] = NodesElem1[cont].Y();
								aux++;                              // aux > 3 if the two elements are the same one (in fact aux == 9)
							}
						}
					} // End finding nodes

					if (aux < 2) { std::cout << " Something wrong with the elements " << std::endl; }        // Must have at least 2 shared nodes
					double length = 0;

					// Computation of the l_char
					if (aux < 3) {                                                                           // It is not an edge element --> The 2 elements are not equal
						length = pow((pow(Xcoord[0] - Xcoord[1], 2) + pow(Ycoord[0] - Ycoord[1], 2)), 0.5);  // Length of the edge between 2 elements
						l_char = length;                                                                     // Currently the characteristic length is the edge length (can be modified)
					}
					else {  // Edge Element
						double ElementArea = Elem1.GetGeometry().Area();
						l_char = sqrt(4 * ElementArea / sqrt(3));   // Cervera's Formula
					} // l_char computed

					this->Set_l_char(l_char, cont);  // Storages the l_char of this side
					this->IterationPlus();
				}

				double l_char = this->Get_l_char(cont);
				double n = 0;
				double A = 0;
				n = sigma_c / sigma_t;
				A = 1.00 / (n*n*Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);  	// Exponential softening parameter
				if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }
				damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
				this->Set_NonConvergeddamages(damage, cont);
			}
			return damage;
		}
	//************************************************************************************
	//************************************************************************************
	double ZaratipitoElement::RankineCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, int cont)
	{
		Vector StressVector1 = Elem1.GetValue(STRESS_VECTOR);
		Vector StressVector2 = Elem2.GetValue(STRESS_VECTOR);

		Vector AverageStressVector, PrincipalStressVector;
		AverageStressVector.resize(3);
		PrincipalStressVector.resize(2);

		AverageStressVector = AverageVector(StressVector1, StressVector2);
		PrincipalStressVector = CalculatePrincipalStress(AverageStressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0, c_max = 0.0, c_threshold = 0.0;
		sigma_c = Elem1.GetProperties()[YIELD_STRESS_C];
		sigma_t = Elem1.GetProperties()[YIELD_STRESS_T];
		friction_angle = Elem1.GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
		E = Elem1.GetProperties()[YOUNG_MODULUS];
		Gt = Elem1.GetProperties()[FRAC_ENERGY_T];
		c_max = abs(sigma_t);

		double f, F; /// F = f-c = 0 classical definition of yield surface
		f = GetMaxAbsValue(PrincipalStressVector);
		
		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);
		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamages(cont);
			this->Set_NonConvergeddamages(damage, cont);
		}
		else         // F > 0 Increase Damage
		{
			if (this->GetIteration() < 3) // Computes the l_char on each side only once
			{
				Geometry< Node < 3 > >& NodesElem1 = Elem1.GetGeometry();  // 3 nodes of the Element 1
				Geometry< Node < 3 > >& NodesElem2 = Elem2.GetGeometry();  // "                    " 2
				Vector Xcoord, Ycoord;
				Xcoord.resize(3);
				Ycoord.resize(3);

				// Let's find the two shared nodes between the 2 elements 
				int aux = 0;
				double l_char = 0;
				for (int cont = 0; cont < 3; cont++)
				{
					for (int cont2 = 0; cont2 < 3; cont2++)
					{
						if (NodesElem1[cont].Id() == NodesElem2[cont2].Id())
						{
							Xcoord[aux] = NodesElem1[cont].X();
							Ycoord[aux] = NodesElem1[cont].Y();
							aux++;                              // aux > 3 if the two elements are the same one (in fact aux == 9)
						}
					}
				} // End finding nodes

				if (aux < 2) { std::cout << " Something wrong with the elements " << std::endl; }        // Must have at least 2 shared nodes
				double length = 0;

				// Computation of the l_char
				if (aux < 3) {                                                                           // It is not an edge element --> The 2 elements are not equal
					length = pow((pow(Xcoord[0] - Xcoord[1], 2) + pow(Ycoord[0] - Ycoord[1], 2)), 0.5);  // Length of the edge between 2 elements
					l_char = length;                                                                     // Currently the characteristic length is the edge length (can be modified)
				}
				else {  // Edge Element
					double ElementArea = Elem1.GetGeometry().Area();
					l_char = sqrt(4 * ElementArea / sqrt(3));   // Cervera's Formula
				} // l_char computed

				this->Set_l_char(l_char, cont);  // Storages the l_char of this side
				this->IterationPlus();
			}

			double l_char = this->Get_l_char(cont);
			double n = 0;
			double A = 0;
			A = 1.00 / (Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);  	// Exponential softening parameter
			if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			this->Set_NonConvergeddamages(damage, cont);
		}
		return damage;
	}
	//************************************************************************************
	//************************************************************************************
	// Only in 2D
	double ZaratipitoElement::Calculate_I1_Invariant(double sigma1, double sigma2) { return sigma1 + sigma2; }
	double ZaratipitoElement::Calculate_J2_Invariant(double sigma1, double sigma2) 
		{
			return  (pow((sigma1 - sigma2), 2) + pow(sigma1, 2) + pow(sigma2, 2)) / 6;
		}
	double ZaratipitoElement::Calculate_J3_Invariant(double sigma1, double sigma2, double I1) 
		{
			return	(sigma1 - I1 / 3)*((sigma2 - I1 / 3))*(-I1 / 3);
		}
	double ZaratipitoElement::Calculate_Theta_Angle (double J2, double J3)
		{
			double sint3;
			sint3 = (-3.0*sqrt(3)*J3) / (2 * J2*sqrt(J2));
			if (sint3 < -0.95) { sint3 = -1; }
			if (sint3 > 0.95)  { sint3 = 1; }
			return asin(sint3) / 3;
		}

	//"Smoothed" stress field
	void  ZaratipitoElement::AssignSmoothedStress(Element& Elem)
	{
		WeakPointerVector< Element >& elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
		Vector CurrentVec, S1, S2, S3, Av1, Av2, Av3;
		CurrentVec.resize(3);
		S1.resize(3);
		S2.resize(3);
		S3.resize(3);
		Av1.resize(3);
		Av2.resize(3);
		Av3.resize(3);
		CurrentVec = this->GetValue(STRESS_VECTOR);
		S1 = elem_neigb[0].GetValue(STRESS_VECTOR);
		S2 = elem_neigb[1].GetValue(STRESS_VECTOR);
		S3 = elem_neigb[2].GetValue(STRESS_VECTOR);

		Av1 = AverageVector(CurrentVec, S1);
		Av2 = AverageVector(CurrentVec, S2);
		Av3 = AverageVector(CurrentVec, S3);

		Vector SmoothedStressVector;
		SmoothedStressVector.resize(3);
		SmoothedStressVector = (Av1 + Av2 + Av3)*(1.0 / 3.0);
		
		Elem.SetValue(SMOOTHED_STRESS_VECTOR, SmoothedStressVector);
	}

	//************************************************************************************
	//***************************FUNCTIONS FOR DYNAMICS***********************************
	void ZaratipitoElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
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
	Vector& ZaratipitoElement::CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN)
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

	//************************************************************************************
	//***************FUNCTIONS FOR PERTURBATED CONSTITUTIVE TANGENT TENSOR****************
	double ZaratipitoElement::GetMaxAbsValue(Vector Strain)
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
	
		return V[n-1];
	}

	double ZaratipitoElement::GetMinAbsValue(Vector Strain)
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

	void ZaratipitoElement::PerturbateStrainComponent(Vector StrainVector, Vector& PertubatedStrain, double perturbation, int component)
	{
		PertubatedStrain = StrainVector;
		PertubatedStrain[component] += perturbation;
	}

	double ZaratipitoElement::CalculatePerturbation(Vector StrainVector, int component)
	{
		double Pert = 0.0;
		if (StrainVector[component] != 0.0) Pert = (1e-5)*StrainVector[component];
		else Pert = (1e-5)*GetMinAbsValue(StrainVector);
		if (Pert < this->GetMaxAbsValue(StrainVector)*(1e-10)) { Pert = GetMaxAbsValue(StrainVector)*(1e-10); }

		return Pert;
	}

	void ZaratipitoElement::CalculateTangentTensorModifiedMohrCoulombCriterion(Element& Elem1, Element& Elem2, int cont, Matrix& TangentTensorEdge, const double DamageEdge)
	{
		TangentTensorEdge.resize(3, 3);
		Vector StrainVector1 = Elem1.GetValue(STRAIN_VECTOR);
		Vector StrainVector2 = Elem2.GetValue(STRAIN_VECTOR);
		Vector StressVector1 = Elem1.GetValue(STRESS_VECTOR);
		Vector StressVector2 = Elem2.GetValue(STRESS_VECTOR);
		Vector AverageStressVector, IntegratedAverageStress, PerturbedStrain1, PerturbedStrain2, PerturbedStress1, PerturbedStress2,
			DeltaStress, PrincipalPerturbedStress, IntegratedPerturbedStress, AveragePerturbedPredictiveStress;
		PrincipalPerturbedStress.resize(2);
		PerturbedStress1.resize(3);
		PerturbedStress2.resize(3);
		
		double E = Elem1.GetProperties()[YOUNG_MODULUS];
		double nu = Elem1.GetProperties()[POISSON_RATIO];
		Matrix ConstitutiveMatrix = ZeroMatrix(3, 3);
		this->CalculateConstitutiveMatrix(ConstitutiveMatrix, E, nu);
		AverageStressVector.resize(3);
		AverageStressVector = AverageVector(StressVector1, StressVector2);
		IntegratedAverageStress = CalculateIntegratedStressVector(AverageStressVector, DamageEdge); // Integrated Stress on edge


		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, Gt = 0.0;
		sigma_c = Elem1.GetProperties()[YIELD_STRESS_C];
		sigma_t = Elem1.GetProperties()[YIELD_STRESS_T];
		friction_angle = Elem1.GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
		Gt = Elem1.GetProperties()[FRAC_ENERGY_T];
		double damage = 0.0;

		double K1, K2, K3, Rmorh, R, alpha_r, c_max, theta, c_threshold;
		R = abs(sigma_c / sigma_t);
		Rmorh = pow(tan((3.14159 / 4) + friction_angle / 2), 2);
		alpha_r = R / Rmorh;
		c_max = abs(sigma_c);

		double I1 = 0.0, J2 = 0.0, J3 = 0.0;
		K1 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r)*sin(friction_angle);
		K2 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r) / sin(friction_angle);
		K3 = 0.5*(1 + alpha_r)*sin(friction_angle) - 0.5*(1 - alpha_r);
		double Perturb1 = 0.0, Perturb2 = 0.0;

		for (int PertComp = 0;PertComp < 3;PertComp++)  // In each we perturbe one component of the strain tensor
		{
			KRATOS_WATCH(StrainVector1)
			KRATOS_WATCH(StrainVector2)
			Perturb1 = this->CalculatePerturbation(StrainVector1, PertComp);
			Perturb2 = this->CalculatePerturbation(StrainVector2, PertComp);

			/*KRATOS_WATCH(Perturb2)
			KRATOS_WATCH(Perturb1)*/
			if (Perturb1 == 0.0 || Perturb2 == 0.0) { KRATOS_ERROR << " Pertubation == 0 --> The strain tensor is null "; }
			/*if (Perturb1 == 0.0) Perturb1 = 1e-6;
			if (Perturb2 == 0.0) Perturb2 = 1e-6;*/
			this->PerturbateStrainComponent(StrainVector1, PerturbedStrain1, Perturb1, PertComp);       // Perturbed Strain1
			this->PerturbateStrainComponent(StrainVector2, PerturbedStrain2, Perturb2, PertComp);       // Perturbed Strain2
			this->CalculateStressVector(PerturbedStress1, ConstitutiveMatrix, PerturbedStrain1);        // Perturbed Predictive Stress1
			this->CalculateStressVector(PerturbedStress2, ConstitutiveMatrix, PerturbedStrain2);        // Perturbed Predictive Stress2

			AveragePerturbedPredictiveStress.resize(3);
			AveragePerturbedPredictiveStress = AverageVector(PerturbedStress1, PerturbedStress2);
			PrincipalPerturbedStress = CalculatePrincipalStress(AveragePerturbedPredictiveStress);      // Principal Perturbed Predictive Stresses 

			I1 = Calculate_I1_Invariant(PrincipalPerturbedStress[0], PrincipalPerturbedStress[1]);
			J2 = Calculate_J2_Invariant(PrincipalPerturbedStress[0], PrincipalPerturbedStress[1]);
			J3 = Calculate_J3_Invariant(PrincipalPerturbedStress[0], PrincipalPerturbedStress[1], I1);

			double f, F; /// F = f-c = 0 classical definition of yield surface

			// Check Modified Mohr-Coulomb criterion
			if (PrincipalPerturbedStress[0] == 0 && PrincipalPerturbedStress[1] == 0) { f = 0; }
			else
			{
				theta = Calculate_Theta_Angle(J2, J3);
				f = (2.00*tan(3.14159*0.25 + friction_angle*0.5) / cos(friction_angle))*((I1*K3 / 3) + sqrt(J2)*(K1*cos(theta) - K2*sin(theta)*sin(friction_angle) / sqrt(3)));
			}
			c_threshold = this->Get_threshold(cont);
			F = f - c_threshold;

			if (F <= 0)  // Elastic region --> Damage is constant 
			{
				damage = this->Get_Convergeddamages(cont);
			}
			else        
			{
				double l_char = this->Get_l_char(cont);
				double n = 0;
				double A = 0;
				n = sigma_c / sigma_t;
				A = 1.00 / (n*n*Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);  	// Exponential softening parameter
				if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }
				damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));
			}

			//KRATOS_WATCH(DamageEdge)
			//KRATOS_WATCH(damage)

			IntegratedPerturbedStress = CalculateIntegratedStressVector(AveragePerturbedPredictiveStress, damage);
			DeltaStress.resize(3);

			//KRATOS_WATCH(IntegratedPerturbedStress)
			//KRATOS_WATCH(AverageStressVector)

			noalias(DeltaStress) = (IntegratedPerturbedStress - IntegratedAverageStress);
			//noalias(DeltaStress) = (IntegratedPerturbedStress - AverageStressVector);
			/*	std::cout << "----------" << std::endl;
				std::cout << "----------" << std::endl;
				KRATOS_WATCH(StrainVector1)
				KRATOS_WATCH(PerturbedStrain1)
				std::cout << "----------" << std::endl;
				KRATOS_WATCH(IntegratedAverageStress)
				KRATOS_WATCH(IntegratedPerturbedStress)
				std::cout << "----------" << std::endl;
				KRATOS_WATCH(AverageStressVector)
				KRATOS_WATCH(PerturbedStress1)
				KRATOS_WATCH(damage)
				KRATOS_WATCH(DamageEdge)
				KRATOS_WATCH(PertComp)
				KRATOS_WATCH(Perturb1)
				KRATOS_WATCH(DeltaStress/Perturb1)
				std::cout << "----------" << std::endl;
				std::cout << "----------" << std::endl;
				KRATOS_WATCH(DeltaStress);*/

			for (int row = 0;row < 3;row++) // DeltaStress is the i column of the Tangent Tensor
			{
				//KRATOS_WATCH(cont);
				//KRATOS_WATCH(DeltaStress / Perturb);
				//KRATOS_WATCH(TangentTensorEdge)
				TangentTensorEdge(row, PertComp) = DeltaStress[row] / ((Perturb1 + Perturb2)*0.5);
				//KRATOS_WATCH(TangentTensorEdge)
				//	std::cout << "" << std::endl;
			}
		}

		//std::cout << "********************" << std::endl;
		//KRATOS_WATCH(TangentTensorEdge)
		//std::cout << "" << std::endl;
		//KRATOS_WATCH((1-DamageEdge)*ConstitutiveMatrix)
		//KRATOS_WATCH(IntegratedAverageStress)
		//KRATOS_WATCH(IntegratedPerturbedStress)
		//KRATOS_WATCH(DeltaStress / ((Perturb1 + Perturb2)*0.5))
		//	KRATOS_WATCH(DeltaStress)
		//KRATOS_WATCH(damage)
		//KRATOS_WATCH(Perturb1)
		//KRATOS_WATCH(Perturb2)
		//KRATOS_WATCH(DamageEdge)
		//std::cout << "" << std::endl;
	}

} // Namespace Kratos