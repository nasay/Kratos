//
//   Project Name:        KratosFemToDemApplication $
//   Created by:          $Author:Alejandro Cornejo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                Sept 2016 $
//   Revision:            $Revision:                  0.0 $
//

#if !defined(KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED )
#define  KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED

#include "custom_constitutive/zarate_law.hpp"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "custom_elements/solid_elements/small_displacement_element.hpp"

namespace Kratos
{


	class FemDem3DElement : public SmallDisplacementElement    // Derived Element from SolidMechanics
	{

	public:

		/// Default constructors
		FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry);

		FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

		///Copy constructor
		FemDem3DElement(FemDem3DElement const& rOther);

		/// Destructor.
		virtual ~FemDem3DElement();

		/// Assignment operator.
		FemDem3DElement& operator=(FemDem3DElement const& rOther);


		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;


		Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

		FemDem3DElement()
		{
		}
		
		// *************** Methods Alejandro Cornejo ***************
		//**********************************************************

		void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
		void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
		void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);
		void CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix, const double &rYoungModulus,
			const double &rPoissonCoefficient);

		void CalculateDN_DX(Matrix& rDN_DX, int PointNumber);

		void CalculateInfinitesimalStrain(Vector& rStrainVector, const Matrix& rDN_DX);

		void CalculateStressVector(Vector& rStressVector, const Matrix& rConstitutiveMAtrix, const Vector& rInfinitesimalStrainVector);

		void CalculatePrincipalStresses(Vector& PrincipalStressVector, const Vector StressVector);

		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);
		void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);
		void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo);

		void CalculateLocalSystem
			(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo);

		void AverageVector(Vector& rAverageVector, const Vector& v, const Vector& w);

		void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
			std::vector<Matrix>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void Get2MaxValues(Vector& MaxValues, double a, double b, double c);
		void Get2MinValues(Vector& MaxValues, double a, double b, double c);

		void IntegrateStressDamageMechanics(Vector& rIntegratedStress, 
			double& rDamage, const Vector StrainVector, const Vector StressVector, int cont, double L_char);


		void ModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);
		void RankineCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);
		void DruckerPragerCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);
		void SimoJuCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StrainVector, const Vector& StressVector, int cont, double L_char);
		void RankineFragileLaw(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);

		void TangentModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);

		// Stress Invariants in 3D
		double Calculate_I1_Invariant(Vector StressVector);
		double Calculate_I2_Invariant(const Vector StressVector);
		double Calculate_I3_Invariant(const Vector StressVector);
		void   CalculateDeviatorVector(Vector& rDeviator, const Vector StressVector, const double I1);
		double Calculate_J2_Invariant(const Vector Deviator);
		double Calculate_J3_Invariant(const Vector Deviator);

		void CalculateIntegratedStressVector(Vector& rIntegratedStressVector, const Vector rStressVector, const double Damage)
		{
			rIntegratedStressVector = (1 - Damage) * rStressVector;
		}

		// Lode's angle
		double Calculate_Theta_Angle(double J2, double J3);

		// Converged values
		void   Set_threshold(double af, int cont) { mThresholds[cont] = af; }
		double Get_threshold(int cont) { return mThresholds[cont]; }

		Vector GetThresholds() { return mThresholds; }
		Vector GetDamages() { return mDamages; }

		void   Set_threshold(double af) { mThreshold = af; }
		double Get_threshold() { return mThreshold; }

		void   Set_Convergeddamage(double af) { mDamage = af; }
		double Get_Convergeddamage() { return mDamage; }

		void   SetConverged_f_sigma(double af) { mF_sigma = af; }
		double GetConverged_f_sigma() { return mF_sigma; }

		void   SetConverged_f_sigmas(double af, int cont) { mF_sigmas[cont] = af; }
		double GetConverged_f_sigmas(int cont) { return mF_sigmas[cont]; }

		void   Set_Convergeddamages(double af, int cont) { mDamages[cont] = af; }
		double Get_Convergeddamages(int cont) { return mDamages[cont]; }

		// Non Converged values
		void   Set_NonConvergeddamages(double af, int cont) { mNonConvergedDamages[cont] = af; }
		double Get_NonConvergeddamages(int cont) { return mNonConvergedDamages[cont]; }

		void   Set_NonConvergeddamage(double af) { mNonConvergedDamage = af; }
		double Get_NonConvergeddamage() { return mNonConvergedDamage; }

		void   Set_NonConvergedf_sigma(double af, int cont) { mNonConvergedFsigmas[cont] = af; }
		double Get_NonConvergedf_sigma(int cont) { return mNonConvergedFsigmas[cont]; }

		void   Set_NonConvergedf_sigma(double af) { mNonConvergedFsigma = af; }
		double Get_NonConvergedf_sigma() { return mNonConvergedFsigma; }

		void   ResetNonConvergedVars()
		{
			this->Set_NonConvergeddamage(0.0);
			this->Set_NonConvergedf_sigma(0.0);

			for (int cont = 0;cont < 6;cont++)
			{
				this->Set_NonConvergeddamages(0, cont);
				this->Set_NonConvergedf_sigma(0, cont);
			}
		}

		// Characteristic length Calculations
		void   Set_l_char(double af, int cont) { mL_char[cont] = af; }
		double Get_l_char(int cont) { return mL_char[cont]; }
		void CalculateLchar();

		//void   SetJ(double af) { mJac = af; }
		//double GetJ() { return mJac; }

		// Auxiliar functions...
		void IterationPlus() { iteration++; }
		int  GetIteration() { return iteration; }
		void SetToZeroIteration() { iteration = 0; }
		//void AssignSmoothedStress(Element& Elem);

		void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
		Vector& CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN);

		// Functions to calculate the Constitutive tangent tensor by numerical derivation
		double GetMaxValue(Vector Strain);
		double GetMaxAbsValue(Vector Strain);
		double GetMinAbsValue(Vector Strain);
		//void PerturbateStrainComponent(const Vector& rStrainVector, Vector& PertubatedStrain, const double& perturbation, int component);
		//double CalculatePerturbation(const Vector& StrainVector, int component);
		//void CalculateTangentTensor(Matrix& rTangentTensor, const Vector& StrainVector, const Vector& IntegratedStressVector ,int cont, double L_char);

		void SetStressVector(Vector toStressVector) { toStressVector.resize(6); mStressVector = toStressVector; }
		Vector GetStressVector() { return mStressVector; }

		void SetStrainVector(Vector toStrainVector) { toStrainVector.resize(6); mStrainVector = toStrainVector; }
		Vector GetStrainVector() { return mStrainVector; }

		void SetIntegratedStressVector(Vector toIntegratedStressVector) { toIntegratedStressVector.resize(6); mIntegratedStressVector = toIntegratedStressVector; }
		Vector GetIntegratedStressVector() { return mIntegratedStressVector; }

		void SetBMatrix(Matrix toBMatrix) {  mB = toBMatrix; }
		Matrix GetBMatrix(){ return mB; }

		void CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX);

		//void SetIntegrationCoefficient(double tomIntegrationCoefficient){ mIntegrationCoefficient = tomIntegrationCoefficient;}
		//double GetIntegrationCoefficient(){ return mIntegrationCoefficient; }

		// Fills mEdgeNeighboursContainer
		void ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo);

		// Storages mEdgeNeighboursContainer
		void SaveEdgeNeighboursContainer(std::vector <std::vector<Element*>> toSave) {mEdgeNeighboursContainer = toSave;}
		std::vector<Element*> GetEdgeNeighbourElements(int edge){return mEdgeNeighboursContainer[edge];}

		void CalculateAverageStressOnEdge(Vector& AverageVector, const std::vector<Element*> VectorOfElems);
		void CalculateAverageStrainOnEdge(Vector& AverageVector, const std::vector<Element*> VectorOfElems);
		void AddDEMContactForces(Vector& NodalRHS);

		void SetNodeIndexes(Matrix& M) // Defines the numbering of the edges with the corresponding nodes
		{
			M.resize(6,2);

			M(0,0) = 0; 
			M(0,1) = 1;
			M(1,0) = 0;
			M(1,1) = 2;
			M(2,0) = 0;
			M(2,1) = 3;
			M(3,0) = 1;
			M(3,1) = 2;
			M(4,0) = 1;
			M(4,1) = 3;
			M(5,0) = 2;
			M(5,1) = 3;
		}

		double CalculateElementalDamage(const Vector& EdgeDamages);
		
	private:
		int iteration = 0;

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

		// Each component == Each edge
		Vector mF_sigmas = ZeroVector(6);     // Equivalent stress
		Vector mThresholds = ZeroVector(6);   // Stress mThreshold on edge

		double mThreshold = 0.0;
		double mF_sigma   = 0.0;

		Vector mDamages = ZeroVector(6);     // Converged mDamage on each edge
		double mDamage = 0.0;                // Converged mDamage

		Vector mNonConvergedDamages = ZeroVector(6);    // mDamages on edges of "i" iteration
		Vector mNonConvergedFsigmas  = ZeroVector(6);   // Equivalent stress of "i" iteration

		double mNonConvergedFsigma = 0.0;
		double mNonConvergedDamage = 0.0;       // mDamage of the element of "i" iteration

		Vector mL_char = ZeroVector(6);  // Characteristic length on each edge

		Vector mStressVector = ZeroVector(voigt_size);
		Vector mStrainVector = ZeroVector(voigt_size);
		Vector mIntegratedStressVector = ZeroVector(voigt_size);
		Matrix mB = ZeroMatrix(voigt_size, dimension*number_of_nodes);

		//double mJac = 0.0;
		//double mIntegrationCoefficient = 0.0;

		// Vector to storage the neigh elements sharing a certain edge
		std::vector <std::vector<Element*>> mEdgeNeighboursContainer;



	}; // Class FemDem3DElement
	
}// Namespace Kratos
#endif // KRATOS_FEMDEM3D_ELEMENT_H_INCLUDED  defined 