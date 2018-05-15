//
//   Project Name:        KratosFemToDemApplication $
//   Created by:          $Author:Alejandro Cornejo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                Sept 2016 $
//   Revision:            $Revision:                  0.0 $
//

#if !defined(KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED )
#define  KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED

#include "custom_constitutive/zarate_law.hpp"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "custom_elements/solid_elements/small_displacement_element.hpp"

namespace Kratos
{


	class AleCornVelElement : public SmallDisplacementElement    // Derived Element from SolidMechanics
	{

	public:

		/// Default constructors
		AleCornVelElement(IndexType NewId, GeometryType::Pointer pGeometry);

		AleCornVelElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

		///Copy constructor
		AleCornVelElement(AleCornVelElement const& rOther);

		/// Destructor.
		virtual ~AleCornVelElement();

		/// Assignment operator.
		AleCornVelElement& operator=(AleCornVelElement const& rOther);


		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;


		Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

		AleCornVelElement()
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

		void CalculatePrincipalStress(Vector& PrincipalStressVector, const Vector StressVector);

		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);
		void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);
		void CalculateLocalSystem
			(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo);

		void AverageVector(Vector& rAverageVector, const Vector& v, const Vector& w);

		void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void Get2MaxValues(Vector& MaxValues, double a, double b, double c);
		void Get2MinValues(Vector& MaxValues, double a, double b, double c);

		void IntegrateStressDamageMechanics(Vector& rIntegratedStress, 
			double& Damage, const Vector StrainVector, const Vector StressVector, int cont, double L_char);


		void ModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);
		void RankineCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);
		void DruckerPragerCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);
		void SimoJuCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StrainVector, const Vector& StressVector, int cont, double L_char);
		void RankineFragileLaw(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);

		void TangentModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& Damage, const Vector& StressVector, int cont, double L_char);

		// Stress Invariants in 2D
		double Calculate_I1_Invariant(double sigma1, double sigma2);
		double Calculate_J2_Invariant(double sigma1, double sigma2);
		double Calculate_J3_Invariant(double sigma1, double sigma2, double I1);

		void CalculateIntegratedStressVector(Vector& rIntegratedStressVector, const Vector rStressVector, const double Damage)
		{
			Vector res = ZeroVector(3);
			res = (1 - Damage) * rStressVector;
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

			for (int cont = 0;cont < 3;cont++)
			{
				this->Set_NonConvergeddamages(0, cont);
				this->Set_NonConvergedf_sigma(0, cont);
			}
		}

		// Characteristic length Calculations
		void   Set_l_char(double af, int cont) { mL_char[cont] = af; }
		double Get_l_char(int cont) { return mL_char[cont]; }
		double CalculateLchar(AleCornVelElement* CurrentElement, const Element& NeibElement, int cont);

		void   SetJ(double af) { mJac = af; }
		double GetJ() { return mJac; }

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
		void PerturbateStrainComponent(const Vector& rStrainVector, Vector& PertubatedStrain, const double perturbation, int component);
		double CalculatePerturbation(const Vector& StrainVector, int component);
		void CalculateTangentTensor(Matrix& rTangentTensor, const Vector& StrainVector, const Vector& IntegratedStressVector ,int cont, double L_char);

		void SetStressVector(Vector toStressVector) { toStressVector.resize(3); mStressVector = toStressVector; }
		Vector GetStressVector() { return mStressVector; }

		void SetStrainVector(Vector toStrainVector) { toStrainVector.resize(3); mStrainVector = toStrainVector; }
		Vector GetStrainVector() { return mStrainVector; }

		void SetIntegratedStressVector(Vector toIntegratedStressVector) { toIntegratedStressVector.resize(3); mIntegratedStressVector = toIntegratedStressVector; }
		Vector GetIntegratedStressVector() { return mIntegratedStressVector; }

		void SetBMatrix(Matrix toBMatrix) { toBMatrix.resize(3, 6); mB = toBMatrix; }
		Matrix GetBMatrix(){ return mB; }

		void CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX);

		void SetIntegrationCoefficient(double tomIntegrationCoefficient){ mIntegrationCoefficient = tomIntegrationCoefficient;}
		double GetIntegrationCoefficient(){ return mIntegrationCoefficient; }

		void AleCornVelElement::SetValueOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void AleCornVelElement::SetValueOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo);
		
	private:
		int iteration = 0;

		// Each component == Each edge
		Vector mF_sigmas = ZeroVector(3);   // Mohr-Coulomb equivalent stress
		Vector mThresholds = ZeroVector(3);   // Stress mThreshold on edge

		double mThreshold = 0.0;
		double mF_sigma   = 0.0;

		Vector mDamages = ZeroVector(3);     // Converged mDamage on each edge
		double mDamage = 0.0;                            // Converged mDamage

		Vector mNonConvergedDamages = ZeroVector(3);    // mDamages on edges of "i" iteration
		Vector mNonConvergedFsigmas  = ZeroVector(3);   // Equivalent stress of "i" iteration

		double mNonConvergedFsigma = 0.0;
		double mNonConvergedDamage = 0.0;       // mDamage of the element of "i" iteration

		Vector mL_char = ZeroVector(3);  // Characteristic length on each edge

		Vector mStressVector = ZeroVector(3);
		Vector mStrainVector = ZeroVector(3);
		Vector mIntegratedStressVector = ZeroVector(3);
		Matrix mB = ZeroMatrix(3, 6);

		double mJac = 0.0;
		double mIntegrationCoefficient = 0.0;

	};// Class AleCornVelElement
	
}// Namespace Kratos
#endif // KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED  defined 