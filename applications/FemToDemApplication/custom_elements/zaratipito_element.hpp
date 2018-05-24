//
//   Project Name:        KratosFemToDemApplication $
//   Created by:          $Author:Alejandro Cornejo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                Sept 2016 $
//   Revision:            $Revision:                  0.0 $
//

#if !defined(KRATOS_ZARATIPITO_ELEMENT_H_INCLUDED )
#define  KRATOS_ZARATIPITO_ELEMENT_H_INCLUDED

#include "custom_constitutive/zarate_law.hpp"
#include "includes/constitutive_law.h"

#include "includes/element.h"
//#include "custom_elements/linear_solid_element.hpp"
#include "custom_elements/solid_elements/small_displacement_element.hpp"

namespace Kratos
{

	class ZaratipitoElement : public SmallDisplacementElement    // Derived Element from SolidMechanics
	{

	public:

		/// Default constructors
		ZaratipitoElement(IndexType NewId, GeometryType::Pointer pGeometry);

		ZaratipitoElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

		///Copy constructor
		ZaratipitoElement(ZaratipitoElement const& rOther);

		/// Destructor.
		virtual ~ZaratipitoElement();

		/// Assignment operator.
		ZaratipitoElement& operator=(ZaratipitoElement const& rOther);


		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;


		Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

		ZaratipitoElement()
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

		void CalculateStressVector(Vector& rStressVector, Matrix& rConstitutiveMAtrix, Vector& rInfinitesimalStrainVector);

		Vector CalculatePrincipalStress(const Vector StressVector);

		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		void CalculateLocalSystem
			(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo);

		Vector AverageVector(const Vector v, const Vector w);

		void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		Vector Get2MaxValues(double a, double b, double c);


		// Yield criteria + damage function
		double ModifiedMohrCoulombCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, int cont);
		double SimoJuCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, Matrix& ConstitutiveMatrix, int cont);
		double DruckerPragerCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, int cont);
		double RankineCriterionCalculateDamageOnEdge(Element& Elem1, Element& Elem2, int cont);

		// Stress Invariants in 2D
		double Calculate_I1_Invariant(double sigma1, double sigma2);
		double Calculate_J2_Invariant(double sigma1, double sigma2);
		double Calculate_J3_Invariant(double sigma1, double sigma2, double I1);

		Vector CalculateIntegratedStressVector(Vector StressVector, double damage)
		{
			Vector res = ZeroVector(3);
			res = (1 - damage)*StressVector;
			return res;
		}

		// Lode's angle
		double Calculate_Theta_Angle(double J2, double J3); 

		// Converged values
		void   Set_threshold(double af, int cont) { threshold[cont] = af; }
		double Get_threshold(int cont) { return threshold[cont]; }

		void   Set_Convergeddamage(double af) { damage = af; }
		double Get_Convergeddamage() { return damage; }

		void   SetConverged_f_sigma(double af, int cont) { f_sigma[cont] = af; }
		double GetConverged_f_sigma(int cont) { return f_sigma[cont]; }

		void   Set_Convergeddamages(double af, int cont) { damages[cont] = af; }
		double Get_Convergeddamages(int cont) { return damages[cont]; }

		// Non Converged values
		void   Set_NonConvergeddamages(double af, int cont) { NonConvergeddamages[cont] = af; }
		double Get_NonConvergeddamages(int cont) { return NonConvergeddamages[cont]; }

		void   Set_NonConvergeddamage(double af) { NonConvergeddamage= af; }
		double Get_NonConvergeddamage() { return NonConvergeddamage; }

		void   Set_NonConvergedf_sigma(double af, int cont) { NonConverged_f_sigma[cont] = af; }
		double Get_NonConvergedf_sigma(int cont) { return NonConverged_f_sigma[cont]; }

		void   ResetNonConvergedVars() 
			{	
				this->Set_NonConvergeddamage(0);

				for (int cont = 0;cont < 3;cont++)
				{
					this->Set_NonConvergeddamages(0, cont);
					this->Set_NonConvergedf_sigma(0, cont);
				}
			}


		// Verificacion albertinho
		Vector CalculateFintFromKu(Matrix LHS)
		{

			const unsigned int number_of_nodes = GetGeometry().PointsNumber();
			const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			Vector u;
			u = ZeroVector(6);


				for (unsigned int i = 0; i < number_of_nodes; i++)
				{
					array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

					u[2 * i] = Displacement[0];
					u[2 * i + 1] = Displacement[1];
				}

				return prod(LHS, u);
		}



		// Characteristic length Calculations
		void   Set_l_char(double af, int cont) { l_char[cont] = af; }
		double Get_l_char(int cont) { return l_char[cont]; }

		// Auxiliar functions...
		void IterationPlus() { iteration++; }
		int  GetIteration() { return iteration; }
		void SetToZeroIteration() { iteration = 0; }
		void AssignSmoothedStress(Element& Elem);

		void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
		Vector& CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN);

		// Functions to calculate the Constitutive tangent tensor by numerical derivation
		double GetMaxAbsValue(Vector Strain);
		double GetMinAbsValue(Vector Strain);
		void PerturbateStrainComponent(Vector StrainVector, Vector& PertubatedStrain, double perturbation, int component);
		double CalculatePerturbation(Vector StrainVector, int component);
		void CalculateTangentTensorModifiedMohrCoulombCriterion(Element& Elem1, Element& Elem2, int cont, Matrix& TangentTensorEdge, const double DamageEdge);


	private:
		int iteration = 0;

		// Each component == Each edge
		double f_sigma [3]   = { 0, 0, 0 };     // Mohr-Coulomb equivalent stress
		double threshold [3] = { 0, 0, 0 };     // Stress Threshold
		double damages [3]   = { 0, 0, 0 };     // Converged damage on each edge
		double damage = 0;                        // Converged damage

		double NonConvergeddamages [3]  = { 0, 0, 0 };  // Damages on edges of "i" iteration
		double NonConverged_f_sigma [3] = { 0, 0, 0 };  // Equivalent stress of "i" iteration
		double NonConvergeddamage = 0;                  // Damage of the element of "i" iteration

		double l_char[3] = { 0, 0, 0 };  // Characteristic length on each edge

	};  // Class ZaratipitoElement

}  // Namespace Kratos

#endif // KRATOS_ZARATIPITO_ELEMENT_H_INCLUDED  defined 