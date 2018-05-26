// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Natalia Saiapova
//                   Philipp Bucher
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "geometries/quadrilateral_2d_4.h"
#include "utilities/variable_utils.h"
#include "processes/structured_mesh_generator_process.h"

#include "custom_strategies/custom_convergencecriterias/general_residual_criteria.h"


/*
Things we will test:
- For both residual- & solutionupdate-based => here we can think if it is really necessary for both, maybe for one of them just one test if they have a common baseclass
    - For relative and absolute convergence
        - For different types of variables:
            - No Input (e.g. what the Displacement and the Residual Criterion do atm)
            - One Array3D + Double Variable (what VelPrCriteria does)
            - Two Double Variables
            - Two Array3D Variables (What the criteria in StrucutralMechanics do atm)
            - ... ?
=> makes at least 2*2*(4+x) tests
*/

namespace Kratos
{
    namespace Testing
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;


        typedef typename SparseSpaceType::DataType TDataType;

        typedef typename SparseSpaceType::MatrixType TSystemMatrixType;
        typedef typename SparseSpaceType::VectorType TSystemVectorType;
        typedef typename SparseSpaceType::MatrixPointerType TSystemMatrixPointerType;
        typedef typename SparseSpaceType::VectorPointerType TSystemVectorPointerType;

        typedef typename LocalSpaceType::MatrixType TLocalMatrixType;
        typedef typename LocalSpaceType::VectorType TLocalVectorType;
        typedef typename LocalSpaceType::MatrixPointerType TLocalMatrixPointerType;
        typedef typename LocalSpaceType::VectorPointerType TLocalVectorPointerType;


        typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
        typedef Kratos::unique_ptr<ConvergenceCriteriaType> ConvergenceCriteriaPointerType;

        typedef GeneralResidualCriteria<SparseSpaceType, LocalSpaceType> GenConvergenceCriteriaType;
        typedef Kratos::unique_ptr<GenConvergenceCriteriaType> GenConvergenceCriteriaPointerType;

        typedef ModelPart::DofsArrayType DofsArrayType;

        typedef Node<3> NodeType;
        typedef Kratos::unique_ptr<NodeType> NodeUniquePointerType;

        // This function initializes the test depending on the passed configuration
        // TODO for now the dofs are hard-coded, this should be made variable
        void SetUpTest(ModelPart& rModelPart,
                       DofsArrayType& rDofSet,
                       TSystemVectorType& rSystemVec_Dx,
                       TSystemVectorType& rSystemVec_b,
                       const std::size_t NumDivisions)
        {
            rModelPart.AddNodalSolutionStepVariable(ROTATION);
            rModelPart.AddNodalSolutionStepVariable(VELOCITY);
            rModelPart.AddNodalSolutionStepVariable(PRESSURE);

            Node<3>::Pointer p_point1(new Node<3>(1, 0.00, 0.00, 0.00));
            Node<3>::Pointer p_point2(new Node<3>(2, 0.00, 10.00, 0.00));
            Node<3>::Pointer p_point3(new Node<3>(3, 10.00, 10.00, 0.00));
            Node<3>::Pointer p_point4(new Node<3>(4, 10.00, 0.00, 0.00));

            Quadrilateral2D4<Node<3> > geometry(p_point1, p_point2, p_point3, p_point4);

            Parameters mesher_parameters(R"(
            {
                "number_of_divisions" : 0,
                "element_name"        : "Element2D3N",
                "create_skin_sub_model_part": false
            }  )");

            mesher_parameters["number_of_divisions"].SetInt(NumDivisions);

            StructuredMeshGeneratorProcess(geometry, rModelPart, mesher_parameters).Execute();

            const std::size_t system_size = rModelPart.NumberOfNodes();

            KRATOS_INFO("System Size") << system_size << std::endl;

            for (auto& r_node : rModelPart.Nodes())
            {
                r_node.AddDof(ROTATION_X);
                r_node.AddDof(ROTATION_Y);
                r_node.AddDof(ROTATION_Z);
                r_node.AddDof(VELOCITY_X);
                r_node.AddDof(VELOCITY_Z);
                r_node.AddDof(VELOCITY_Y);
                r_node.AddDof(PRESSURE);
            }

            VariableUtils var_utils();
            // var_utils.AddDof(ROTATION, rModelPart);
            // VariableUtils().AddDof(VELOCITY, rModelPart);
            // VariableUtils().AddDof(PRESSURE, rModelPart);

            const std::size_t num_dofs_per_node = 7; // TODO hard coded for now, to be made variable

            std::vector< Dof<double>::Pointer > DoF;
            DoF.reserve(system_size * num_dofs_per_node);
            for (auto& r_node : rModelPart.Nodes())
            {
                DoF.push_back(r_node.pGetDof(ROTATION_X));
                DoF.push_back(r_node.pGetDof(ROTATION_Y));
                DoF.push_back(r_node.pGetDof(ROTATION_Z));
                DoF.push_back(r_node.pGetDof(VELOCITY_X));
                DoF.push_back(r_node.pGetDof(VELOCITY_Y));
                DoF.push_back(r_node.pGetDof(VELOCITY_Z));
                DoF.push_back(r_node.pGetDof(PRESSURE));
            }

            rDofSet.clear();
            rDofSet.reserve(DoF.size());
            for (auto it= DoF.begin(); it!= DoF.end(); it++)
                rDofSet.push_back( it->get() );

            rDofSet.Sort();

            // Initializing the Solution Vector
            SparseSpaceType::Resize(rSystemVec_Dx, system_size*num_dofs_per_node);
            SparseSpaceType::SetToZero(rSystemVec_Dx);
            SparseSpaceType::Resize(rSystemVec_b, system_size*num_dofs_per_node);
            SparseSpaceType::SetToZero(rSystemVec_b);
        }

        // This function sets the solution on the dofs and the system vectors
        // Usually this would be done by setting up and solving the system
        // Natasha please implement this, think abt sth smart to set an apropriate solution
        void SetSolution(DofsArrayType& rDofSet,
                         TSystemVectorType& rSystemVec_Dx,
                         TSystemVectorType& rSystemVec_b)
        {
            for (auto& r_dof : rDofSet)
                r_dof.GetSolutionStepValue() = 5.123; // TODO so far this is a random value

            const std::size_t system_size = SparseSpaceType::Size(rSystemVec_Dx);
            KRATOS_ERROR_IF(system_size != SparseSpaceType::Size(rSystemVec_b))
                << "Sytsem Vector sizes are inconsistent!" << std::endl;

            for (std::size_t i=0; i< system_size; ++i)
            {
                rSystemVec_Dx[i] = -5.147;
                rSystemVec_b[i] = 565.147;
            }
        }

        KRATOS_TEST_CASE_IN_SUITE(GeneralResudialConvergenceCriteriaTest1, KratosStructuralMechanicsFastSuite)
        {
        KRATOS_INFO("CONVERGENCE_CRITERIA") << "Some regular print" << std::endl;
            DofsArrayType dofs_array;
            const std::size_t num_divisions = 10;

            // Natasha even though I am initializing both b and Dx the implementation of the criteria should
            // for now only consider b!
            // In a second step we will make it more general that it can be used with either b or Dx
            TSystemMatrixPointerType p_dummy_system_matrix = SparseSpaceType::CreateEmptyMatrixPointer();
            TSystemVectorPointerType p_system_vector_b = SparseSpaceType::CreateEmptyVectorPointer();
            TSystemVectorPointerType p_system_vector_Dx = SparseSpaceType::CreateEmptyVectorPointer();

            TSystemMatrixType& r_system_matrix = *p_dummy_system_matrix;
            TSystemVectorType& r_system_vector_b = *p_system_vector_b;
            TSystemVectorType& r_system_vector_Dx = *p_system_vector_Dx;

            ModelPart dummy_model_part("dummy");
            SetUpTest(dummy_model_part, dofs_array, r_system_vector_Dx, r_system_vector_b, num_divisions);

            const TDataType NewRatioTolerance = 1e-7;
            const TDataType AlwaysConvergedNorm = 1e-5;

            ConvergenceCriteriaPointerType p_conv_crit = Kratos::make_unique<GenConvergenceCriteriaType>(
                NewRatioTolerance, AlwaysConvergedNorm);

            const std::size_t num_solution_steps = 5; // Natasha this should be used by "SetSolution"
            // in order to set a solution that will make the problem converge

            p_conv_crit->Initialize(dummy_model_part);

            bool is_converged = false;

            std::size_t counter = 0;

            while(!is_converged)
            {
                SetSolution(dofs_array, r_system_vector_Dx, r_system_vector_b);
                p_conv_crit->InitializeSolutionStep(dummy_model_part,
                                                    dofs_array,
                                                    r_system_matrix,
                                                    r_system_vector_Dx,
                                                    r_system_vector_b);

                is_converged = p_conv_crit->PostCriteria(dummy_model_part,
                                                    dofs_array,
                                                    r_system_matrix,
                                                    r_system_vector_Dx,
                                                    r_system_vector_b);

                counter += 1;

                if (counter > 10) // TODO this is an intermediate solution until the criteria is implemented
                    is_converged = true;
            }
        }

    } // namespace Testing
}  // namespace Kratos.
