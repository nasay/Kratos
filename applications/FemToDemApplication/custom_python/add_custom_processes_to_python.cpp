//   
//   Project Name:           
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $

// External includes
#include <pybind11/pybind11.h>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

//#include "processes/find_elements_neighbours_process.h"
#include "custom_processes/adaptive_mesh_refinement_process.hpp"
#include "custom_processes/mapping_variables_process.hpp"  
#include "custom_processes/stress_to_nodes_process.hpp"  

namespace Kratos
{

	namespace Python
	{
		
		using namespace pybind11;

		void AddCustomProcessesToPython(pybind11::module& m)
		{
			typedef Process                           ProcessBaseType;
			typedef AdaptiveMeshRefinementProcess     AdaptiveMeshRefinementProcessType;
			typedef MappingVariablesProcess           MappingVariablesProcessType;
			typedef StressToNodesProcess              StressToNodesProcessType;


			// class_<FindElementalNeighboursProcess, Process>(m,"FindElementalNeighboursProcess")
			// 	.def(init<ModelPart&, int, unsigned int>())
			// 	.def("Execute", &FindElementalNeighboursProcess::Execute)
			// 	;

			// Adaptive Mesh Refinement Process
			class_<AdaptiveMeshRefinementProcessType, Process>(m,"AdaptiveMeshRefinementProcess")
				.def(init < ModelPart&, std::string, std::string, std::string, std::string, double, int >())
				.def("Execute", &AdaptiveMeshRefinementProcessType::Execute)
				;
	
				
			// Mapping Variables Process
			class_<MappingVariablesProcessType, Process>(m,"MappingVariablesProcess")
				.def(init < ModelPart&,ModelPart&, std::string, std::string >())
				.def("Execute", &MappingVariablesProcessType::Execute)
				;

			// Stress extrapolation to Nodes
			class_<StressToNodesProcessType, Process>(m,"StressToNodesProcess")
				.def(init < ModelPart&, unsigned int >())
				.def("Execute", &StressToNodesProcessType::Execute)
				;

		}

	}  // namespace Python.

} // Namespace Kratos