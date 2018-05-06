//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <pybind11/pybind11.h>
//#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"

//#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "python/add_mesh_to_python.h"


//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

//hardening laws
#include "custom_constitutive/zarate_law.hpp"

namespace Kratos
{
	

	namespace Python
	{

		void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
		{
			typedef ConstitutiveLaw                  ConstitutiveLawBaseType;

		    class_<ZarateLaw, ZarateLaw::Pointer>(m, "ZarateLaw")
        	.def(init<>());
		}

	}  // namespace Python.
}  // namespace Kratos.