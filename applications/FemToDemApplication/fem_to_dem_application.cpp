//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Alejandro Cornejo Vel�zquez
//


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"

#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/serializer.h"

#include "fem_to_dem_application.h"
#include "fem_to_dem_application_variables.h"


namespace Kratos {

KratosFemToDemApplication::KratosFemToDemApplication(): KratosApplication("FemToDemApplication"),
//mZaratipitoElement(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mAleCornVelElement(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mFemDem3DElement(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mRomFemDem3DElement(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}




void KratosFemToDemApplication::Register() 
	{


 	// calling base class register to register Kratos components
 	KratosApplication::Register();
	
	//REGISTER VARIABLES FEM2DEM
	KRATOS_REGISTER_VARIABLE(DAMAGE_EDGE1)
	KRATOS_REGISTER_VARIABLE(DAMAGE_EDGE2)
	KRATOS_REGISTER_VARIABLE(DAMAGE_EDGE3)
	KRATOS_REGISTER_VARIABLE(DAMAGE_ELEMENT)
	KRATOS_REGISTER_VARIABLE(STRESS_VECTOR)
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C)
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T)
	KRATOS_REGISTER_VARIABLE(FRAC_ENERGY_T)
	KRATOS_REGISTER_VARIABLE(FRAC_ENERGY_C)
	KRATOS_REGISTER_VARIABLE(ITER)
	KRATOS_REGISTER_VARIABLE(STRESS_VECTOR_INTEGRATED)
	KRATOS_REGISTER_VARIABLE(THRESHOLD)
	KRATOS_REGISTER_VARIABLE(SMOOTHED_STRESS_VECTOR)
	KRATOS_REGISTER_VARIABLE(YIELD_SURFACE)
	KRATOS_REGISTER_VARIABLE(STRAIN_VECTOR)
	KRATOS_REGISTER_VARIABLE(SMOOTHING)
	KRATOS_REGISTER_VARIABLE(IS_DAMAGED)
	KRATOS_REGISTER_VARIABLE(TANGENT_CONSTITUTIVE_TENSOR)
	//KRATOS_REGISTER_VARIABLE(CHARACTERISTIC_LENGTH)
	KRATOS_REGISTER_VARIABLE(MESH_REFINED)
	KRATOS_REGISTER_VARIABLE(IS_DYNAMIC)
	KRATOS_REGISTER_VARIABLE(STRESS_THRESHOLD)
	KRATOS_REGISTER_VARIABLE(INTEGRATION_COEFFICIENT)
	KRATOS_REGISTER_VARIABLE(MAPPING_PROCEDURE)
	KRATOS_REGISTER_VARIABLE(INITIAL_THRESHOLD)
	KRATOS_REGISTER_VARIABLE(IS_DEM)
	KRATOS_REGISTER_VARIABLE(DEM_RADIUS)
	KRATOS_REGISTER_VARIABLE(DEM_GENERATED)
	KRATOS_REGISTER_VARIABLE(INACTIVE_NODE)
	KRATOS_REGISTER_VARIABLE(NUMBER_OF_ACTIVE_ELEMENTS)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_APPLIED)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_X)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_Y)
	KRATOS_REGISTER_VARIABLE(NODAL_FORCE_Z)
	KRATOS_REGISTER_VARIABLE(NODAL_STRESS_VECTOR)
	KRATOS_REGISTER_VARIABLE(EQUIVALENT_NODAL_STRESS)
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EQUIVALENT_NODAL_STRESS_GRADIENT)
	KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUXILIAR_GRADIENT)

	KRATOS_REGISTER_VARIABLE(STRAIN_TENSOR);
	KRATOS_REGISTER_VARIABLE(STRESS_TENSOR);
	KRATOS_REGISTER_VARIABLE(STRESS_TENSOR_INTEGRATED);
	
	// Composite
	KRATOS_REGISTER_VARIABLE(CONCRETE_STRESS_TENSOR);
	KRATOS_REGISTER_VARIABLE(STEEL_STRESS_TENSOR);
	KRATOS_REGISTER_VARIABLE(CONCRETE_STRESS_VECTOR);
	KRATOS_REGISTER_VARIABLE(STEEL_STRESS_VECTOR);
	KRATOS_REGISTER_VARIABLE(YOUNG_MODULUS_STEEL);
	KRATOS_REGISTER_VARIABLE(DENSITY_STEEL);
	KRATOS_REGISTER_VARIABLE(POISSON_RATIO_STEEL);
	KRATOS_REGISTER_VARIABLE(STEEL_VOLUMETRIC_PART);
	KRATOS_REGISTER_VARIABLE(CONCRETE_STRESS_TENSOR_INTEGRATED);
	
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C_STEEL);
	KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T_STEEL);
	KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_STEEL);
	KRATOS_REGISTER_VARIABLE(PLASTIC_DISSIPATION_CAPAP);
	KRATOS_REGISTER_VARIABLE(EQUIVALENT_STRESS_VM);
	
	
	//Register element
	//KRATOS_REGISTER_ELEMENT("ZaratipitoElement", mZaratipitoElement)
	KRATOS_REGISTER_ELEMENT("AleCornVelElement", mAleCornVelElement)
	KRATOS_REGISTER_ELEMENT("FemDem3DElement", mFemDem3DElement)
	KRATOS_REGISTER_ELEMENT("RomFemDem3DElement", mRomFemDem3DElement)
			
	//Register Constitutive Laws
	Serializer::Register("ZarateLaw", mZarateLaw);

	

	}

}  // namespace Kratos.
