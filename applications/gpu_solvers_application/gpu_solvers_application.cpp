/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


// System includes


// External includes 
#include "cublas.h"


// Project includes
#include "includes/define.h"
#include "gpu_solvers_application.h"


namespace Kratos
{

	KratosGPUSolversApplication::KratosGPUSolversApplication()
	{
		// Initialize cuBlas
		if (cublasInit() != CUBLAS_STATUS_SUCCESS)
			printf("*** Warning! cuBlas failed to initialize! ***\n");
	}
	
	KratosGPUSolversApplication::~KratosGPUSolversApplication()
	{
		// Shutdown cuBlas
		if (cublasShutdown() != CUBLAS_STATUS_SUCCESS)
			printf("*** Warning! cuBlas failed to shutdown! ***\n");
	}

	void KratosGPUSolversApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosGPUSolversApplication... " << std::endl;

	}


}  // namespace Kratos.


