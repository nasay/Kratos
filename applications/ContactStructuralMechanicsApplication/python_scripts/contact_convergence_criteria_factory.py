## This script collects the available convergence criteria to be used in the ContactStructuralMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KM backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics as KM

# Check that applications were imported in the main script
KM.CheckRegisteredApplications("StructuralMechanicsApplication")
KM.CheckRegisteredApplications("ContactStructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        self.echo_level = convergence_criterion_parameters["echo_level"].GetInt()
        self.convergence_criterion_name = convergence_criterion_parameters["convergence_criterion"].GetString()
        if "contact" in self.convergence_criterion_name:
            D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
            D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
            R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
            R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
            CD_RT = convergence_criterion_parameters["contact_displacement_relative_tolerance"].GetDouble()
            CD_AT = convergence_criterion_parameters["contact_displacement_absolute_tolerance"].GetDouble()
            CR_RT = convergence_criterion_parameters["contact_residual_relative_tolerance"].GetDouble()
            CR_AT = convergence_criterion_parameters["contact_residual_absolute_tolerance"].GetDouble()
            condn_convergence_criterion = convergence_criterion_parameters["condn_convergence_criterion"].GetBool()
            self.fancy_convergence_criterion = convergence_criterion_parameters["fancy_convergence_criterion"].GetBool()
            self.print_convergence_criterion = convergence_criterion_parameters["print_convergence_criterion"].GetBool()
            self.mortar_type = convergence_criterion_parameters["mortar_type"].GetString()
            ensure_contact = convergence_criterion_parameters["ensure_contact"].GetBool()
            self.gidio_debug = convergence_criterion_parameters["gidio_debug"].GetBool()

            if(self.echo_level >= 1):
                KM.Logger.PrintInfo("::[Mechanical Solver]:: ", "CONVERGENCE CRITERION : " + self.convergence_criterion_name)

            if (self.fancy_convergence_criterion == True):
                self.table = KM.TableStreamUtility()

            if("contact_displacement_criterion" in self.convergence_criterion_name):
                if (self.fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact, self.table, self.print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)

            elif("contact_residual_criterion" in self.convergence_criterion_name):
                if (self.fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.table, self.print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)

            elif("contact_mixed_criterion" in self.convergence_criterion_name):
                if (self.fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.table, self.print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = CSMA.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)

            elif("contact_and_criterion" in self.convergence_criterion_name):
                if (self.fancy_convergence_criterion == True):
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, self.table, self.print_convergence_criterion)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.table, self.print_convergence_criterion)
                else:
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)

                Displacement.SetEchoLevel(self.echo_level)
                Residual.SetEchoLevel(self.echo_level)
                self.mechanical_convergence_criterion = KM.AndCriteria(Residual, Displacement)

            elif("contact_or_criterion" in self.convergence_criterion_name):
                if (self.fancy_convergence_criterion == True):
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, self.table, self.print_convergence_criterion)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, self.table, self.print_convergence_criterion)
                else:
                    Displacement = CSMA.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT,ensure_contact)
                    Residual = CSMA.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)

                Displacement.SetEchoLevel(self.echo_level)
                Residual.SetEchoLevel(self.echo_level)
                self.mechanical_convergence_criterion = KM.OrCriteria(Residual, Displacement)

            # Adding the mortar criteria
            Mortar = self.GetMortarCriteria()

            if (self.fancy_convergence_criterion == True):

                if (condn_convergence_criterion == True):
                    # Construct the solver
                    import eigen_solver_factory
                    settings_max = KM.Parameters("""
                    {
                        "solver_type"             : "power_iteration_highest_eigenvalue_solver",
                        "max_iteration"           : 10000,
                        "tolerance"               : 1e-9,
                        "required_eigen_number"   : 1,
                        "verbosity"               : 0,
                        "linear_solver_settings"  : {
                            "solver_type"             : "SuperLUSolver",
                            "max_iteration"           : 500,
                            "tolerance"               : 1e-9,
                            "scaling"                 : false,
                            "verbosity"               : 0
                        }
                    }
                    """)
                    eigen_solver_max = eigen_solver_factory.ConstructSolver(settings_max)
                    settings_min = KM.Parameters("""
                    {
                        "solver_type"             : "power_iteration_eigenvalue_solver",
                        "max_iteration"           : 10000,
                        "tolerance"               : 1e-9,
                        "required_eigen_number"   : 1,
                        "verbosity"               : 0,
                        "linear_solver_settings"  : {
                            "solver_type"             : "SuperLUSolver",
                            "max_iteration"           : 500,
                            "tolerance"               : 1e-9,
                            "scaling"                 : false,
                            "verbosity"               : 0
                        }
                    }
                    """)
                    eigen_solver_min = eigen_solver_factory.ConstructSolver(settings_min)

                    condition_number_utility = KM.ConditionNumberUtility(eigen_solver_max, eigen_solver_min)
                else:
                    condition_number_utility = None

                self.mechanical_convergence_criterion = CSMA.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar, self.table, self.print_convergence_criterion, condition_number_utility)
            else:
                self.mechanical_convergence_criterion = CSMA.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar)
            self.mechanical_convergence_criterion.SetEchoLevel(self.echo_level)
            self.mechanical_convergence_criterion.SetActualizeRHSFlag(True)

        else: # Standard criteria (same as structural mechanics application)
            # Construction of the class convergence_criterion
            import convergence_criteria_factory
            base_mechanical_convergence_criterion = convergence_criteria_factory.convergence_criterion(convergence_criterion_parameters)

            # Adding the mortar criteria
            Mortar = self.GetMortarCriteria(False)
            if ("ALMContact" in self.mortar_type or "MeshTying" in self.mortar_type):
                self.mechanical_convergence_criterion = KM.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            else:
                self.mechanical_convergence_criterion = base_mechanical_convergence_criterion.mechanical_convergence_criterion

    def GetMortarCriteria(self, include_table = True):
        # Adding the mortar criteria
        if (self.mortar_type == "ALMContactFrictionless"):
            if (self.fancy_convergence_criterion is True and include_table is True):
                Mortar = CSMA.ALMFrictionlessMortarConvergenceCriteria(self.table, self.print_convergence_criterion, self.gidio_debug)
            else:
                Mortar = CSMA.ALMFrictionlessMortarConvergenceCriteria()
        elif (self.mortar_type == "ALMContactFrictionlessComponents"):
            if (self.fancy_convergence_criterion is True and include_table is True):
                Mortar = CSMA.ALMFrictionlessComponentsMortarConvergenceCriteria(self.table, self.print_convergence_criterion, self.gidio_debug)
            else:
                Mortar = CSMA.ALMFrictionlessComponentsMortarConvergenceCriteria()
        elif (self.mortar_type == "ALMContactFrictional"):
            if (self.fancy_convergence_criterion is True and include_table is True):
                Mortar = CSMA.ALMFrictionalMortarConvergenceCriteria(self.table, self.print_convergence_criterion, self.gidio_debug)
            else:
                Mortar = CSMA.ALMFrictionalMortarConvergenceCriteria()
        elif ("MeshTying" in self.mortar_type):
            if (self.fancy_convergence_criterion is True and include_table is True):
                Mortar = CSMA.MeshTyingMortarConvergenceCriteria(self.table)
            else:
                Mortar = CSMA.MeshTyingMortarConvergenceCriteria()

        Mortar.SetEchoLevel(self.echo_level)

        return Mortar
