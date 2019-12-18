#include <iostream>
#include <string>
#include <libconfig.h++>
#include <sys/stat.h>

#include "physics/cvode_gas.h"
#include "Tools/Common/common.h"


//----------------------------------------------------------------------------
// Convenience function to check if file exists
//----------------------------------------------------------------------------

inline bool file_exists(const std::string& name)
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

//----------------------------------------------------------------------------
// structs for initialization from libconfig file
//----------------------------------------------------------------------------
struct mechanism_options_t
{
  std::string filename;
  std::string name;
  void print_self()
  {
    std::cout << "Mechanism options:"
      << "\n + filename: "   << this->filename
      << "\n + name: " << this->name
      << std::endl;
  }
};
//----------------------------------------------------------------------------
struct integrator_options_t
{
  int num_steps;
  double integration_time_per_step;
  int solver_type;
  int jacobian_type;
  double relative_tolerance;
  double absolute_tolerance;
  void print_self()
  {
    std::cout << "Integrator options:"
      << "\n + num_steps: "                 << this->num_steps
      << "\n + integration_time_per_step: " << this->integration_time_per_step
      << "\n + solver_type: "               << this->solver_type
      << "\n + jacobian_type: "             << this->jacobian_type
      << "\n + relative_tolerance: "        << this->relative_tolerance
      << "\n + absolute_tolerance: "        << this->absolute_tolerance
      << std::endl;
  }
};
//----------------------------------------------------------------------------
struct initial_condition_t
{
  double temperature_ref;
  double pressure_ref;
  std::vector<int> non_zero_species_index;
  std::vector<double> non_zero_species_mass_fraction;
  void print_self()
  {
    std::cout << "Initial condition:"
      << "\n + temperature_ref "   << this->temperature_ref
      << "\n + pressure_ref "   << this->pressure_ref
      << "\n + non-zero species (" << this->non_zero_species_index.size()
      <<") mass fractions:";
    for (int sdx = 0; sdx < this->non_zero_species_index.size(); sdx++)
    {
      std::cout << "\n   * " << sdx << " : "
                << this->non_zero_species_mass_fraction[sdx];
    }
    std::cout << std::endl;
  }
};
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Function to read settings from libconfig file
//----------------------------------------------------------------------------
void Read_libconfig_file(
  const std::string input_filename_In,
  mechanism_options_t & mechanism_options_Out,
  integrator_options_t & integrator_options_Out,
  initial_condition_t & initial_condition_Out)
{

  // flowControl while (for error handling)
  int flowControl = 1;
  while ( flowControl == 1 )
  {

  // Check if file exists
  if ( not file_exists(input_filename_In) )
  {
    common::errorHandler->SendError(_FILE_LINE + 
        " - input filename '" + input_filename_In + "' does not exist");
    break; // Exit flowControl with error
  }

  libconfig::Config input_options;
  input_options.readFile(input_filename_In.c_str());
  const libconfig::Setting & settings = input_options.getRoot();

  try
  {

    // 1) Mechanism
    // ------------
    const libconfig::Setting &mechanism_settings_In = settings["mechanism"];
    if ( not mechanism_settings_In.lookupValue(
               "file", mechanism_options_Out.filename) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - mechanism 'file' must be specified");
      break; // Exit flowControl with error
    }
    if ( not mechanism_settings_In.lookupValue(
               "name", mechanism_options_Out.name) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - mechanism 'name' must be specified");
      break; // Exit flowControl with error
    }

    // 2) Integrator_options
    // ---------------------
    const libconfig::Setting &integrator_settings_In = settings["integrator"];
    if ( not integrator_settings_In.lookupValue(
               "num_steps", integrator_options_Out.num_steps) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - integrator 'num_steps' must be specified");
      break; // Exit flowControl with error
    }
    if ( not integrator_settings_In.lookupValue(
               "integration_time_per_step",
               integrator_options_Out.integration_time_per_step) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - integrator 'integration_time_per_step' must be specified");
      break; // Exit flowControl with error
    }
    std::string integrator_solver_type;
    std::string integrator_jacobian_type;
    if ( not integrator_settings_In.lookupValue(
               "solver_type", integrator_solver_type) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - integrator 'solver_type' must be specified");
      break; // Exit flowControl with error
    }
    if ( not integrator_settings_In.lookupValue(
               "jacobian_type", integrator_jacobian_type) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - integrator 'jacobian_type' must be specified");
      break; // Exit flowControl with error
    }
    if ( integrator_solver_type == "DIAG" )
    {
      integrator_options_Out.solver_type = math::DIAG;
      if ( integrator_jacobian_type == "JAC" )
      {
        integrator_options_Out.jacobian_type = math::JAC;
      }
      else
      {
        common::errorHandler->SendError(_FILE_LINE + 
            " - 'DIAG' solver type can only take 'JAC' jacobian type");
        break; // Exit flowControl with error
      }
    }
    else if ( integrator_solver_type == "BAND" )
    {
      integrator_options_Out.solver_type = math::BAND;
      if ( integrator_jacobian_type == "NOJAC" )
      {
        integrator_options_Out.jacobian_type = math::NOJAC;
      }
      else
      {
        common::errorHandler->SendError(_FILE_LINE + 
            " - 'DIAG' solver type can only take 'NOJAC' jacobian type");
        break; // Exit flowControl with error
      }
    }
    else if ( integrator_solver_type == "DENSE" )
    {
      integrator_options_Out.solver_type = math::DENSE;
      if ( integrator_jacobian_type == "JAC" )
      {
        integrator_options_Out.jacobian_type = math::JAC;
      }
      else if ( integrator_jacobian_type == "NOJAC" )
      {
        integrator_options_Out.jacobian_type = math::NOJAC;
      }
      else
      {
        common::errorHandler->SendError(_FILE_LINE + 
          " - Unknown integrator jacobian type '" +
          integrator_jacobian_type + "'");
        break; // Exit flowControl with error
      }
    }
    else
    {
      common::errorHandler->SendError(_FILE_LINE + 
        " - Unknown integrator solver type '" + integrator_solver_type + "'");
      break; // Exit flowControl with error
    }

    // Tolerances
    if ( not integrator_settings_In.lookupValue(
               "relative_tolerance",
               integrator_options_Out.relative_tolerance) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - integrator 'relative_tolerance' must be specified");
      break; // Exit flowControl with error
    }
    if ( not integrator_settings_In.lookupValue(
               "absolute_tolerance",
               integrator_options_Out.absolute_tolerance) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - integrator 'absolute_tolerance' must be specified");
      break; // Exit flowControl with error
    }

    // 3) Initial condition
    // --------------------
    const libconfig::Setting & initial_condition_settings_In =
      settings["initial_condition"];
    if ( not initial_condition_settings_In.lookupValue(
               "pressure_ref", initial_condition_Out.pressure_ref) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - initial_condition 'pressure_ref' must be specified");
      break; // Exit flowControl with error
    }
    if ( not initial_condition_settings_In.lookupValue(
               "temperature_ref", initial_condition_Out.temperature_ref) )
    {
      common::errorHandler->SendError(_FILE_LINE + 
          " - initial_condition 'temperature_ref' must be specified");
      break; // Exit flowControl with error
    }
    const libconfig::Setting & indexes =
      initial_condition_settings_In["non_zero_species"]["indexes"];
    for (int idx = 0; idx < indexes.getLength(); idx++)
    {
      initial_condition_Out.non_zero_species_index.push_back(indexes[idx]);
    }
    const libconfig::Setting & mass_fractions =
      initial_condition_settings_In["non_zero_species"]["mass_fractions"];
    for (int idx = 0; idx < mass_fractions.getLength(); idx++)
    {
      initial_condition_Out.non_zero_species_mass_fraction.push_back(
        mass_fractions[idx]);
    }

  }
  catch (const libconfig::SettingNotFoundException &nfex)
  {
    // Ignore.
  }
  catch (const libconfig::ParseException &nfex)
  {
    common::errorHandler->SendError(_FILE_LINE +
      " - Exception caught with 'what':" + nfex.what());
    throw;
  }
  catch (const libconfig::SettingTypeException &nfex)
  {
    common::errorHandler->SendError(_FILE_LINE +
      " - Exception caught with 'what':" + nfex.what());
    throw;
  }

  // Continue flow normally (no unrecoverable errors up to this point)
  flowControl++;
  }

}

//----------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  // flowControl while (for error handling)
  int flowControl = 1;
  while ( flowControl == 1 )
  {

  std::cout << "Stand alone test." << std::endl;
  if (argc != 2)
  {
    common::errorHandler->SendError(_FILE_LINE +
      " - Usage: program input_filename");
  }
  std::string input_filename = argv[1];
  initial_condition_t initial_condition;
  mechanism_options_t mechanism_options;
  integrator_options_t integrator_options;
  Read_libconfig_file(
    input_filename,
    mechanism_options,
    integrator_options,
    initial_condition);
  if ( common::errorHandler->GetErrorStatus() )
  {
    common::errorHandler->SendError(" - Reading input file " + input_filename);
    break; // Exit flowControl with error
  }

  // report read options
  mechanism_options.print_self();
  integrator_options.print_self();
  initial_condition.print_self();

  // Initialize cvode_gas
  physics::CVode_gas * cvode_gas = physics::CVode_gas::New(); 
  cvode_gas->Set_input(mechanism_options.filename, mechanism_options.name);
  cvode_gas->Set_integrator_parameters(
    integrator_options.solver_type + integrator_options.jacobian_type,
    integrator_options.relative_tolerance,
    integrator_options.absolute_tolerance);
  cvode_gas->Init();
  const int num_species = cvode_gas->Get_nscalar();
  std::cout << "num_species = " << num_species << std::endl;
  std::vector<double> transport(num_species, 0.0);

  double pressure    = initial_condition.pressure_ref;
  double temperature = initial_condition.temperature_ref;
  // Initialize species with non-zero mass fractions in the mixture
  std::vector<double> species_mass_fractions(num_species,0.0);
  for (int idx = 0 ;
       idx < initial_condition.non_zero_species_index.size();
       idx++)
  {
    species_mass_fractions[initial_condition.non_zero_species_index[idx]] =
      initial_condition.non_zero_species_mass_fraction[idx];
  }
  for (int tdx = 0; tdx < integrator_options.num_steps; tdx++)
  {
    std::cout << "step = " << tdx << std::endl;
    cvode_gas->Set_transport_vector(transport);
    std::cout << "setting state via T, p, and scalars" << std::endl;
    cvode_gas->Set_state_temperature_pressure_scalars(
      temperature,
      pressure,
      species_mass_fractions.data());
    cvode_gas->zeroD_integrator(integrator_options.integration_time_per_step);
    pressure    = cvode_gas->Get_pressure();
    temperature = cvode_gas->Get_temperature();
    cvode_gas->Get_mass_fractions(&species_mass_fractions[0]) ;
    // Report outcome
    std::cout << "pressure = "    << pressure    << std::endl;
    std::cout << "temperature = " << temperature << std::endl;
    for (int sdx = 0; sdx < num_species; sdx++)
    {
      std::cout << "mass_fraction[" << sdx << "] = "
                << species_mass_fractions[sdx]
                << std::endl;
    }
  }
  cvode_gas->Delete();

  // Continue flow normally (no unrecoverable errors up to this point)
  flowControl++;
  }
  return 0;

}
