/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 3/21/22 3:17 PM                                                                         *
 * @Last Modified by: pravendra12                                                                    *
 * @Last Modified time: 11/28/24 8:30 PM                                                            *
 **************************************************************************************************/

/*! \file  Home.h
 *  \brief File for the Home class implementation.
 */

#include "Home.h"

namespace api
{

  void Print(const Parameter &parameter)
  {
    cout << "Parameters" << endl;
    cout << "simulation_method: " << parameter.method << endl;

    if (parameter.method == "KineticMcFirstMpi" || parameter.method == "KineticMcFirstOmp" ||
        parameter.method == "KineticMcChainOmp" || parameter.method == "KineticMcChainMpi" ||
        parameter.method == "KineticMcChainOmpi")
    {

      cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
           << endl;
      cout << "config_filename: " << parameter.config_filename_ << endl;
      cout << "map_filename: " << parameter.map_filename_ << endl;
      cout << "supercell_size: " << parameter.supercell_size_ << endl;
      cout << "lattice_param: " << parameter.lattice_param_ << endl;
      cout << "structure_type: " << parameter.structure_type_ << endl;
      cout << "log_dump_steps: " << parameter.log_dump_steps_ << endl;
      cout << "config_dump_steps: " << parameter.config_dump_steps_ << endl;
      cout << "maximum_steps: " << parameter.maximum_steps_ << endl;
      cout << "max_cluster_size: " << parameter.max_cluster_size_ << endl;
      cout << "max_bond_order: " << parameter.max_bond_order_ << endl;
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;
      cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_
           << endl;
      cout << "temperature: " << parameter.temperature_ << endl;
      cout << "restart_steps: " << parameter.restart_steps_ << endl;
      cout << "restart_energy: " << parameter.restart_energy_ << endl;
      cout << "restart_time: " << parameter.restart_time_ << endl;
      cout << "time_temperature_filename: " << parameter.time_temperature_filename_
           << endl;
      cout << "rate_corrector: " << parameter.rate_corrector_ << endl;
      cout << "vacancy_trajectory: " << parameter.vacancy_trajectory_
           << endl;
    }

    else if (parameter.method == "SimulatedAnnealing")
    {
      // General configuration information
      cout << "config_filename: "
           << parameter.config_filename_ << endl;
      cout << "json_coefficients_filename: "
           << parameter.json_coefficients_filename_ << endl;

      // Lattice and structure information
      cout << "lattice_param: "
           << parameter.lattice_param_ << endl;
      cout << "structure_type: "
           << parameter.structure_type_ << endl;

      // Simulation parameters
      cout << "maximum_steps: "
           << parameter.maximum_steps_ << endl;
      cout << "log_dump_steps: "
           << parameter.log_dump_steps_ << endl;
      cout << "config_dump_steps: "
           << parameter.config_dump_steps_ << endl;
      cout << "restart_steps: "
           << parameter.restart_steps_ << endl;
      cout << "restart_energy: "
           << parameter.restart_energy_ << endl;

      // Cluster expansion settings
      cout << "max_cluster_size: "
           << parameter.max_cluster_size_ << endl;
      cout << "max_bond_order: "
           << parameter.max_bond_order_ << endl;
      // Supercell size used for CE model fitting
      cout << "supercell_size: "
           << parameter.supercell_size_ << endl;

      // Cutoffs
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;

      // Temperature settings
      cout << "initial_temperature: "
           << parameter.initial_temperature_ << endl;
    }

    else if (parameter.method == "CanonicalMcSerial" || parameter.method == "CanonicalMcOmp")
    {
      cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
           << endl;
      cout << "config_filename: " << parameter.config_filename_ << endl;
      cout << "map_filename: " << parameter.map_filename_ << endl;
      cout << "supercell_size: " << parameter.supercell_size_ << endl;
      cout << "lattice_param: " << parameter.lattice_param_ << endl;
      cout << "structure_type: " << parameter.structure_type_ << endl;
      cout << "max_cluster_size: " << parameter.max_cluster_size_ << endl;
      cout << "max_bond_order: " << parameter.max_bond_order_ << endl;
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;
      cout << "log_dump_steps: " << parameter.log_dump_steps_ << endl;
      cout << "config_dump_steps: " << parameter.config_dump_steps_ << endl;
      cout << "maximum_steps: " << parameter.maximum_steps_ << endl;
      cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_
           << endl;
      cout << "temperature: " << parameter.temperature_ << endl;
      cout << "restart_steps: " << parameter.restart_steps_ << endl;
      cout << "restart_energy: " << parameter.restart_energy_ << endl;
      // cout << "initial_temperature: " << parameter.initial_temperature_ << endl;
      // cout << "decrement_temperature: " << parameter.decrement_temperature_ << endl;
    }

    else if (parameter.method == "Ansys" || parameter.method == "Reformat")
    {
      cout << "initial_steps: " << parameter.initial_steps_ << endl;
      cout << "increment_steps: " << parameter.increment_steps_ << endl;
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;
      cout << "log_type: " << parameter.log_type_ << endl;
      cout << "config_type: " << parameter.config_type_ << endl;

      cout << "Parameters for CE Encoding" << endl;
      cout << "extract_encoding: " << parameter.extract_encoding_ << endl;
      cout << "max_bond_order: " << parameter.max_bond_order_ << endl;
      cout << "max_cluster_size: " << parameter.max_cluster_size_ << endl;
    }
  }

  void Run(const Parameter &parameter)
  {
    if (parameter.method == "CanonicalMcSerial")
    {
      auto canonical_mc_serial = api::BuildCanonicalMcSerialFromParameter(parameter);
      canonical_mc_serial.Simulate();
    }
    else if (parameter.method == "KineticMcChainOmpi")
    {
      auto kinetic_mc_chain_ompi = api::BuildKineticMcChainOmpiFromParameter(parameter);
      cout << "Built successfully" << endl;
      kinetic_mc_chain_ompi.Simulate();
    }
    else if (parameter.method == "KineticMcFirstMpi")
    {
      auto kinetic_mc_first_mpi = api::BuildKineticMcFirstMpiFromParameter(parameter);
      kinetic_mc_first_mpi.Simulate();
    }
    else if (parameter.method == "Ansys")
    {
      auto iterator = api::BuildIteratorFromParameter(parameter);
      iterator.RunAnsys();
    }
    else if (parameter.method == "SimulatedAnnealing")
    {
      auto simulated_annealing = api::BuildSimulatedAnnealingFromParameter(parameter);
      simulated_annealing.Simulate();
    }
  }

  mc::CanonicalMcSerial BuildCanonicalMcSerialFromParameter(const Parameter &parameter)
  {
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);

      config.UpdateNeighborList(parameter.cutoffs_);
    }
    else
    {
      // Need to figure this out
      config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
    }

    auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                      parameter.lattice_param_,
                                                      "Mo",
                                                      parameter.structure_type_);

    // supercell_config.UpdateNeighborList(parameter.cutoffs_);

    // Getting element set from the configuration
    auto atom_vector = config.GetAtomVector();
    set<Element> element_set(atom_vector.begin(), atom_vector.end());

    // Print the elements in the set
    cout << "element_set: ";
    for (const auto &element : element_set)
    {
      cout << element.GetElementString() << " ";
    }
    cout << endl;

    cout << "Finish config reading. Start CMC." << endl;


    return mc::CanonicalMcSerial{config,
                                 supercell_config,
                                 parameter.log_dump_steps_,
                                 parameter.config_dump_steps_,
                                 parameter.maximum_steps_,
                                 parameter.thermodynamic_averaging_steps_,
                                 parameter.restart_steps_,
                                 parameter.restart_energy_,
                                 parameter.temperature_,
                                 ceParams};
  }

  mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter
                                                                  &parameter)
  {
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);

      config.UpdateNeighborList(parameter.cutoffs_);
    }
    else
    {
      // Need to figure this out
      config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
    }

    auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                      parameter.lattice_param_,
                                                      "X",
                                                      parameter.structure_type_);

    supercell_config.UpdateNeighborList(parameter.cutoffs_);

    // Getting element set from the configuration
    auto atom_vector = config.GetAtomVector();
    set<Element> element_set(atom_vector.begin(), atom_vector.end());

    // Print the elements in the set
    cout << "element_set: ";
    for (const auto &element : element_set)
    {
      cout << element.GetElementString() << " ";
    }
    cout << endl;

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    cout << "Finish config reading. Start CMC." << endl;


    return mc::CanonicalMcSerial{config,
                                 supercell_config,
                                 parameter.log_dump_steps_,
                                 parameter.config_dump_steps_,
                                 parameter.maximum_steps_,
                                 parameter.thermodynamic_averaging_steps_,
                                 parameter.restart_steps_,
                                 parameter.restart_energy_,
                                 parameter.temperature_,
                                 ceParams};
  }

  mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter
                                                                  &parameter)
  {
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);

      config.UpdateNeighborList(parameter.cutoffs_);
    }
    else
    {
      // Need to figure this out
      config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
    }

    auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                      parameter.lattice_param_,
                                                      "X",
                                                      parameter.structure_type_);

    supercell_config.UpdateNeighborList(parameter.cutoffs_);

    // Getting element set from the configuration
    auto atom_vector = config.GetAtomVector();
    set<Element> element_set(atom_vector.begin(), atom_vector.end());

    // Print the elements in the set
    cout << "element_set: ";
    for (const auto &element : element_set)
    {
      cout << element.GetElementString() << " ";
    }
    cout << endl;

=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

    cout << "Finish config reading. Start KMC." << endl;

    return mc::KineticMcChainOmpi{config,
                                  supercell_config,
                                  parameter.log_dump_steps_,
                                  parameter.config_dump_steps_,
                                  parameter.maximum_steps_,
                                  parameter.thermodynamic_averaging_steps_,
                                  parameter.restart_steps_,
                                  parameter.restart_energy_,
                                  parameter.restart_time_,
                                  parameter.temperature_,
                                  ceParams,
                                  parameter.time_temperature_filename_,
                                  parameter.rate_corrector_,
                                  parameter.vacancy_trajectory_};
  }

  mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter(const Parameter
                                                                &parameter)
  {
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);

      config.UpdateNeighborList(parameter.cutoffs_);
    }
    else
    {
      // Need to figure this out
      config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
    }

    auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                      parameter.lattice_param_,
                                                      "X",
                                                      parameter.structure_type_);

    supercell_config.UpdateNeighborList(parameter.cutoffs_);

    // Getting element set from the configuration
    auto atom_vector = config.GetAtomVector();
    set<Element> element_set(atom_vector.begin(), atom_vector.end());

    // Print the elements in the set
    cout << "element_set: ";
    for (const auto &element : element_set)
    {
      cout << element.GetElementString() << " ";
    }
    cout << endl;

    cout << "Finish config reading. Start KMC." << endl;

    return mc::KineticMcFirstMpi{config,
                                 supercell_config,
                                 parameter.log_dump_steps_,
                                 parameter.config_dump_steps_,
                                 parameter.maximum_steps_,
                                 parameter.thermodynamic_averaging_steps_,
                                 parameter.restart_steps_,
                                 parameter.restart_energy_,
                                 parameter.restart_time_,
                                 parameter.temperature_,
                                 ceParams,
                                 parameter.time_temperature_filename_,
                                 parameter.rate_corrector_,
                                 parameter.vacancy_trajectory_};
  }

  ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter)
  {

    return ansys::Traverse{parameter.initial_steps_,
                           parameter.increment_steps_,
                           parameter.cutoffs_,
                           parameter.log_type_,
                           parameter.config_type_,
                           parameter.extract_encoding_,
                           parameter.max_bond_order_,
                           parameter.max_cluster_size_};
  }

  mc::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter
                                                              &parameter)
  {
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    Config config;
    try
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);

      config.UpdateNeighborList(parameter.cutoffs_);
    }

    catch (...)
    { // Catch all other exceptions
      cerr << "File cannot be read" << endl;
      exit(EXIT_FAILURE);
    }

    auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                      parameter.lattice_param_,
                                                      "X", // can be any element
                                                      parameter.structure_type_);

    supercell_config.UpdateNeighborList(parameter.cutoffs_);

    // Getting element set from the configuration
    auto atom_vector = config.GetAtomVector();
    set<Element> element_set(atom_vector.begin(), atom_vector.end());

    // Print the elements in the set
    cout << "element_set: ";
    for (const auto &element : element_set)
    {
      cout << element.GetElementString() << " ";
    }
    cout << endl;

    cout << "Finish config reading. Start SA." << endl;

    return mc::SimulatedAnnealing{config,
                              supercell_config,
                              parameter.log_dump_steps_,
                              parameter.config_dump_steps_,
                              parameter.maximum_steps_,
                              parameter.restart_steps_,
                              parameter.restart_energy_,
                              parameter.initial_temperature_,
                              ceParams};
  }

}