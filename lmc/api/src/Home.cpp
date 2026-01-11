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
      cout << "enable_sro: " << parameter.enable_sro_ << endl;
      cout << "enable_b2_order_param: " << parameter.enable_b2_order_param_ << endl;
      cout << "enable_b2_cluster_ansys: " << parameter.enable_b2_cluster_ansys_ << endl;
    }
  }

  void Run(const Parameter &parameter)
  {
    if (parameter.method == "CanonicalMcSerial")
    {
      api::RunCanonicalMcSerialFromParameter(parameter);
    }
    else if (parameter.method == "KineticMcChainOmpi")
    {
      api::RunKineticMcChainOmpiFromParameter(parameter);
    }
    else if (parameter.method == "KineticMcFirstMpi")
    {
      api::RunKineticMcFirstMpiFromParameter(parameter);
    }
    else if (parameter.method == "KineticMcFirstOmp")
    {
      api::RunKineticMcFirstOmpFromParameter(parameter);
    }
    else if (parameter.method == "Ansys")
    {
      auto iterator = api::BuildIteratorFromParameter(parameter);
      iterator.RunAnsys();
    }
    else if (parameter.method == "SimulatedAnnealing")
    {
      api::RunSimulatedAnnealingFromParameter(parameter);
    }
  }

  void RunCanonicalMcSerialFromParameter(const Parameter &parameter)
  {
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);
    }

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    config.UpdateNeighborList({maxClusterCutoff});

    // Generate a small config for declaring symCE
    const size_t supercellSize = 2;

    auto primConfig = Config::GenerateSupercell(
        supercellSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictor symCEEnergyPredictor(
        ceParams,
        config,
        primConfig);

    // Again update the neighbor list
    config.UpdateNeighborList(parameter.cutoffs_);

    // Declare LVFE Predictor
    // LVFEPredictor lvfePredictor(
    //     ceParams,
    //     config);

    // Declare Energy Predictor
    EnergyPredictor energyChangePredictor(
        symCEEnergyPredictor);

    cout << "Finish config reading. Start CMC." << endl;

    mc::CanonicalMcSerial canonicalMcSerial(config,
                                            parameter.log_dump_steps_,
                                            parameter.config_dump_steps_,
                                            parameter.maximum_steps_,
                                            parameter.thermodynamic_averaging_steps_,
                                            parameter.restart_steps_,
                                            parameter.restart_energy_,
                                            parameter.temperature_,
                                            energyChangePredictor);

    canonicalMcSerial.Simulate();

    cout << "Simulation Completed" << endl;
  }

  void RunKineticMcChainOmpiFromParameter(const Parameter
                                              &parameter)
  {
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);
    }

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    config.UpdateNeighborList({maxClusterCutoff});

    // Generate a small config for declaring symCE
    const size_t supercellSize = 2;

    auto primConfig = Config::GenerateSupercell(
        supercellSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictor symCEEnergyPredictor(
        ceParams,
        config,
        primConfig);

    // Again update the neighbor list
    config.UpdateNeighborList(parameter.cutoffs_);

    // Declare LVFE Predictor
    // LVFEPredictor lvfePredictor(
    //     ceParams,
    //     config);

    // Declare Energy Predictor
    EnergyPredictor energyChangePredictor(
        symCEEnergyPredictor);

    // Declare KRA Predictor
    KRAPredictor eKRAPredictor(
        ceParams,
        config);

    // Declare Vacacny migration predictor
    VacancyMigrationPredictor vacancyMigrationPredictor(
        eKRAPredictor,
        energyChangePredictor);

    cout << "Finish config reading. Start KMC." << endl;

    mc::KineticMcChainOmpi kmcChainOmpi(config,
                                        parameter.log_dump_steps_,
                                        parameter.config_dump_steps_,
                                        parameter.maximum_steps_,
                                        parameter.thermodynamic_averaging_steps_,
                                        parameter.restart_steps_,
                                        parameter.restart_energy_,
                                        parameter.restart_time_,
                                        parameter.temperature_,
                                        vacancyMigrationPredictor,
                                        parameter.time_temperature_filename_,
                                        parameter.rate_corrector_,
                                        parameter.vacancy_trajectory_);

    kmcChainOmpi.Simulate();

    cout << "Simulation Completed" << endl;
  }

  void RunKineticMcFirstMpiFromParameter(const Parameter
                                             &parameter)
  {
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);
    }

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    config.UpdateNeighborList({maxClusterCutoff});

    // Generate a small config for declaring symCE
    const size_t supercellSize = 2;

    auto primConfig = Config::GenerateSupercell(
        supercellSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictor symCEEnergyPredictor(
        ceParams,
        config,
        primConfig);

    // Again update the neighbor list
    config.UpdateNeighborList(parameter.cutoffs_);

    // Declare LVFE Predictor
    // LVFEPredictor lvfePredictor(
    //     ceParams,
    //     config);

    // Declare Energy Predictor
    EnergyPredictor energyChangePredictor(
        symCEEnergyPredictor);

    // Declare KRA Predictor
    KRAPredictor eKRAPredictor(
        ceParams,
        config);

    // Declare Vacacny migration predictor
    VacancyMigrationPredictor vacancyMigrationPredictor(
        eKRAPredictor,
        energyChangePredictor);

    cout << "Finish config reading. Start KMC." << endl;

    mc::KineticMcFirstMpi kmcFirstMpi(config,
                                      parameter.log_dump_steps_,
                                      parameter.config_dump_steps_,
                                      parameter.maximum_steps_,
                                      parameter.thermodynamic_averaging_steps_,
                                      parameter.restart_steps_,
                                      parameter.restart_energy_,
                                      parameter.restart_time_,
                                      parameter.temperature_,
                                      vacancyMigrationPredictor,
                                      parameter.time_temperature_filename_,
                                      parameter.rate_corrector_,
                                      parameter.vacancy_trajectory_);

    kmcFirstMpi.Simulate();

    cout << "Simulation Completed" << endl;
  }

  void RunKineticMcFirstOmpFromParameter(const Parameter
                                             &parameter)
  {
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);
    }

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    config.UpdateNeighborList({maxClusterCutoff});

    // Generate a small config for declaring symCE
    const size_t supercellSize = 2;

    auto primConfig = Config::GenerateSupercell(
        supercellSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictor symCEEnergyPredictor(
        ceParams,
        config,
        primConfig);

    // Again update the neighbor list
    config.UpdateNeighborList(parameter.cutoffs_);

    // Declare LVFE Predictor
    // LVFEPredictor lvfePredictor(
    //     ceParams,
    //     config);

    // Declare Energy Predictor
    EnergyPredictor energyChangePredictor(
        symCEEnergyPredictor);

    // Declare KRA Predictor
    KRAPredictor eKRAPredictor(
        ceParams,
        config);

    // Declare Vacacny migration predictor
    VacancyMigrationPredictor vacancyMigrationPredictor(
        eKRAPredictor,
        energyChangePredictor);

    cout << "Finish config reading. Start KMC." << endl;

    mc::KineticMcFirstOmp kmcFirstOmp(config,
                                      parameter.log_dump_steps_,
                                      parameter.config_dump_steps_,
                                      parameter.maximum_steps_,
                                      parameter.thermodynamic_averaging_steps_,
                                      parameter.restart_steps_,
                                      parameter.restart_energy_,
                                      parameter.restart_time_,
                                      parameter.temperature_,
                                      vacancyMigrationPredictor,
                                      parameter.time_temperature_filename_,
                                      parameter.rate_corrector_,
                                      parameter.vacancy_trajectory_);

    kmcFirstOmp.Simulate();

    cout << "Simulation Completed" << endl;
  }

  ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter)
  {
    AnsysFlags ansysFlags;
    ansysFlags.SRO = parameter.enable_sro_;
    ansysFlags.B2OrderParam = parameter.enable_b2_order_param_;
    ansysFlags.B2ClusterAnsys = parameter.enable_b2_cluster_ansys_;

    return ansys::Traverse{parameter.initial_steps_,
                           parameter.increment_steps_,
                           parameter.cutoffs_,
                           ansysFlags,
                           parameter.log_type_,
                           parameter.config_type_};
  }

  void RunSimulatedAnnealingFromParameter(const Parameter
                                              &parameter)
  {
    Config config;
    if (parameter.map_filename_.empty())
    {
      // Generalized function to read configuration
      // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
      config = Config::ReadConfig(parameter.config_filename_);
    }

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    config.UpdateNeighborList({maxClusterCutoff});

    // Generate a small config for declaring symCE
    const size_t supercellSize = 2;

    auto primConfig = Config::GenerateSupercell(
        supercellSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictor symCEEnergyPredictor(
        ceParams,
        config,
        primConfig);

    // Again update the neighbor list
    config.UpdateNeighborList(parameter.cutoffs_);

    // Declare LVFE Predictor
    // LVFEPredictor lvfePredictor(
    //     ceParams,
    //     config);

    // Declare Energy Predictor
    EnergyPredictor energyChangePredictor(
        symCEEnergyPredictor);

    cout << "Finish config reading. Start SA." << endl;

    mc::SimulatedAnnealing simulatedAnnealing(config,
                                              parameter.log_dump_steps_,
                                              parameter.config_dump_steps_,
                                              parameter.maximum_steps_,
                                              parameter.restart_steps_,
                                              parameter.restart_energy_,
                                              parameter.initial_temperature_,
                                              energyChangePredictor);

    simulatedAnnealing.Simulate();
    cout << "Simulation Completed" << endl;
  }
}