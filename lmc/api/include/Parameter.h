/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 3/21/22 3:17 PM                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 11/28/24 8:10 PM                                                          *
 **************************************************************************************************/

/*! \file  Parameter.h
 *  \brief File for the Parameter Struct declaration.
 */

#ifndef LMC_API_INCLUDE_PARAMETER_H_
#define LMC_API_INCLUDE_PARAMETER_H_

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <Eigen/Dense>

namespace api {

struct Parameter {
 public:

  /*!
   * \brief Constructs a Parameter object and parses command-line arguments.
   * \param : argc Number of command-line arguments.
   * \param : argv Array of command-line argument strings.
   */
  Parameter(int argc, char *argv[]);

  /*!
   * \brief Constructs a Parameter object and reads parameters from a file.
   * \param param_filename : Parameter file name.
   */
  explicit Parameter(const std::string &param_filename);

  /*!
   * \brief Parses command-line arguments and initializes parameters.
   * \param argc : Number of command-line arguments.
   * \param argv : Array of command-line argument strings.
   */
  void ParseArgs(int argc, char *argv[]);

  /*!
   * \brief Reads and initializes parameters from a file.
   * \param param_filename : Parameter file name.
   */
  void ReadParam(const std::string &param_filename);

  /// Member variables

  /// Filename containing simulation parameters.
  std::string parameters_filename{};

  /// Method used for the simulation (e.g. CanonicalMcSerial, KineticMcChainOmpi).
  std::string method{};

  /// Configuration file name.
  std::string config_filename_{};

  /// Filename for the mapping file.
  std::string map_filename_{};

  /// JSON file containing pre-trained coefficients.
  std::string json_coefficients_filename_{};

  /// Filename for time-temperature data.
  std::string time_temperature_filename_{};

  /// Type of logging.
  std::string log_type_{};

  /// Type of configuration format (e.g., POSCAR, XYZ)       
  /// **************Need To Check Once **********
  std::string config_type_{};

  /// Type of structure (e.g., BCC, FCC).
  std::string structure_type_{};

  /// Number of steps after which the log file is updated.
  unsigned long long int log_dump_steps_{};

  /// Number of steps after which the configuration is dumped.
  unsigned long long int config_dump_steps_{};

  /// Maximum number of steps for the simulation.
  unsigned long long int maximum_steps_{};

  /// Vector to store vacancy trajectory
  Eigen::RowVector3d vacancy_trajectory_{};

  /// Maximum size of clusters to be considered.
  size_t max_cluster_size_{};

  /// Maximum bond order to be considered.
  size_t max_bond_order_{};

  /// Number of steps used for thermodynamic averaging.
  unsigned long long int thermodynamic_averaging_steps_{};

  /// Temperature of the simulation in Kelvin.
  double temperature_{};

  /// Initial Temperature in Kelvin.
  double initial_temperature_{};

  /// Decrement in temperature for annealing.
  double decrement_temperature_{};

  /// Size of the supercell used for training.
  size_t supercell_size_{};

  /// Lattice parameter of the structure.
  double lattice_param_{};

  /// Steps at which the simulation restarts.
  unsigned long long int restart_steps_{};

  /// Energy at the restart step.
  double restart_energy_{};

  /// Time at the restart step.
  double restart_time_{};

  /// Enables or disables rate correction.
  bool rate_corrector_{};

  /// Cutoff distances for cluster interactions.
  std::vector<double> cutoffs_{};

  /// Extract the local environment
  bool extract_encoding_{};
  
  /// Initial number of steps for the simulation.
  unsigned long long int initial_steps_{};

  /// Incremental steps for specific operations.
  unsigned long long int increment_steps_{};

  /// Criteria for the smallest clusters.
  size_t smallest_cluster_criteria_{};

  /// Criteria for solvent bonds.
  size_t solvent_bond_criteria_{};

  /// Scaling factor for specific parameters.
  size_t factor_{};

  /// Element representing the solvent in the simulation.
  std::string solvent_element_{};

  /// List of solute elements in the simulation.
  std::vector<std::string> solute_element_set_{};

  /// Number of each solute element in the system.
  std::vector<size_t> solute_number_set_{};

  /// Early stopping criteria based on the number of steps.
  unsigned long long int early_stop_steps_{};


  // Flags for ansys
  bool enable_sro_{};
  bool enable_b2_order_param_{};
  bool enable_b2_cluster_ansys_{};

  bool enable_b2_cluster_atom_map_{};


  // For widom insertion
  // element which will be inserted
  std::string ghost_element_{};
};
} // namespace api

#endif //LMC_API_INCLUDE_PARAMETER_H_
