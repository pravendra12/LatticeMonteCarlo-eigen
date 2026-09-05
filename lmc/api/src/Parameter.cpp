/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 3/21/22 3:17 PM                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 11/28/24 8:20 PM                                                          *
 **************************************************************************************************/

/*! \file  Parameter.cpp
 *  \brief File for the Parameter Struct implementation.
 */

#include "Parameter.h"
#include "Utility.h"
#include <algorithm>

namespace api
{

  Parameter::Parameter(int argc, char *argv[])
  {
    ParseArgs(argc, argv);
    ReadParam(parameters_filename);
  }

  Parameter::Parameter(const std::string &param_filename)
  {
    ReadParam(param_filename);
  }

  void Parameter::ParseArgs(int argc, char *argv[])
  {
    for (int i = 0; i < argc; i++)
    {
      if (!std::strcmp(argv[i], "--p") || !std::strcmp(argv[i], "-p"))
        parameters_filename = std::string(argv[++i]);
    }
  }

  void Parameter::ReadParam(const std::string &param_filename)
  {
    if (param_filename.empty())
    {
      return;
    }
    std::ifstream ifs(param_filename, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw std::runtime_error("Cannot open " + param_filename);
    }
    std::string buffer;
    while (std::getline(ifs, buffer))
    {
      if (buffer.empty())
      {
        continue;
      }
      if (buffer[0] == '#')
      {
        continue;
      }
      std::vector<std::string> segs(split(buffer, " "));
      if (segs[0] == "simulation_method")
      {
        method = std::string(segs[1]);
      }
      else if (segs[0] == "config_filename")
      {
        config_filename_ = std::string(segs[1]);
      }
      else if (segs[0] == "map_filename")
      {
        map_filename_ = std::string(segs[1]);
      }
      else if (segs[0] == "json_coefficients_filename")
      {
        json_coefficients_filename_ = std::string(segs[1]);
      }
      else if (segs[0] == "time_temperature_filename")
      {
        time_temperature_filename_ = std::string(segs[1]);
      }
      else if (segs[0] == "log_type")
      {
        log_type_ = std::string(segs[1]);
      }
      else if (segs[0] == "config_type")
      {
        config_type_ = std::string(segs[1]);
      }
      else if (segs[0] == "structure_type")
      {
        structure_type_ = std::string(segs[1]);
      }
      else if (segs[0] == "log_dump_steps")
      {
        log_dump_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "config_dump_steps")
      {
        config_dump_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "maximum_steps")
      {
        maximum_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "max_cluster_size")
      {
        max_cluster_size_ = stoul(segs[1]);
      }
      else if (segs[0] == "max_bond_order")
      {
        max_bond_order_ = stoul(segs[1]);
      }
      else if (segs[0] == "thermodynamic_averaging_steps")
      {
        thermodynamic_averaging_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "temperature")
      {
        temperature_ = stod(segs[1]);
      }
      else if (segs[0] == "initial_temperature")
      {
        initial_temperature_ = stod(segs[1]);
      }
      else if (segs[0] == "decrement_temperature")
      {
        decrement_temperature_ = stod(segs[1]);
      }
      else if (segs[0] == "cutoffs")
      {
        cutoffs_.clear();
        std::transform(segs.begin() + 1,
                       segs.end(),
                       std::back_inserter(cutoffs_),
                       [](const auto &cutoff)
                       { return stod(cutoff); });
      }
      else if (segs[0] == "vacancy_trajectory")
      {
        if (segs.size() != 4)
        { // 1 for the keyword "cutoffs" + 3 values
          throw std::runtime_error("Invalid number of elements for Vacancy Trajectory. Expected 3 values.");
        }
        vacancy_trajectory_ << stod(segs[1]), stod(segs[2]), stod(segs[3]);
      }
      else if (segs[0] == "supercell_size")
      {
        supercell_size_ = stoul(segs[1]);
      }
      else if (segs[0] == "lattice_param")
      {
        lattice_param_ = stod(segs[1]);
      }
      else if (segs[0] == "initial_steps")
      {
        initial_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "increment_steps")
      {
        increment_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "smallest_cluster_criteria")
      {
        smallest_cluster_criteria_ = stoul(segs[1]);
      }
      else if (segs[0] == "solvent_bond_criteria")
      {
        solvent_bond_criteria_ = stoul(segs[1]);
      }
      else if (segs[0] == "restart_steps")
      {
        restart_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "rate_corrector")
      {
        std::string bool_string = std::string(segs[1]);
        if (bool_string == "true")
        {
          rate_corrector_ = true;
        }
        else
        {
          rate_corrector_ = false;
        }
      }
      else if (segs[0] == "restart_energy")
      {
        restart_energy_ = stod(segs[1]);
      }
      else if (segs[0] == "restart_time")
      {
        restart_time_ = stod(segs[1]);
      }
      else if (segs[0] == "factor")
      {
        factor_ = stoul(segs[1]);
      }
      else if (segs[0] == "solvent_element")
      {
        solvent_element_ = std::string(segs[1]);
      }
      else if (segs[0] == "solute_element_set")
      {
        solute_element_set_.clear();
        std::copy(segs.begin() + 1, segs.end(),
                  std::back_inserter(solute_element_set_));
      }
      else if (segs[0] == "solute_number_set")
      {
        solute_number_set_.clear();
        std::transform(segs.begin() + 1,
                       segs.end(),
                       std::back_inserter(solute_number_set_),
                       [](const auto &number)
                       { return stoul(number); });
      }
      else if (segs[0] == "early_stop_steps")
      {
        early_stop_steps_ = stoull(segs[1]);
      }
      else if (segs[0] == "extract_encoding")
      {
        std::string bool_string = std::string(segs[1]);
        if (bool_string == "true")
        {
          extract_encoding_ = true;
        }
        else
        {
          extract_encoding_ = false;
        }
      }
      else if (segs[0] == "enable_sro")
      {
        std::string bool_string = std::string(segs[1]);
        if (bool_string == "true")
        {
          enable_sro_ = true;
        }
        else
        {
          enable_sro_ = false;
        }
      }
      else if (segs[0] == "enable_b2_order_param")
      {
        std::string bool_string = std::string(segs[1]);
        if (bool_string == "true")
        {
          enable_b2_order_param_ = true;
        }
        else
        {
          enable_b2_order_param_ = false;
        }
      }
      else if (segs[0] == "enable_b2_cluster_ansys")
      {
        std::string bool_string = std::string(segs[1]);
        if (bool_string == "true")
        {
          enable_b2_cluster_ansys_ = true;
        }
        else
        {
          enable_b2_cluster_ansys_ = false;
        }
      }
      else if (segs[0] == "enable_b2_cluster_atom_map")
      {
        std::string bool_string = std::string(segs[1]);
        if (bool_string == "true")
        {
          enable_b2_cluster_atom_map_ = true;
        }
        else
        {
          enable_b2_cluster_atom_map_ = false;
        }
      }

      else if (segs[0] == "ghost_element")
      {
        ghost_element_ = std::string(segs[1]);
      }
    }
    ifs.close();
  }
} // api
