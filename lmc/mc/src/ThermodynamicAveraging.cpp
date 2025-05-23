/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/30/23 3:09 PM                                                          *
 **************************************************************************************************/

#include "ThermodynamicAveraging.h"
#include <cmath>
#include "Constants.hpp"
#include <iostream>

namespace mc
{
  mc::ThermodynamicAveraging::ThermodynamicAveraging(size_t size)
      : size_(size)
  {}

  void ThermodynamicAveraging::AddEnergy(double value)
  {
    if (energy_list_.size() == size_)
    {
      // std::cout << "I am here in inside Add Energy" << std::endl;
      sum_ -= energy_list_.front();
      energy_list_.pop_front();
    }
    // std::cout << "i m here at loc3" << std::endl;
    energy_list_.push_back(value);
    sum_ += value;
  }

  double ThermodynamicAveraging::GetAverage() const
  {
    if (energy_list_.empty())
    {
      return 0;
    }
    return sum_ / static_cast<double>(energy_list_.size());
  }

  double ThermodynamicAveraging::GetThermodynamicAverage(double beta) const
  {
    const auto average = GetAverage();
    double partition = 0.0;
    double thermodynamic_average_energy = 0.0;
    
    for (auto energy : energy_list_)
    {
      energy -= average;
      const double exp_value = std::exp(-energy * beta);
      thermodynamic_average_energy += energy * exp_value;
      partition += exp_value;
    }
    return thermodynamic_average_energy / partition + average;
  }
} // mc
