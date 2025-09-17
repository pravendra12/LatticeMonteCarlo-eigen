/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcOmp.h
 *  @brief File for CanonicalMcOmp class implementation.
 */

#include "CanonicalMcSerial.h"
#include <utility>
#include <chrono>
#include <omp.h>
namespace mc
{
  CanonicalMcSerial::CanonicalMcSerial(Config config,
                                       unsigned long long int logDumpSteps,
                                       unsigned long long int configDumpStep,
                                       unsigned long long int maximumSteps,
                                       unsigned long long int thermodynamicAveragingSteps,
                                       unsigned long long int restartSteps,
                                       double restartEnergy,
                                       double temperature,
                                       EnergyPredictor &energyChangePredictor)
      : CanonicalMcAbstract(move(config),
                            logDumpSteps,
                            configDumpStep,
                            maximumSteps,
                            thermodynamicAveragingSteps,
                            restartSteps,
                            restartEnergy,
                            temperature,
                            energyChangePredictor)
  {
    if (world_size_ != 1)
    {
      std::cout << "Must use 1 processes. Terminating...\n"
                << std::endl;
      MPI_Finalize();
      exit(0);
    }
#pragma omp parallel default(none) shared(std::cout)
    {
#pragma omp master
      {
        std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
      }
    }
  }

  void CanonicalMcSerial::Simulate()
  {
    while (steps_ <= maximumSteps_)
    {

      // auto latticeIdJumpPair = GenerateLatticeIdJumpPair();
      // One site would need to be a vacant site 
      // This can be made conditional based on the type of hamiltonian being used
    
      auto latticeIdJumpPair = GenerateVacancyLatticeIdJumpPair();

      auto dE = energyChangePredictor_.GetEnergyChange(config_,
                                                       latticeIdJumpPair);
      Dump();

      SelectEvent(latticeIdJumpPair, dE);

      ++steps_;
    }
  }
} // mc
