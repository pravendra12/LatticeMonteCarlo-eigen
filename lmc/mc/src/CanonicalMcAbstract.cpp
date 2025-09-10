/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcAbstract.cpp
 *  @brief File for CanonicalMcAbstract class implementation.
 */

#include "CanonicalMcAbstract.h"

namespace mc
{
  CanonicalMcAbstract::CanonicalMcAbstract(Config config,
                                           Config supercellConfig,
                                           const unsigned long long int logDumpSteps,
                                           const unsigned long long int configDumpSteps,
                                           const unsigned long long int maximumSteps,
                                           const unsigned long long int thermodynamicAveragingSteps,
                                           const unsigned long long int restartSteps,
                                           const double restartEnergy,
                                           const double temperature,
                                           const ClusterExpansionParameters &ceParams)
      : McAbstract(move(config),
                   supercellConfig,
                   logDumpSteps,
                   configDumpSteps,
                   maximumSteps,
                   thermodynamicAveragingSteps,
                   restartSteps,
                   restartEnergy,
                   0,
                   temperature,
                   "cmc_log.txt"),
        energyChangePredictor_(
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
          predictorFilename, 
          config, 
          supercellConfig, 
          elementSet
        ),
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
=======
            ceParams,
            config,
            supercellConfig),
>>>>>>> Stashed changes
        atomIndexSelector_(0, config_.GetNumAtoms() - 1)
  {
  }

  pair<size_t, size_t> CanonicalMcAbstract::GenerateLatticeIdJumpPair()
  {
    size_t latticeId1, latticeId2;
    do
    {
      latticeId1 = atomIndexSelector_(generator_);
      latticeId2 = atomIndexSelector_(generator_);
    } while (config_.GetElementOfLattice(latticeId1) == config_.GetElementOfLattice(latticeId2));
    return {latticeId1, latticeId2};
  }

  void CanonicalMcAbstract::Dump() const
  {
    if (is_restarted_)
    {
      is_restarted_ = false;
      return;
    }
    if (steps_ == 0)
    {
      ofs_ << "steps\ttemperature\tenergy\tenergy_per_atom\taverage_energy" << endl;
    }
    if (steps_ % configDumpSteps_ == 0)
    {
      config_.WriteConfig(to_string(steps_) + ".cfg.gz", config_);
    }
    if (steps_ == maximumSteps_)
    {
      config_.WriteConfig("end.cfg.gz", config_);
    }
    unsigned long long int logDumpSteps;
    if (steps_ > 10 * logDumpSteps_)
    {
      logDumpSteps = logDumpSteps_;
    }
    else
    {
      logDumpSteps = static_cast<unsigned long long int>(
          pow(10, static_cast<unsigned long long int>(log10(
                                                          steps_ + 1) -
                                                      1)));
      logDumpSteps = max(logDumpSteps, static_cast<unsigned long long int>(1));
      logDumpSteps = min(logDumpSteps, logDumpSteps_);
    }
    if (steps_ % logDumpSteps == 0)
    {
      ofs_ << steps_ << '\t'
           << temperature_ << '\t'
           << energy_ << '\t'
           << energy_ / static_cast<double>(config_.GetNumAtoms()) << "\t"
           << thermodynamicAveraging_.GetThermodynamicAverage(beta_)
           << endl;
    }
  }
  void CanonicalMcAbstract::SelectEvent(const pair<size_t, size_t> &lattice_id_jump_pair,
                                        const double dE)
  {
    if (dE < 0)
    {
      config_.LatticeJump(lattice_id_jump_pair);
      energy_ += dE;
      absolute_energy_ += dE;

      thermodynamicAveraging_.AddEnergy(dE);
    }
    else
    {
      double possibility = exp(-dE * beta_);
      double random_number = unitDistribution_(generator_);
      if (random_number < possibility)
      {
        config_.LatticeJump(lattice_id_jump_pair);
        energy_ += dE;
        absolute_energy_ += dE;

        thermodynamicAveraging_.AddEnergy(dE);
      }
    }
  }
} // mc