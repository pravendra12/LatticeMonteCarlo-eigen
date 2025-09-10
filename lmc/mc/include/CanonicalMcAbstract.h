/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcAbstract.h
 *  @brief File for CanonicalMcAbstract class declaration.
 */

#ifndef LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#define LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_

#include <random>
#include <mpi.h>
#include <omp.h>
#include <utility>
#include <chrono>
#include "McAbstract.h"
#include "EnergyPredictor.h"

using namespace std;

namespace mc
{
  class CanonicalMcAbstract : public McAbstract
  {
  public:
    /**
     * @brief Constructor for the Canonical Monte Carlo abstract class.
     *
     * Initializes the Canonical Monte Carlo simulation with the given configuration,
     * training setup, and simulation parameters.
     *
     * @param config Configuration used for the simulation.
     * @param supercellConfig Configuration used to train the Cluster Expansion model.
     * @param logDumpSteps Number of steps between logging simulation progress.
     * @param configDumpSteps Number of steps between dumping the configuration.
     * @param maximumSteps Maximum number of Monte Carlo steps to run.
     * @param thermodynamicAveragingSteps Number of steps to use for thermodynamic averaging.
     * @param restartSteps Step index for restarting the simulation.
     * @param restartEnergy Energy value to resume the simulation from.
     * @param temperature Constant simulation temperature (in Kelvin).
     * @param elementSet Set of elements present in the configuration.
     * @param predictorFilename Path to the JSON file containing cluster interaction coefficients.
     */

    CanonicalMcAbstract(Config config,
                        Config supercellConfig,
                        unsigned long long int logDumpSteps,
                        unsigned long long int configDumpSteps,
                        unsigned long long int maximumSteps,
                        unsigned long long int thermodynamicAveragingSteps,
                        unsigned long long int restartSteps,
                        double restartEnergy,
                        double temperature,
                        const ClusterExpansionParameters &ceParams);

    /*!
     * @brief Starts the simulation process.
     */
    void Simulate() override = 0;

  protected:
    // void UpdateTemperature();

    /*!
     * @brief Dumps the current simulation state.
     */
    virtual void Dump() const;

    /*! @brief Generates a lattice Id jump pair.
     *  \returns A pair of lattice Ids.
     */
    pair<size_t, size_t> GenerateLatticeIdJumpPair();

    /*!
     * @brief Selects an event based on the lattice ID pair and energy change.
     *
     * This function determines whether to accept or reject the proposed event,
     * based on Metropolis Algorithm.
     *
     * \param lattice_id_jump_pair Pair of lattice Ids.
     * \param dE                   Energy change associated with the event.
     */
    void SelectEvent(const pair<size_t, size_t> &lattice_id_jump_pair,
                     double dE);

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    // const ClusterExpansionParameters ceParams_;




=======
>>>>>>> Stashed changes
    /** @brief Energy change estimator
     */
    // EnergyPredictor energyChangePredictor_;
    PotentialEnergyEstimator energyChangePredictor_;
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
    /** @brief Energy change estimator
     */
    // EnergyPredictor energyChangePredictor_;
    // PotentialEnergyEstimator energyChangePredictor_;

    EnergyPredictor energyChangePredictor_;
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

    /** @brief Random Lattice Id Generator
     */
    mutable uniform_int_distribution<size_t> atomIndexSelector_;
  };
}

#endif // LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
