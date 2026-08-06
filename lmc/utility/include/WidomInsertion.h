#ifndef LMC_UTILITY_INCLUDE_WIDOMINSERTION_H_
#define LMC_UTILITY_INCLUDE_WIDOMINSERTION_H_

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "Config.h"
#include "SymmetricCEPredictor.h"

namespace fs = std::filesystem;
using namespace std;

/*

// config : can be equilibrium config or random config depends on what quantiy want to be computed
// but in priciple should be equilibuirm config
// symCEPredictor : To compute local formation energy of insertion

*/
struct WidomInsertionRecord
{
  size_t latticeId;
  Element originalElement;
  double localFormationEnergyChange;
  vector<string> firstNNElements;
  vector<string> secondNNElements;
  vector<string> thirdNNElements;
};

/**
 * @brief Perform Widom-style ghost insertion over all lattice sites.
 *
 * For each lattice site i:
 *   1) Save the original element at site i
 *   2) Temporarily set site i to `ghostElement`
 *   3) Compute the local formation energy for site i using `symCEPredictor`
 *   4) Restore the original element at site i
 *
 * This routine does NOT accumulate Boltzmann weights or averages by itself;
 * it only evaluates per-site local formation energies for a ghost substitution.
 *
 * @param config
 *   Configuration to probe. It is modified in-place during evaluation but is
 *   restored to its original state before moving to the next site.
 *
 * @param symCEPredictor
 *   Cluster-expansion-based predictor used to compute the local formation energy.
 *
 * @param ghostElement
 *   The element/species used as the ghost insertion (e.g., vacancy "X").
 *
 * @return std::vector<WidomInsertionRecord>
 *   Per-site Widom insertion results, in the same order as the lattice sweep
 *   (`latticeId = 0 .. numSites-1`). Each record contains:
 *     - `latticeId`
 *     - `originalElement`
 *     - `localFormationEnergy`
 */

vector<WidomInsertionRecord> WidomInsertion(
    Config &config,
    SymmetricCEPredictor &symCEPredictor,
    Element ghostElement = Element("X"));

// This is a helper funciton which iterates over all the files after CMC simulation
void WidomInsertionHelper(
    unsigned long long int initialSteps,
    unsigned long long int incrementSteps,
    SymmetricCEPredictor &symCEPredictor,
    Element ghostElement = Element("X"));

#endif // LMC_UTILITY_INCLUDE_WIDOMINSERTION_H_
