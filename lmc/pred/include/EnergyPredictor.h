#ifndef LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

#include "Config.h"
#include "PrintUtility.h"
#include "ClusterExpansionParameters.h"
#include "SymmetricCEPredictor.h"
#include "LVFEPredictor.h"

using namespace std;

class EnergyPredictor
{
public:
  EnergyPredictor(
      SymmetricCEPredictor &symCEEnergyPredictor,
      LVFEPredictor *lvfePredictor = nullptr);

  // Applied for any pair
  double GetEnergyChange(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair);

  // Specifically for <vacancy, atom> pair
  double GetEnergyChangeWithVacancy(
    const Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair);

private:
  SymmetricCEPredictor &symCEEnergyPredictor_;
  LVFEPredictor *lvfePredictor_;

  const vector<string> &allowedElements_;
};

#endif // LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_