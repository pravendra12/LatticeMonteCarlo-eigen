/*! \file  main.cpp
 *  \brief File for the main function.
 */

#include "Config.h"
#include "KRAPredictor.h"
#include "ClusterExpansionParameters.h"
#include "VacancyMigrationPredictor.h"
#include "LVFEPredictor.h"
#include "EnergyPredictor.h"

int main()
{
  auto cfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");

  // Randomly assign Mo or Ta to each atom
  size_t numAtoms = cfg.GetNumAtoms();
  std::vector<std::string> elements = {"Mo", "Ta"};

  for (size_t i = 0; i < numAtoms; ++i)
  {
    // Pick randomly between "Mo" and "Ta"
    std::string elem = elements[std::rand() % elements.size()];
    cfg.SetElementOfAtom(i, Element(elem));
  }

  ClusterExpansionParameters ceParams("/home/pravendra3/Documents/LatticeMonteCarlo-eigen/script/coefficientFile_MoTa_V3.1.json");
  cfg.UpdateNeighborList({ceParams.GetMaxClusterCutoff()});

  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");


  SymmetricCEPredictor symCEEnergyPredictor(
      ceParams,
      cfg,
      primCfg);
  cfg.UpdateNeighborList({3, 4, 5});

  pair<size_t, size_t> jumpPair = {0, cfg.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]};

  cfg.SetElementOfLattice(0, Element("X"));
  // cfg.SetElementOfLattice(jumpPair.first, Element("Mo"));
  cfg.SetElementOfLattice(jumpPair.second, Element("Ta"));

  KRAPredictor kraPredictor(
      ceParams,
      cfg);

  LVFEPredictor lvfePredictor(
      ceParams,
      cfg);

  EnergyPredictor energyPredictor(
      symCEEnergyPredictor,
      lvfePredictor);

  VacancyMigrationPredictor vacMigPredictor(
      kraPredictor,
      energyPredictor);

  auto start = std::chrono::high_resolution_clock::now();
  double result = energyPredictor.GetEnergyChange(cfg, jumpPair);
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsedEnergy = end - start;
  std::cout << "GetEnergyChange = " << elapsedEnergy.count() << " ms\n";
  cout << result << endl;

  start = std::chrono::high_resolution_clock::now();
  auto barrierDe = vacMigPredictor.GetBarrierAndDeltaE(cfg, jumpPair);
  end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsedBarrier = end - start;
  std::cout << "GetBarrierAndDeltaE = " << elapsedBarrier.count() << " ms\n";
  cout << barrierDe.first << endl;
  cout << barrierDe.second << endl;
}
