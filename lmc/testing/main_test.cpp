#include "Home.h"
#include "Config.h"
#include "ClusterExpansion.h"
#include "SymmetrySpglib.h"
#include "GetEquivalentClusters.h"
#include "SymmetrySpglib.h"
#include "TrainingUtility.h"
#include "ClusterExpansionParameters.h"
#include "LVFEPredictor.h"
#include "KRAPredictor.h"
#include "SymmetricCEPredictor.h"
#include "EnergyPredictor.h"

using namespace std;
using namespace Eigen;
namespace fs = filesystem;

int main()
{
  ClusterExpansionParameters ceParams(
      "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/script/coefficientFile_MoTa_V3.1.json",
      true);

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

  // Seed random generator once
  std::srand(std::time(nullptr));

  // Generate supercell
  auto cfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");
  /*
  cfg.UpdateNeighborList({13});

  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");

  EnergyPredictor symCEEnergyPredictor(
      ceParams,
      cfg,
      primCfg);
  */

  // exit(1);

  cfg.UpdateNeighborList({3, 4, 5});

  // Randomly assign Mo or Ta to each atom
  size_t numAtoms = cfg.GetNumAtoms();
  std::vector<std::string> elements = {"Mo", "Ta"};

  for (size_t i = 0; i < numAtoms; ++i)
  {
    // Pick randomly between "Mo" and "Ta"
    std::string elem = elements[std::rand() % elements.size()];
    cfg.SetElementOfAtom(i, Element(elem));
  }

  cfg.SetElementOfLattice(0, Element("Ta"));
  cfg.SetElementOfLattice(100, Element("Ta"));
  cfg.SetElementOfLattice(111, Element("Ta"));

  cfg.SetElementOfLattice(cfg.GetCentralAtomLatticeId(), Element("X"));

  LVFEPredictor lvfePredictor(ceParams,
                              cfg);

  size_t vacancyId = cfg.GetVacancyLatticeId();
  auto firstNN = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1); // 1NN

  for (auto neighborId : firstNN)
  {
    cfg.SetElementOfLattice(neighborId, Element("Ta")); // change element to Ta
  }

  size_t neighborId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0];
  auto jumpPair = std::make_pair(vacancyId, neighborId);

  // LVFE: Effective vacancy formation energy
  double evfeBefore = lvfePredictor.GetEffectiveVacancyFormationEnergy(cfg, vacancyId);
  cout << "Effective Vacancy Formation Energy (initial): " << evfeBefore << endl;

  // LVFE: Migration energy (forward)
  double deForwardBefore = lvfePredictor.GetDeForVacancyMigration(cfg, jumpPair);
  cout << "Migration Energy (forward, before jump): " << deForwardBefore << endl;

  // Perform the jump
  cfg.LatticeJump(jumpPair);

  // LVFE: Migration energy (backward)
  double deBackwardAfter = lvfePredictor.GetDeForVacancyMigration(cfg, std::make_pair(neighborId, vacancyId));
  cout << "Migration Energy (backward, after jump): " << deBackwardAfter << endl;

  // LVFE: Effective vacancy formation energy at the new site
  double evfeAfter = lvfePredictor.GetEffectiveVacancyFormationEnergy(cfg, neighborId);
  cout << "Effective Vacancy Formation Energy (after jump): " << evfeAfter << endl;

  // KRA predictor
  KRAPredictor kraPredictor(ceParams, cfg);

  // Forward barrier using KRA
  double kRAForward = kraPredictor.GetKRA(cfg, jumpPair);
  cout << "KRA Forward Barrier: " << kRAForward << endl;

  // Backward barrier using KRA
  double kRABackward = kraPredictor.GetKRA(cfg, std::make_pair(neighborId, vacancyId));
  cout << "KRA Backward Barrier: " << kRABackward << endl;

  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");

  // auto cfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");
  cfg.UpdateNeighborList({13});

  SymmetricCEPredictor symCEEnergyPredictor(
      ceParams,
      cfg,
      primCfg);

  EnergyPredictor energyPredictor(
      symCEEnergyPredictor,
      lvfePredictor);

  VacancyMigrationPredictor vacMigrationPredictor(
      kraPredictor,
      energyPredictor);

  cfg.UpdateNeighborList({3, 4, 5});
  

  auto barrierDe = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);

  cout << barrierDe.first << endl;
  cout << barrierDe.second << endl;

  cfg.LatticeJump(jumpPair);

  barrierDe = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);

  cout << barrierDe.first << endl;
  cout << barrierDe.second << endl;

  /// Loop test for the vacMigrationPredictor
  auto elementVac = cfg.GetElementOfLattice(vacancyId);

  vector<size_t> loopingPath;

  cout << elementVac.GetElementString() << endl;

  for (auto nnId1 : cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1))
  {
    auto elem1 = cfg.GetElementOfLattice(nnId1);

    // Iterate over neighbors of nnId1
    for (auto nnId2 : cfg.GetNeighborLatticeIdVectorOfLattice(nnId1, 1))
    {
      if (nnId2 == vacancyId)
        continue; // don’t loop back immediately

      // Neighbors of nnId2
      for (auto nnId3 : cfg.GetNeighborLatticeIdVectorOfLattice(nnId2, 1))
      {
        if (nnId3 == nnId1 || nnId3 == vacancyId)
          continue;

        // Finally, check if it connects back to vacancy
        for (auto nnId4 : cfg.GetNeighborLatticeIdVectorOfLattice(nnId3, 1))
        {
          if (nnId4 == vacancyId)
          {
            // Found a closed loop
            cout << "Loop: "
                 << vacancyId << " -> "
                 << nnId1 << " -> "
                 << nnId2 << " -> "
                 << nnId3 << " -> "
                 << vacancyId << " (all same element)" << endl;

            loopingPath = {vacancyId, nnId1, nnId2, nnId3, vacancyId};

            goto next;
          }
        }
      }
    }
  }

next:
  print1DVector(loopingPath);

  for (int i = 1; i < loopingPath.size() - 1; i++)
  {
    cfg.SetElementOfLattice(loopingPath[i], Element("Ta"));
    cout << loopingPath[i] << endl;
  }

  vector<pair<size_t, size_t>> jumpPairs = {
      {loopingPath[0], loopingPath[1]},
      {loopingPath[1], loopingPath[2]},
      {loopingPath[2], loopingPath[3]},
      {loopingPath[3], loopingPath[0]}};

  double energy = 0;

  double energyFromSwap = 0;

  for (const auto &jumpPair : jumpPairs)
  {
    cout << "Jump pair: (" << jumpPair.first << " -> " << jumpPair.second << ")\n";

    // Barrier and ΔE before jump
    auto barrierDeBefore = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);
    cout << "Before jump: Barrier = " << barrierDeBefore.first
         << ", ΔE = " << barrierDeBefore.second << endl;

    auto dEFromSwap = lvfePredictor.GetDeSwap(cfg, jumpPair);

    cout << "From dESwap: " << dEFromSwap << endl;

    energy += barrierDeBefore.second;
    energyFromSwap += dEFromSwap;

    // Perform the jump
    cfg.LatticeJump(jumpPair);

    // Barrier and ΔE after jump
    auto barrierDeAfter = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);
    cout << "After jump:  Barrier = " << barrierDeAfter.first
         << ", ΔE = " << barrierDeAfter.second << endl;

    cout << "From dESwap: " << lvfePredictor.GetDeSwap(cfg, jumpPair) << endl;

    // Total energy after jump

    cout
        << "---------------------------------------------\n";
  }
  cout << "Total energy after jump: " << energy << endl;

  vacancyId = cfg.GetVacancyLatticeId();
  cout << lvfePredictor.GetEffectiveVacancyFormationEnergy(cfg,
                                                           vacancyId);

  cfg.SetElementOfLattice(vacancyId, Element("Ta"));

  cout << "Ta : " << symCEEnergyPredictor.ComputeLocalFormationEnergyOfSite(cfg, vacancyId) << endl;

  cfg.SetElementOfLattice(vacancyId, Element("Mo"));

  cout << "Mo : " << symCEEnergyPredictor.ComputeLocalFormationEnergyOfSite(cfg, vacancyId) << endl;

  cfg.SetElementOfLattice(vacancyId, Element("X"));

  cout << vacancyId << endl;

  cfg.LatticeJump(jumpPair);

  vacancyId = cfg.GetVacancyLatticeId();

  cfg.SetElementOfLattice(vacancyId, Element("Mo"));

  cout << vacancyId << endl;
  cout << cfg.GetElementOfLattice(vacancyId).GetElementString() << endl;

  cout << cfg.GetElementOfLattice(jumpPair.first).GetElementString() << endl;
  cout << cfg.GetElementOfLattice(jumpPair.second).GetElementString() << endl;

  //
  cout << "From dE swap " << endl;

  auto dEFromSwap = symCEEnergyPredictor.GetDeSwap(
      cfg,
      jumpPair);

  cout << "From dE swap : " << dEFromSwap << endl;

  // const swap

  cout << "From dE swap Const" << endl;

  auto dEFromSwapConst = symCEEnergyPredictor.GetDeSwapConst(
      cfg,
      jumpPair);

  cout << "From dE swap Const: " << dEFromSwapConst << endl;

  cout << cfg.GetElementOfLattice(jumpPair.first).GetElementString() << endl;
  cout << cfg.GetElementOfLattice(jumpPair.second).GetElementString() << endl;

  cfg.SetElementOfLattice(loopingPath[0], Element("X"));
  cout << cfg.GetElementOfLattice(loopingPath[0]).GetElementString() << endl;

  energy = 0;

  energyFromSwap = 0;

  double energyFromEnergyPredictor = 0;


  for (int i = 1; i < loopingPath.size() - 1; i++)
  {
    cfg.SetElementOfLattice(loopingPath[i], Element("Ta"));
    cout << loopingPath[i] << endl;
  }

  for (const auto &jumpPair : jumpPairs)
  {
    cout << "Jump pair: (" << jumpPair.first << " -> " << jumpPair.second << ")\n";

    // Barrier and ΔE before jump
    auto barrierDeBefore = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);
    cout << "Before jump: Barrier = " << barrierDeBefore.first
         << ", ΔE = " << barrierDeBefore.second << endl;

    auto dEFromSwap = lvfePredictor.GetDeSwap(cfg, jumpPair);

    cout << "From dESwap: " << dEFromSwap << endl;

    auto dETrue = energyPredictor.GetEnergyChange(cfg, jumpPair);

    cout << "From EnergyPredictor: " << dETrue << endl;

    energyFromEnergyPredictor += dETrue;

    energy += barrierDeBefore.second;
    energyFromSwap += dEFromSwap;

    // Perform the jump
    cfg.LatticeJump(jumpPair);

    // Barrier and ΔE after jump
    auto barrierDeAfter = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);
    cout << "After jump:  Barrier = " << barrierDeAfter.first
         << ", ΔE = " << barrierDeAfter.second << endl;

    cout << "From dESwap: " << lvfePredictor.GetDeSwap(cfg, jumpPair) << endl;

    dETrue = energyPredictor.GetEnergyChange(cfg, jumpPair);

    cout << "From EnergyPredictor: " << dETrue << endl;

    // Total energy after jump

    cout
        << "---------------------------------------------\n";
  }

  cout << "Total energy after jump: " << energyFromEnergyPredictor << endl;

  //  cfg.UpdateNeighborList({13.1});
  //
  //  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");
  //
  //  EnergyPredictor symCEEnergyPredictor(
  //      ceParams,
  //      cfg,
  //      primCfg);
  //
  // IterateDirectoryToExtractData();
}

/*
int main()
{

  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");
  vector<string> allowedElements = {"Mo", "Ta"};
  vector<double> clusterCutoffs = {13, 8, 6};

  auto supercellCfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");

  supercellCfg.UpdateNeighborList({3});
  auto centralId = supercellCfg.GetCentralAtomLatticeId();

  auto nnId = supercellCfg.GetNeighborLatticeIdVectorOfLattice(centralId, 1)[0];

  // Number of Ta atoms you want to add
  int numTa = 1000; // example

  // Get total number of sites
  int numSites = supercellCfg.GetNumAtoms();

  // Random number generator
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> dis(0, numSites - 1);

  // Keep track of which sites are already replaced
  vector<size_t> replacedSites;

  // Randomly replace Mo with Ta
  for (int i = 0; i < numTa; ++i)
  {
    size_t site;
    do
    {
      site = dis(gen);
    } while (find(replacedSites.begin(), replacedSites.end(), site) != replacedSites.end());

    replacedSites.push_back(site);
    supercellCfg.SetElementOfLattice(site, Element("Ta")); // assuming such a function exists
  }

  vector<double> maxCutoffVector = {
      *max_element(clusterCutoffs.begin(), clusterCutoffs.end())};
  supercellCfg.UpdateNeighborList(maxCutoffVector);

  string predictorFile = "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/icetECIs.json";
  EnergyPredictor symCEEnergyPredictor(
      predictorFile,
      supercellCfg,
      primCfg,
      allowedElements,
      clusterCutoffs);

  auto energyBefore = symCEEnergyPredictor.GetEnergyOfSiteWithVacancy(supercellCfg,
                                                                 centralId);

  auto elvfeBefore = symCEEnergyPredictor.GetEffectiveVacancyFormationEnergy(
      supercellCfg,
      centralId);

  cout << "After JUMP" << endl;
  supercellCfg.SetElementOfLattice(centralId, supercellCfg.GetElementOfLattice(nnId));

  auto energyAfter = symCEEnergyPredictor.GetEnergyOfSiteWithVacancy(supercellCfg,
                                                                nnId);

  auto elvfeAfter = symCEEnergyPredictor.GetEffectiveVacancyFormationEnergy(supercellCfg, nnId);


  cout << "Energy Before : " << energyBefore << endl;

  cout << "Energy After  : " << energyAfter << endl;

  cout << "ELVFE Before  : " << elvfeBefore << endl;

  cout << "ELVFE After   : " << elvfeAfter << endl;
  cout << "Change in Energy : " << energyAfter - energyBefore << endl;
  cout << "Change in ELVFE : " << elvfeAfter - elvfeBefore << endl;

  return 0;
}
  */