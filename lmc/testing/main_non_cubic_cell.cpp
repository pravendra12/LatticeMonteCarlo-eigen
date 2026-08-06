/*! \file  main.cpp
 *  \brief File for the main function.
 */

#include "Home.h"
#include "Home.h"
#include <iostream>
#include <vector>
#include <string>
#include <cstddef>
#include <random>

using namespace std;

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>

using namespace std;

int main()
{
  // ============================================================
  // Input parameters
  // ============================================================

  const filesystem::path inputConfigPath =
      "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/bin/MoTa_10x10x10.cfg.gz";

  const filesystem::path coefficientFilePath =
      "/home/pravendra3/Documents/coefficients/"
      "coefficientFile_MoTa_V3.2.json";

  const filesystem::path outputDirectory =
      "testingoC12/MoTa_10x10x10_cubic";

  const filesystem::path jumpParameterOutputPath =
      outputDirectory / "i.testing_v2_oC12";

  const filesystem::path energyOutputPath =
      outputDirectory / "o.testing_v2_oC12_lmc";

  const string outputConfigPrefix =
      "MoTa_10x10x10_cubic";

  constexpr size_t primitiveSupercellSize = 2;
  constexpr double latticeParameter = 3.2;
  constexpr int numberOfJumpAttempts = 1000;

  // ============================================================
  // Create output directory
  // ============================================================

  filesystem::create_directories(outputDirectory);

  // ============================================================
  // Read configuration
  // ============================================================

  // auto cfg = Config::ReadCfg(inputConfigPath.string());

  auto cfg = Config::GenerateAlloySupercell(10, 3.2, "BCC", {"Mo", "Ta"}, {50, 50}, 1);

  cout << cfg.GetBasis() << endl;
  cout << cfg.GetNumAtoms() << endl;

  auto centralLatticeId = cfg.GetCentralAtomLatticeId();

  cout << centralLatticeId << endl;
  cout << cfg.GetRelativePositionOfLattice(centralLatticeId).transpose() << endl;
  cout << cfg.GetCartesianPositionOfLattice(centralLatticeId).transpose() << endl;


  // ============================================================
  // Initialize cluster expansion
  // ============================================================

  ClusterExpansionParameters ceParams(
      coefficientFilePath.string());

  const double maxClusterCutoff =
      ceParams.GetMaxClusterCutoff();

  cfg.UpdateNeighborList({maxClusterCutoff});

  cout << "Done with neighbour list update" << endl;

  auto primConfig = Config::GenerateSupercell(
      primitiveSupercellSize,
      latticeParameter,
      "X",
      "BCC");

  SymmetricCEPredictor symCEEnergyPredictor(
      ceParams,
      cfg,
      primConfig);

  cout << symCEEnergyPredictor.ComputeEnergyOfConfig(cfg)
       << endl;

  EnergyPredictor energyPredictor(
      symCEEnergyPredictor);

  // ============================================================
  // Open output files
  // ============================================================

  ofstream jumpParameterFile(
      jumpParameterOutputPath);

  if (!jumpParameterFile)
  {
    cerr << "Failed to open: "
         << jumpParameterOutputPath << endl;

    return 1;
  }

  ofstream energyOutputFile(
      energyOutputPath);

  if (!energyOutputFile)
  {
    cerr << "Failed to open: "
         << energyOutputPath << endl;

    return 1;
  }

  jumpParameterFile
      << "# jump_number lattice_id_1 lattice_id_2\n";

  energyOutputFile << setprecision(17);

  // ============================================================
  // Random-number generator
  // ============================================================

  random_device randomDevice;
  mt19937_64 randomGenerator(
      randomDevice());

  uniform_int_distribution<size_t> latticeDistribution(
      0,
      cfg.GetNumLattices() - 1);

  // ============================================================
  // Perform jump attempts
  // ============================================================

  for (int i = 0; i < numberOfJumpAttempts; ++i)
  {
    const size_t latticeId1 =
        latticeDistribution(randomGenerator);

    const size_t latticeId2 =
        latticeDistribution(randomGenerator);

    if (latticeId1 == latticeId2)
    {
      continue;
    }

    if (cfg.GetElementOfLattice(latticeId1) !=
        cfg.GetElementOfLattice(latticeId2))
    {
      const pair<size_t, size_t> jumpPair = {
          latticeId1,
          latticeId2};

      // Write the jump parameters before performing the jump.
      jumpParameterFile
          << i << ' '
          << latticeId1 << ' '
          << latticeId2 << '\n';

      const double dE =
          energyPredictor.GetEnergyChange(
              cfg,
              jumpPair);

      cout << "---------- "
           << i
           << " ----------------"
           << endl;

      cout << setprecision(17)
           << "dE: "
           << dE
           << endl;

      energyOutputFile
          << "------- "
          << i
          << " --------\n";

      energyOutputFile
          << dE
          << "\n\n";

      const filesystem::path initialConfigPath =
          outputDirectory /
          (outputConfigPrefix +
           "_0_jump_" +
           to_string(i) +
           ".cfg.gz");

      const filesystem::path finalConfigPath =
          outputDirectory /
          (outputConfigPrefix +
           "_1_jump_" +
           to_string(i) +
           ".cfg.gz");

      Config::WriteConfig(
          initialConfigPath.string(),
          cfg);

      cfg.LatticeJump(jumpPair);

      Config::WriteConfig(
          finalConfigPath.string(),
          cfg);

      cout << endl;
      cout << endl;
    }
  }

  cout << "Jump parameters written to: "
       << jumpParameterOutputPath
       << endl;

  cout << "Energy changes written to: "
       << energyOutputPath
       << endl;

  return 0;
}

/*
int main(int argc, char *argv[])
{

 if (argc == 1)
 {
   std::cout << "No input parameter filename." << std::endl;
   return 1;
 }
 api::Parameter parameter(argc, argv);
 api::Print(parameter);
 api::Run(parameter);
}*/
