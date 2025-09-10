#include "EnergyPredictor.h"

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
#include "PrintUtility.h"

// This should just take the config and predictorFilename and
// Rest should be read from the param file

EnergyPredictor::EnergyPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &config) : maxClusterSize_(ceParams.GetMaxClusterSizeCE()),
                            maxBondOrder_(ceParams.GetMaxBondOrderCE()),
                            ecis_(ceParams.GetECIs()),
                            atomicBasis_(
                                ceParams.GetElementSetCE(),
                                ceParams.GetBasisType()),
                            equivalentClustersEncoding_(
                                GetEquivalentClustersEncoding(
                                    config,
                                    maxBondOrder_,
                                    maxClusterSize_))
{
  auto atomVector = config.GetAtomVector();
  set<Element> configElementSet(atomVector.begin(), atomVector.end());

  auto elementSet = ceParams.GetElementSetCE();

  for (auto ele : configElementSet)
  {
    if (elementSet.find(ele) == elementSet.end())
    {
      cerr << "Error in `EnergyPredictor`: element " << ele
           << " is not in the allowed element set!\n";

      cerr << "This cluster expansion is defined only for the following elements: ";
      for (const auto &allowedEle : elementSet)
        cerr << allowedEle << " ";
      cerr << endl;

      exit(1);
    }
  }
}

double EnergyPredictor::ComputeEnergyOfSite(
=======
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
EnergyPredictor::EnergyPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &supercellConfig,
    const Config &primitiveConfig) : nAtoms_(supercellConfig.GetNumAtoms()),
                                     allowedElements_(
                                         ceParams.GetAllowedElements()),
                                     clusterCutoffs_(
                                         ceParams.GetClusterCutoffs()),
                                     symCE_(supercellConfig,
                                            primitiveConfig,
                                            allowedElements_,
                                            clusterCutoffs_),
                                     localOrbitsEncoding_(
                                         symCE_.GetLocalOrbitsEncoding()),
                                     chemicalPotentialsMap_(
                                         ceParams.GetChemicalPotentialsMap()),
                                     elementCountMap_(
                                         GetElementCountMap(
                                             supercellConfig)),
                                     ecis_(
                                         ceParams.GetECIs("symCE"))

{
  const int width = 80;
  const int labelWidth = 40;
  const int valueWidth = width - labelWidth - 2;

  // Header
  cout << string(width, '-') << "\n";
  cout << setw((width + 22) / 2) << right << "Energy Predictor Info" << "\n";
  cout << string(width, '-') << "\n";

  // Table
  cout << left << setw(labelWidth) << "Number of atoms:"
            << right << setw(valueWidth) << nAtoms_ << "\n";

  cout << left << setw(labelWidth) << "Allowed elements:"
            << right << setw(valueWidth);
  for (const auto &ele : allowedElements_)
    cout << ele << " ";
  cout << "\n";

  cout << left << setw(labelWidth) << "Cluster cutoffs:"
            << right << setw(valueWidth);
  for (const auto &cut : clusterCutoffs_)
    cout << cut << " ";
  cout << "\n";

  cout << left << setw(labelWidth) << "Local orbits encoding size:"
            << right << setw(valueWidth) << localOrbitsEncoding_.size() << "\n";

  cout << left << setw(labelWidth) << "Chemical potentials:"
            << "\n";
  for (const auto &[ele, mu] : chemicalPotentialsMap_)
  {
    cout << "  " << left << setw(labelWidth - 2) << ele
              << right << setw(valueWidth) << mu << "\n";
  }

  cout << left << setw(labelWidth) << "Sym CE ECIs size:"
            << right << setw(valueWidth) << ecis_.size() << "\n";

  cout << string(width, '-') << "\n\n";

}

// Returns total formation energy
/*
  Formation energy per atom:

      Ef_per_atom = (E_total - (N_A * μ_A + N_B * μ_B)) / (N_A + N_B)

  - E_total: total energy of the supercell
  - N_A, N_B: number of atoms of type A and B
  - μ_A, μ_B: reference chemical potentials of elements A and B

  The cluster expansion model predicts Ef_per_atom,
  but for most calculations (e.g., energy differences or KMC jumps),
  we are interested in the total energy E_total.
*/

double EnergyPredictor::GetTotalFormationEnergy(
    const vector<double> &clusterVector)
{

  double formationEnergy = 0; // per atom

  // print1DVector(ecis_);

  // cout << ecis_.size() << endl;
  // cout << clusterVector.size() << endl;

  for (int i = 0; i < clusterVector.size(); i++)
  {
    formationEnergy += clusterVector[i] * ecis_[i];
  }

  double totalFormationEnergy = nAtoms_ * formationEnergy;

  return totalFormationEnergy;
}

double EnergyPredictor::GetTotalEnergy(
    const vector<double> &clusterVector)
{
  auto totalFormationEnergy = GetTotalFormationEnergy(
      clusterVector);

  double muContribution = 0;
  for (const auto &[element, count] : elementCountMap_)
  {
    // Chemical Potential Value
    double muVal = chemicalPotentialsMap_.at(element);

    muContribution += count * muVal;
  }

  // Etotal = (N_A + N_B)*Eformation + (N_A * μ_A + N_B * μ_B))
  double totalEnergy = totalFormationEnergy + muContribution;

  return totalEnergy;
}

// This will value will be higher than the formation energy for whole
// supercell due to the constant term contribution to be 1
// But that does not matter as we only dE which can be calculated using
// total formation energy change which is same as dE

double EnergyPredictor::ComputeLocalFormationEnergyOfSite(
>>>>>>> Stashed changes
    const Config &config,
    const size_t &latticeId)
{
  // Can be Stored for each site or may be not
  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
      config,
      latticeId,
      maxBondOrder_);

  // Need to include the site as well to be consistent with the encoding
  canonicalSortedLatticeIds.emplace_back(latticeId);

  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis_,
      canonicalSortedLatticeIds,
      equivalentClustersEncoding_);

  // E = J.Φ_α
  double energyValue = ecis_.dot(correlationVector) / config.GetNumAtoms();

  return energyValue;
}

double EnergyPredictor::ComputeEnergyOfConfig(
    const Config &config)
{
  double totalEnergy = 0;

  for (size_t latticeId = 0; latticeId < config.GetNumLattices(); latticeId++)
  {
    totalEnergy += ComputeEnergyOfSite(config, latticeId);
  }

  return totalEnergy;
}

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
double EnergyPredictor::GetDeMigration(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  auto canonicalSortedLatticeIdsSite1 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdJumpPair.first,
      maxBondOrder_);

  auto canonicalSortedLatticeIdsSite2 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdJumpPair.second,
      maxBondOrder_);

  cout << "Site1 : ";
  cout << latticeIdJumpPair.first << endl;
  print1DVector(canonicalSortedLatticeIdsSite1);
  cout << "Site2 : ";
  cout << latticeIdJumpPair.second << endl;

  print1DVector(canonicalSortedLatticeIdsSite2);

  return 0;
}

=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
double EnergyPredictor::GetDeSwap(
    Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  // Before Swap

  auto nnLatticeIds = config.GetNeighboringLatticeIdSetOfPair(latticeIdJumpPair, maxBondOrder_);
  nnLatticeIds.insert(latticeIdJumpPair.first);
  nnLatticeIds.insert(latticeIdJumpPair.second);

  double energyBeforeSwap = 0;

  for (const auto latticeId : nnLatticeIds)
  {
    energyBeforeSwap += ComputeEnergyOfSite(config, latticeId);
  }

  config.LatticeJump(latticeIdJumpPair);

  // After Swap
  double energyAfterSwap = 0;

  for (const auto latticeId : nnLatticeIds)
  {
    energyAfterSwap += ComputeEnergyOfSite(config, latticeId);
  }
  
  auto dE = energyAfterSwap - energyBeforeSwap;

  config.LatticeJump(latticeIdJumpPair);

  return dE;
}

unordered_map<string, int> EnergyPredictor::GetElementCountMap(
    const Config &supercellConfig)
{
  unordered_map<string, int> elementCountMap;
  for (auto ele : supercellConfig.GetAtomVector())
  {
    elementCountMap[ele.GetElementString()]++;
  }

  return elementCountMap;
}
