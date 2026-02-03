#include "SymmetricCEPredictor.h"

SymmetricCEPredictor::SymmetricCEPredictor(
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
                                         ceParams.GetECIs("symCE")),
                                     symmetricallySortedLatticeIdsVectorMap_(
                                         move(GetSymmetricallySortedLatticeIdsVectorMap(
                                             supercellConfig)))

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

  cout << left << setw(labelWidth) << "Cluster vector size (with constant term):"
       << right << setw(valueWidth) << (localOrbitsEncoding_.size() + 1) << "\n";

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

const vector<string> &SymmetricCEPredictor::GetAllowedElements() const
{
  return allowedElements_;
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

double SymmetricCEPredictor::GetTotalFormationEnergy(
    const vector<double> &clusterVector)
{

  double formationEnergy = 0; // per atom

  for (int i = 0; i < clusterVector.size(); i++)
  {
    formationEnergy += clusterVector[i] * ecis_[i];
  }

  double totalFormationEnergy = nAtoms_ * formationEnergy;

  return totalFormationEnergy;
}

double SymmetricCEPredictor::GetTotalEnergy(
    const vector<double> &clusterVector)
{
  auto totalFormationEnergy = GetTotalFormationEnergy(
      clusterVector);

  double muContribution = 0;
  for (const auto &[element, count] : elementCountMap_)
  {
    // Chemical Potential Value
    auto muVal = chemicalPotentialsMap_.at(element);

    muContribution += count * muVal;
  }

  // Etotal = (N_A + N_B)*Eformation + (N_A * μ_A + N_B * μ_B))
  double totalEnergy = totalFormationEnergy + muContribution;

  return totalEnergy;
}

double SymmetricCEPredictor::GetTotalEnergy(
    const vector<double> &clusterVector,
    const unordered_map<string, int> &elementCountMap)
{
  auto totalFormationEnergy = GetTotalFormationEnergy(
      clusterVector);

  double muContribution = 0;
  for (const auto &[element, count] : elementCountMap)
  {
    // Chemical Potential Value
    auto muVal = chemicalPotentialsMap_.at(element);

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

double SymmetricCEPredictor::ComputeLocalFormationEnergyOfSite(
    const Config &config,
    const size_t &latticeId)
{
  /*
  // Can be Stored for each site or may be not
  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
      config,
      latticeId,
      1);

  // Need to include the site as well to be consistent with the encoding
  canonicalSortedLatticeIds.emplace_back(latticeId);
  */

  // symmetricallySortedLatticeIdsVectorMap_ is indexed using latticeId
  auto canonicalSortedLatticeIds = symmetricallySortedLatticeIdsVectorMap_[latticeId];

  auto clusterVector = symCE_.GetLocalClusterVector(
      config,
      canonicalSortedLatticeIds,
      localOrbitsEncoding_);

  // cout << clusterVector.size() << endl;

  // E = J.Φ_α
  double energyValue = GetTotalFormationEnergy(clusterVector);

  return energyValue;
}

double SymmetricCEPredictor::ComputeEnergyOfConfig(
    const Config &config)
{
  auto clusterVector = symCE_.GetClusterVectorForConfig(config);

  double totalEnergy = GetTotalEnergy(clusterVector);

  return totalEnergy;
}

double SymmetricCEPredictor::GetDeSwap(
    Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  // Before Swap

  auto energyBeforeSwap = ComputeLocalFormationEnergyOfSite(
                              config,
                              latticeIdJumpPair.first) +
                          ComputeLocalFormationEnergyOfSite(
                              config,
                              latticeIdJumpPair.second);

  config.LatticeJump(latticeIdJumpPair);

  // After Swap
  double energyAfterSwap = ComputeLocalFormationEnergyOfSite(
                               config,
                               latticeIdJumpPair.first) +
                           ComputeLocalFormationEnergyOfSite(
                               config,
                               latticeIdJumpPair.second);

  auto dE = energyAfterSwap - energyBeforeSwap;

  config.LatticeJump(latticeIdJumpPair);

  return dE;
}

double SymmetricCEPredictor::ComputeLocalFormationEnergyForPair(
    const Config &config,
    const pair<size_t, size_t> &latticeIdPair, // here the elements will be assigned in the same order
    // This is expected to get desired behaviour
    // latticeIdPairElements.first <-> latticeIdPair.first
    // latticeIdPairElements.second <-> latticeIdPair.second
    const pair<Element, Element> &latticeIdPairElements)
{

  // 1 AT THE END OF VARIABLE CORRESPONDS TO latticeIdPair.first
  // 2 AT THE END OF VARIABLE CORRESPONDS TO latticeIdPair.second

  // Sorted Lattice Ids
  // Can be Stored for each site or may be not
  /*
  auto canonicalSortedLatticeIds1 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdPair.first,
      1);
  canonicalSortedLatticeIds1.emplace_back(latticeIdPair.first);
  */

  // symmetricallySortedLatticeIdsVectorMap_ is indexed using latticeId
  auto canonicalSortedLatticeIds1 = symmetricallySortedLatticeIdsVectorMap_[latticeIdPair.first];

  /*
  auto canonicalSortedLatticeIds2 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdPair.second,
      1);
  canonicalSortedLatticeIds2.emplace_back(latticeIdPair.second);
  */

  auto canonicalSortedLatticeIds2 = symmetricallySortedLatticeIdsVectorMap_[latticeIdPair.second];

  auto clusterVectorPair = symCE_.GetLocalClusterVectorForPair(
      config,
      latticeIdPair,
      latticeIdPairElements,
      canonicalSortedLatticeIds1,
      canonicalSortedLatticeIds2,
      localOrbitsEncoding_);


  if (ecis_.size() != clusterVectorPair.first.size() ||
      ecis_.size() != clusterVectorPair.second.size())
  {
    ostringstream msg;
    msg << "Error in `SymmetricCEPredictor::ComputeLocalFormationEnergyForPair`: "
        << "ECI size (" << ecis_.size() << ") does not match "
        << "clusterVectorPair.first size (" << clusterVectorPair.first.size() << ") or "
        << "clusterVectorPair.second size (" << clusterVectorPair.second.size() << ").";
    throw runtime_error(msg.str());
  }

  // E = J.Φ_α
  double formationEnergy = 0;

  for (int i = 0; i < ecis_.size(); i++)
  {
    // Contribution of first site
    formationEnergy += clusterVectorPair.first[i] * ecis_[i];
    // Contribution of second site
    formationEnergy += clusterVectorPair.second[i] * ecis_[i];
  }

  double totalFormationEnergy = nAtoms_ * formationEnergy;

  return totalFormationEnergy;
}

// Const Object
double SymmetricCEPredictor::GetDeSwapConst(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  auto element1 = config.GetElementOfLattice(latticeIdJumpPair.first);
  auto element2 = config.GetElementOfLattice(latticeIdJumpPair.second);

  // Before Swap
  // element1 <-> latticeIdJumpPair.first
  // element2 <-> latticeIdJumpPair.second

  pair<Element, Element> latticeIdPairElements = {element1, element2};

  // Function Swaps it as
  // element1 <-> latticeIdJumpPair.second
  // element2 <-> latticeIdJumpPair.first

  double dEValue = GetDeSwapConst(
      config,
      latticeIdJumpPair,
      latticeIdPairElements);

  return dEValue;
}

double SymmetricCEPredictor::GetDeSwapConst(
    const Config &config,
    const pair<size_t, size_t> &latticeIdPair,
    const pair<Element, Element> &latticeIdPairElements)
{
  // Before Swap
  // element1 <-> latticeIdPair.first
  // element2 <-> latticeIdPair.second

  auto energyBeforeSwap = ComputeLocalFormationEnergyForPair(
      config,
      latticeIdPair,
      latticeIdPairElements);

  // After Swap
  // element1 will move to latticeIdPair.second
  // element2 will move to latticeIdPair.first

  pair<Element, Element> afterSwapPairElements = {
      latticeIdPairElements.second,
      latticeIdPairElements.first};

  // element1 <
  double energyAfterSwap = ComputeLocalFormationEnergyForPair(
      config,
      latticeIdPair,
      afterSwapPairElements);

  auto dE = energyAfterSwap - energyBeforeSwap;

  return dE;
}

unordered_map<string, int> SymmetricCEPredictor::GetElementCountMap(
    const Config &supercellConfig)
{
  unordered_map<string, int> elementCountMap;
  for (auto ele : supercellConfig.GetAtomVector())
  {
    elementCountMap[ele.GetElementString()]++;
  }

  return elementCountMap;
}

vector<vector<size_t>> SymmetricCEPredictor::GetSymmetricallySortedLatticeIdsVectorMap(
    const Config &supercellConfig)
{
  // Expecting supercellConfig to be updated in the following way
  // supercellConfig.UpdateNeighbours({maxClusterCutoff})
  const size_t maxBondOrder = 1;

  const size_t numLattices = supercellConfig.GetNumLattices();

  vector<vector<size_t>> symmetricallySortedLatticeIdsVectorMap;
  symmetricallySortedLatticeIdsVectorMap.reserve(numLattices);

  size_t numSortedIds = supercellConfig.GetNeighborLatticeIdVectorOfLattice(0, 1).size() + 1;
  vector<size_t> canonicalSortedLatticeIds;
  canonicalSortedLatticeIds.reserve(numSortedIds);

  // Iterate over all the sites and get store the sorted lattice Ids
  for (size_t latticeId = 0; latticeId < numLattices; latticeId++)
  {
    canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
        supercellConfig,
        latticeId,
        maxBondOrder);

    canonicalSortedLatticeIds.emplace_back(latticeId);

    symmetricallySortedLatticeIdsVectorMap.emplace_back(move(canonicalSortedLatticeIds));
  }

  return symmetricallySortedLatticeIdsVectorMap;
}
