/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-01
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 *******************************************************************************/

/*! @file KRAPredictor.cpp
    @brief File contains implementation of Kinetically Resolved Activation Barrier
 */

#include "KRAPredictor.h"

KRAPredictor::KRAPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &config) : maxBondOrder_(ceParams.GetMaxBondOrder("kra")),
                            maxClusterSize_(
                                ceParams.GetMaxClusterSize("kra")),
                            referenceJumpDirection_(
                                ceParams.GetReferenceJumpDirection()),
                            atomicBasis_(
                                ceParams.GetElementSet("kra"),
                                ceParams.GetBasisType()),
                            canonicalReferenceMap_(
                                GetCenteredNeighborsAlongJumpDirection(
                                    config,
                                    maxBondOrder_,
                                    referenceJumpDirection_)),
                            symmetryOperations_(
                                GetSymmetryOperations(
                                    config)),
                            encodedOrbitsForPair_(
                                GetLocalEncodedOrbitsForPair(
                                    config,
                                    maxBondOrder_,
                                    maxClusterSize_,
                                    referenceJumpDirection_,
                                    true)),
                            KECIsMap_(
                                ceParams.GetKECIsMap())

{
  // Precompute the canonical sorted lattice Ids
  // And PairToIndexMap
  GetSymmetricallySortedLatticeIdsVectorMap(config);

  // Print KRAPredictorInfo
  PrintKRAPredictorInfo();
}

double KRAPredictor::GetKRA(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  const auto &elementFirst = config.GetElementOfLattice(latticeIdJumpPair.first);
  const auto &elementSecond = config.GetElementOfLattice(latticeIdJumpPair.second);

  // At least one element must be a vacancy
  if (elementFirst != Element("X") && elementSecond != Element("X"))
  {
    throw runtime_error("Error in `KRAPredictor::GetKRA`: Jump pair must contain at least one vacancy.");
  }

  Element migratingElement = (elementFirst == Element("X"))
                                 ? elementSecond
                                 : elementFirst;

  // Canonical Pair
  size_t id1 = latticeIdJumpPair.first;
  size_t id2 = latticeIdJumpPair.second;

  pair<size_t, size_t> canonicalJumpPair = (id1 < id2)
                                               ? make_pair(id1, id2)
                                               : make_pair(id2, id1);

  /*
  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForPair(
      config,
      canonicalJumpPair,
      maxBondOrder_,
      canonicalReferenceMap_,
      symmetryOperations_);
  */

  size_t indexOfPair = pairToIndexHashMap_.at(canonicalJumpPair);
  auto canonicalSortedLatticeIds = symmetricallySortedLatticeIdsVectorMap_[indexOfPair];
  
  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis_,
      canonicalSortedLatticeIds,
      encodedOrbitsForPair_);

  // EKRA = J.Φ_α

  VectorXd kecis = KECIsMap_.at(migratingElement);

  if (kecis.size() != correlationVector.size())
  {
    throw runtime_error(
        "Error in `KRAPredictor::GetKRA`: migratingElement `" + migratingElement.GetElementString() +
        "` has kecis size (" + to_string(kecis.size()) +
        ") which does not match correlationVector size (" + to_string(correlationVector.size()) +
        ").");
  }

  double eKraValue = kecis.dot(correlationVector);

  return eKraValue;
}

void KRAPredictor::GetSymmetricallySortedLatticeIdsVectorMap(
    const Config &config)
{
  const size_t numLattices = config.GetNumLattices();

  auto numNeighboursPerSite = config.GetNeighborLatticeIdVectorOfLattice(0, 1).size();
  const size_t numUniquePairs = (numLattices * numNeighboursPerSite) / 2;

  // Declared as private members
  symmetricallySortedLatticeIdsVectorMap_.reserve(numUniquePairs);
  pairToIndexHashMap_.reserve(numUniquePairs);

  size_t numSortedSites = canonicalReferenceMap_.size();
  vector<size_t> canonicalSortedLatticeIds;
  canonicalSortedLatticeIds.reserve(numSortedSites);

  size_t index = 0;

  // Iterate over all the sites and get store the sorted lattice Ids
  for (size_t id1 = 0; id1 < numLattices; id1++)
  {
    for (const size_t id2 : config.GetNeighborLatticeIdVectorOfLattice(id1, 1))
    {
      pair<size_t, size_t> canonicalJumpPair = (id1 < id2)
                                                   ? make_pair(id1, id2)
                                                   : make_pair(id2, id1);

      if (pairToIndexHashMap_.find(canonicalJumpPair) != pairToIndexHashMap_.end())
      {
        // Pair already exists
        continue;
      }

      pairToIndexHashMap_[canonicalJumpPair] = index;
      index += 1;

      canonicalSortedLatticeIds = GetCanonicalSortedSitesForPair(
          config,
          canonicalJumpPair,
          maxBondOrder_,
          canonicalReferenceMap_,
          symmetryOperations_);

      symmetricallySortedLatticeIdsVectorMap_.emplace_back(move(canonicalSortedLatticeIds));
    }
  }
}

void KRAPredictor::PrintKRAPredictorInfo() const
{
  const int width = 80;
  const int labelWidth = 40;
  const int valueWidth = width - labelWidth - 2; // 2 for spacing

  // Header
  cout << string(width, '-') << "\n";
  cout << setw((width + 22) / 2) << right << "KRA Predictor Info" << "\n";
  cout << string(width, '-') << "\n";

  // Table of values
  cout << left << setw(labelWidth) << "Max bond order:"
       << right << setw(valueWidth) << maxBondOrder_ << "\n";

  cout << left << setw(labelWidth) << "Max cluster size:"
       << right << setw(valueWidth) << maxClusterSize_ << "\n";

  cout << left << setw(labelWidth) << "Reference jump direction:"
       << right << setw(valueWidth - 3)
       << "(" << referenceJumpDirection_.transpose() << ")" << "\n";

  cout << left << setw(labelWidth) << "Canonical reference map size:"
       << right << setw(valueWidth) << canonicalReferenceMap_.size() << "\n";

  cout << left << setw(labelWidth) << "Number of symmetry operations:"
       << right << setw(valueWidth) << symmetryOperations_.size() << "\n";

  cout << left << setw(labelWidth) << "Encoded orbits for pair size:"
       << right << setw(valueWidth) << encodedOrbitsForPair_.size() << "\n";

  cout << left << setw(labelWidth) << "Size of precomputed vector map for jump pairs:"
       << right << setw(valueWidth-10) << "(" << pairToIndexHashMap_.size() << ", "
       << symmetricallySortedLatticeIdsVectorMap_.size() << ")" << "\n";

  cout << left << setw(labelWidth) << "Kinetic ECIs map size:"
       << right << setw(valueWidth) << KECIsMap_.size() << "\n";

  cout << left << setw(labelWidth) << "KECIs per element:" << "\n";
  for (const auto &[element, vec] : KECIsMap_)
  {
    cout << "  " << left << setw(labelWidth - 2)
         << element.GetElementString()
         << right << setw(valueWidth)
         << vec.size() << "\n";
  }

  cout << string(width, '-') << "\n\n";
}