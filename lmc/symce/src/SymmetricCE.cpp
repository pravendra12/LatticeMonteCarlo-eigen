#include "SymmetricCE.h"

SymmetricCE::SymmetricCE(
    const Config &supercellConfig, // No need to update the neighbour list
    const Config &primitiveConfig,
    const vector<string> &allowedElements,
    const vector<double> &clusterCutoffs,
    const double symprec) : supercellConfig_(supercellConfig),
                            symCEInterface_(primitiveConfig,
                                            allowedElements,
                                            clusterCutoffs,
                                            symprec),
                            fractionalPositionTolerance_(symCEInterface_.GetFractionalPositionTolerance()),
                            primitiveOrbitList_(
                                make_shared<OrbitList>(OrbitList(
                                    symCEInterface_.GetPrimitiveStructure(),
                                    symCEInterface_.GetMatrixOfEquivalentSites(),
                                    symCEInterface_.GetNeighbourLists(),
                                    fractionalPositionTolerance_))),
                            supercellStructure_(make_shared<Structure>(
                                ConvertConfigToStructure(
                                    supercellConfig)))

{

  primitiveOrbitList_->removeOrbitsWithInactiveSites();

  clusterSpace_ = ClusterSpace(
      primitiveOrbitList_,
      symCEInterface_.GetPositionTolerance(),
      fractionalPositionTolerance_);

  // So each site will have one neighbourlist and will contains the lattice
  // sites upto the nearest neighbours such that the max distance b/w sites
  // is less than or equal to maxCutoff of clusters
}

vector<vector<vector<int>>> SymmetricCE::GetLocalOrbitsForLatticeSite(
    const size_t &latticeId)
{

  LocalOrbitListGenerator LOLG(
      clusterSpace_.getPrimitiveOrbitList(),
      supercellStructure_,
      symCEInterface_.GetFractionalPositionTolerance());

  Vector3d position = supercellStructure_->positionByIndex(latticeId);

  auto site = clusterSpace_.primitiveStructure().findLatticeSiteByPosition(
      position,
      fractionalPositionTolerance_);

  Vector3i offset = site.unitcellOffset();

  // Orbit List Object
  auto localOrbit = LOLG.getLocalOrbitList(
      offset, true);

  cout << localOrbit.size() << endl;

  // vector<vector<vector<int>>> Contains the orbits<equivalentClusters<latticeClusters>>
  auto localOrbitVector = clusterSpace_.getLocalOrbitClusterVector(
      localOrbit,
      supercellStructure_);

  return localOrbitVector;
}

vector<vector<vector<int>>> SymmetricCE::GetLocalOrbitsEncoding()
{
  auto centralLatticeId = supercellConfig_.GetCentralAtomLatticeId();

  // Since neighbour lists are updated for the maxClusterSize hence only one
  // neighbour list corresponding to each lattice site
  size_t maxBondOrder = 1;

  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
      supercellConfig_,
      centralLatticeId,
      maxBondOrder);

  // Must add the latticeId as the last element to be consistent
  // Also one can add it in the GetCanonicalSortedSitesForSite function
  // But for now lets keep these two separate
  canonicalSortedLatticeIds.emplace_back(centralLatticeId);

  cout << canonicalSortedLatticeIds.size() << endl;

  unordered_map<size_t, int> latticeIdToIndexMap;

  for (int i = 0; i < canonicalSortedLatticeIds.size(); i++)
  {
    latticeIdToIndexMap[canonicalSortedLatticeIds[i]] = i;
  }

  auto localOrbits = GetLocalOrbitsForLatticeSite(centralLatticeId);

  cout << localOrbits.size() << endl;

  vector<vector<vector<int>>> localOrbitsEncoding;
  localOrbitsEncoding.reserve(localOrbits.size());

  for (auto orbit : localOrbits)
  {
    vector<vector<int>> encodedOrbit;
    encodedOrbit.reserve(orbit.size());

    for (auto cluster : orbit)
    {
      vector<int> encodedCluster;
      encodedCluster.reserve(cluster.size());

      for (auto latticeId : cluster)
      {
        encodedCluster.emplace_back(latticeIdToIndexMap.at(latticeId));
      }
      encodedOrbit.emplace_back(encodedCluster);
    }
    localOrbitsEncoding.emplace_back(encodedOrbit);
  }

  return localOrbitsEncoding;
}

// Returns local cluster vector for a lattice Id for which canonicalSortedLatticeIdVector are given
vector<double> SymmetricCE::GetLocalClusterVector(
    const Config &config, // The current configuration
    const vector<size_t> &canonicalSortedLatticeIdVector,
    const vector<vector<vector<int>>> &localOrbitEncoding) const
{

  // Check that orbit lists match
  if (primitiveOrbitList_->size() != localOrbitEncoding.size())
  {
    ostringstream msg;
    msg << "Error in `SymmetricCE::GetLocalClusterVector`: Orbit lists do not match."
        << endl
        << localOrbitEncoding.size() << " >= " << primitiveOrbitList_->size() << endl;
    throw runtime_error(msg.str());
  }

  vector<double> clusterVector;

  // Then we will calculate a local cluster vector.
  // The zerolet (which is 1 in the full cluster vector)
  // can be considered as made up of equal contributions from
  // all sites.
  clusterVector.push_back(1.0);

  // Define some variables before loop for performance reasons
  map<vector<int>, double> clusterCounts;

  vector<int> allowedOccupations;
  vector<int> indicesOfRepresentativeSites;

  // Start to occupy the cluster vector orbit by orbit
  for (size_t orbitIndex = 0; orbitIndex < primitiveOrbitList_->size(); orbitIndex++)
  {

    const Orbit &primitiveOrbit = primitiveOrbitList_->getOrbit(orbitIndex);
    auto encodedSupercellOrbit = localOrbitEncoding[orbitIndex];

    // Count clusters
    // for a given orbit containing the equivalent cluster count the number of atom clusters
    // but seems like they are not counting the clusters canonically but the order or
    // permutation will be determined by the primitive orbit list

    // print2DVector(supercellOrbit);

    // This is basically atom cluster count for each orbit
    clusterCounts.clear();

    // Here instead of storing atomCluster as vector<int>
    // One can store it as AtomClusterType ?
    // But not sure about the order though where A-A-B-A is
    // taken to be equivalent as A-A-A-B or may be not sure about this

    for (auto encodedCluster : encodedSupercellOrbit)
    {
      vector<int> atomCluster;
      atomCluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        auto latticeId = canonicalSortedLatticeIdVector[encodedIdx];
        atomCluster.emplace_back(static_cast<int>(config.GetElementOfLattice(latticeId).GetAtomicIndex()));
      }
      double unit = 1;
      clusterCounts[atomCluster] += unit;
    }

    // Extract allowed occupations (needed to calculate point functions)
    allowedOccupations = primitiveOrbit.representativeCluster().getNumberOfAllowedSpeciesPerSite();

    // cout << "Allowed Occupations : ";
    // print1DVector(allowedOccupations);

    // Extract indices of representative sites
    // (needed to extract internal integer representation of each element)
    indicesOfRepresentativeSites.clear();
    for (const LatticeSite &site : primitiveOrbit.representativeCluster().latticeSites())
    {
      indicesOfRepresentativeSites.push_back(site.index());
    }

    //    cout << "Cluster Counts: " << endl;
    //    for (auto entry : clusterCounts)
    //    {
    //      print1DVector(entry.first);
    //      cout << " : " << entry.second << endl;
    //    }
    //
    //    cout << "Indices of representative sites : ";
    //    print1DVector(indicesOfRepresentativeSites);

    // Loop over all multi-component vectors for this orbit.
    // These are vectors of integers, where the integer represents
    // a cluster function index.
    //
    // Example 1: For a binary alloy we obtain [0, 0] and [0, 0, 0]
    // for pair and triplet terms, respectively.
    //
    // Example 2: For a ternary alloy we obtain [0, 0], [0, 1], [1, 0], and [1, 1]
    // for pairs (and similarly for triplets). However, if the two sites of the pair are
    // equivalent (which, for example, is always the case in systems with one atom
    // in the primitive cell, such as FCC) [0, 1] and [1, 0] are considered equivalent
    // and we will only get one of them.

    for (auto &cvElement : primitiveOrbit.clusterVectorElements())
    {
      //      cout << "MulitComponentVector : ";
      //      print1DVector(cvElement.multiComponentVector);
      //
      //      cout << endl;

      double clusterVectorElement = 0;

      /// Loop over all the counts for this orbit
      for (const auto &clusterCount : clusterCounts)
      {
        /// Loop over all permutations belonging to this multi-component vector
        for (const auto &permutation : cvElement.sitePermutations)
        {
          clusterVectorElement += clusterCount.second * clusterSpace_.evaluateClusterProduct(
                                                            cvElement.multiComponentVector,
                                                            allowedOccupations,
                                                            clusterCount.first,
                                                            indicesOfRepresentativeSites,
                                                            permutation);
        }
      }
      // Usually we could have counted multiplicity by simply adding the number of
      // clusters in the orbit (clusterCount.second), but in the case of
      // local cluster vectors or changes in cluster vectors, we have only counted
      // a subset of the clusters. We therefore use the pre-computed multiplicity.
      clusterVector.push_back(clusterVectorElement / (double)cvElement.multiplicity *
                              (double)primitiveOrbitList_->structure().size() /
                              (double)config.GetNumLattices());
    }
  }
  return clusterVector;
}

// Returns the local cluster vector for a site
//
vector<double> SymmetricCE::GetLocalClusterVector(
    const Config &config, // The current configuration
    const size_t &targetLatticeId,
    const Element &elementToAssign, // Element which need to be assigned to centralLatticeId
    const vector<size_t> &canonicalSortedLatticeIdVector,
    const vector<vector<vector<int>>> &localOrbitEncoding) const
{

  // Check that orbit lists match
  if (primitiveOrbitList_->size() != localOrbitEncoding.size())
  {
    ostringstream msg;
    msg << "Error in `SymmetricCE::GetLocalClusterVector`: Orbit lists do not match."
        << endl
        << localOrbitEncoding.size() << " >= " << primitiveOrbitList_->size() << endl;
    throw runtime_error(msg.str());
  }

  vector<double> clusterVector;

  // Then we will calculate a local cluster vector.
  // The zerolet (which is 1 in the full cluster vector)
  // can be considered as made up of equal contributions from
  // all sites.
  clusterVector.push_back(1.0);

  // Define some variables before loop for performance reasons
  map<vector<int>, double> clusterCounts;

  vector<int> allowedOccupations;
  vector<int> indicesOfRepresentativeSites;

  // Start to occupy the cluster vector orbit by orbit
  for (size_t orbitIndex = 0; orbitIndex < primitiveOrbitList_->size(); orbitIndex++)
  {

    const Orbit &primitiveOrbit = primitiveOrbitList_->getOrbit(orbitIndex);
    auto encodedSupercellOrbit = localOrbitEncoding[orbitIndex];

    // Count clusters
    // for a given orbit containing the equivalent cluster count the number of atom clusters
    // but seems like they are not counting the clusters canonically but the order or
    // permutation will be determined by the primitive orbit list

    // print2DVector(supercellOrbit);

    // This is basically atom cluster count for each orbit
    clusterCounts.clear();

    // Here instead of storing atomCluster as vector<int>
    // One can store it as AtomClusterType ?
    // But not sure about the order though where A-A-B-A is
    // taken to be equivalent as A-A-A-B or may be not sure about this

    for (auto encodedCluster : encodedSupercellOrbit)
    {
      vector<int> atomCluster;
      atomCluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        auto latticeId = canonicalSortedLatticeIdVector[encodedIdx];
        if (latticeId == targetLatticeId)
        {
          atomCluster.emplace_back(static_cast<int>(elementToAssign.GetAtomicIndex()));
        }
        else
        {
          atomCluster.emplace_back(static_cast<int>(config.GetElementOfLattice(latticeId).GetAtomicIndex()));
        }
      }
      double unit = 1;
      clusterCounts[atomCluster] += unit;
    }

    // Extract allowed occupations (needed to calculate point functions)
    allowedOccupations = primitiveOrbit.representativeCluster().getNumberOfAllowedSpeciesPerSite();

    // cout << "Allowed Occupations : ";
    // print1DVector(allowedOccupations);

    // Extract indices of representative sites
    // (needed to extract internal integer representation of each element)
    indicesOfRepresentativeSites.clear();
    for (const LatticeSite &site : primitiveOrbit.representativeCluster().latticeSites())
    {
      indicesOfRepresentativeSites.push_back(site.index());
    }

    //    cout << "Cluster Counts: " << endl;
    //    for (auto entry : clusterCounts)
    //    {
    //      print1DVector(entry.first);
    //      cout << " : " << entry.second << endl;
    //    }
    //
    //    cout << "Indices of representative sites : ";
    //    print1DVector(indicesOfRepresentativeSites);

    // Loop over all multi-component vectors for this orbit.
    // These are vectors of integers, where the integer represents
    // a cluster function index.
    //
    // Example 1: For a binary alloy we obtain [0, 0] and [0, 0, 0]
    // for pair and triplet terms, respectively.
    //
    // Example 2: For a ternary alloy we obtain [0, 0], [0, 1], [1, 0], and [1, 1]
    // for pairs (and similarly for triplets). However, if the two sites of the pair are
    // equivalent (which, for example, is always the case in systems with one atom
    // in the primitive cell, such as FCC) [0, 1] and [1, 0] are considered equivalent
    // and we will only get one of them.

    for (auto &cvElement : primitiveOrbit.clusterVectorElements())
    {
      //      cout << "MulitComponentVector : ";
      //      print1DVector(cvElement.multiComponentVector);
      //
      //      cout << endl;

      double clusterVectorElement = 0;

      /// Loop over all the counts for this orbit
      for (const auto &clusterCount : clusterCounts)
      {
        /// Loop over all permutations belonging to this multi-component vector
        for (const auto &permutation : cvElement.sitePermutations)
        {
          clusterVectorElement += clusterCount.second * clusterSpace_.evaluateClusterProduct(
                                                            cvElement.multiComponentVector,
                                                            allowedOccupations,
                                                            clusterCount.first,
                                                            indicesOfRepresentativeSites,
                                                            permutation);
        }
      }
      // Usually we could have counted multiplicity by simply adding the number of
      // clusters in the orbit (clusterCount.second), but in the case of
      // local cluster vectors or changes in cluster vectors, we have only counted
      // a subset of the clusters. We therefore use the pre-computed multiplicity.
      clusterVector.push_back(clusterVectorElement / (double)cvElement.multiplicity *
                              (double)primitiveOrbitList_->structure().size() /
                              (double)config.GetNumLattices());
    }
  }
  return clusterVector;
}

//
vector<double> SymmetricCE::GetClusterVectorForConfig(
    const Config &config)
{
  auto structure = ConvertConfigToStructure(config);

  auto clusterVector = clusterSpace_.getClusterVector(
      structure, fractionalPositionTolerance_);

  return clusterVector;
}

/*
Returns the local cluster vector centered around a given lattice site, taking
into account different possible elements at the specified centralLatticeId.

This is particularly useful in Kinetic Monte Carlo (KMC) simulations for computing the
change in energy (dE) using the local cluster expansion, for example when modeling
vacancy diffusion in dilute systems.

For instance, if we consider two elements which can occupy (A and B) site1,
the function will return the cluster vector corresponding to each element in the
specified order. These vectors can be used to compute terms like Ev in the
energy difference:

    E(LVFE) = Ev - 1/2 (EA + EB)

Also since one will be taking the average so order will not matter and also this
is just local contribution as we are interested in the dE not the absolute energy
Not sure about this but that is the first feeling
*/

vector<vector<double>> SymmetricCE::GetMultiElementLocalClusterVector(
    const Config &config, // The current configuration
    const size_t &centralLatticeId,
    const vector<Element> &elementVector,
    const vector<size_t> &canonicalSortedLatticeIdVector,
    const vector<vector<vector<int>>> &localOrbitEncoding) const
{
  // Check that orbit lists match
  if (primitiveOrbitList_->size() != localOrbitEncoding.size())
  {
    ostringstream msg;
    msg << "Error in `SymmetricCE::GetLocalClusterVector`: Orbit lists do not match."
        << endl
        << localOrbitEncoding.size() << " >= " << primitiveOrbitList_->size() << endl;
    throw runtime_error(msg.str());
  }

  const size_t numElements = elementVector.size();
  if (numElements == 0)
  {
    // Handle edge case if needed, but assuming at least one element
    throw runtime_error("Error: elementVector cannot be empty.");
  }

  vector<vector<double>> multiElementClusterVectors(numElements);
  for (auto &cv : multiElementClusterVectors)
  {
    // The zerolet (which is 1 in the full cluster vector)
    // can be considered as made up of equal contributions from
    // all sites.
    cv.push_back(1.0);
  }

  // Process each orbit
  for (size_t orbitIndex = 0; orbitIndex < primitiveOrbitList_->size(); ++orbitIndex)
  {
    const Orbit &primitiveOrbits = primitiveOrbitList_->getOrbit(orbitIndex);
    const auto &encodedSupercellOrbits = localOrbitEncoding[orbitIndex];

    // Count clusters for each element variant
    vector<map<vector<int>, double>> clusterCounts(numElements);

    for (const auto &encodedCluster : encodedSupercellOrbits)
    {
      const size_t clusterSize = encodedCluster.size();
      vector<vector<int>> atomClusters(numElements);
      for (auto &ac : atomClusters)
      {
        ac.reserve(clusterSize);
      }

      for (const auto encodedIndex : encodedCluster)
      {
        const size_t latticeId = canonicalSortedLatticeIdVector[encodedIndex];

        if (latticeId == centralLatticeId)
        {
          for (size_t elementIndex = 0; elementIndex < numElements; ++elementIndex)
          {
            const int atomicNumber = static_cast<int>(elementVector[elementIndex].GetAtomicIndex());
            atomClusters[elementIndex].emplace_back(atomicNumber);
          }
        }
        else
        {
          // Use actual configuration for non-central sites
          const int atomicNumber = static_cast<int>(config.GetElementOfLattice(latticeId).GetAtomicIndex());
          for (size_t elementIndex = 0; elementIndex < numElements; ++elementIndex)
          {
            atomClusters[elementIndex].emplace_back(atomicNumber);
          }
        }
      }

      const double unit = 1.0;
      for (size_t elementIndex = 0; elementIndex < numElements; ++elementIndex)
      {
        clusterCounts[elementIndex][atomClusters[elementIndex]] += unit;
      }
    }

    // Extract allowed occupations (needed to calculate point functions)
    const vector<int> allowedOccupations = primitiveOrbits.representativeCluster().getNumberOfAllowedSpeciesPerSite();

    // Extract indices of representative sites
    // (needed to extract internal integer representation of each element)
    vector<int> indicesOfRepresentativeSites;
    indicesOfRepresentativeSites.reserve(primitiveOrbits.representativeCluster().latticeSites().size());
    for (const LatticeSite &site : primitiveOrbits.representativeCluster().latticeSites())
    {
      indicesOfRepresentativeSites.push_back(site.index());
    }

    // Loop over all multi-component vectors for this orbit.
    // These are vectors of integers, where the integer represents
    // a cluster function index.
    //
    // Example 1: For a binary alloy we obtain [0, 0] and [0, 0, 0]
    // for pair and triplet terms, respectively.
    //
    // Example 2: For a ternary alloy we obtain [0, 0], [0, 1], [1, 0], and [1, 1]
    // for pairs (and similarly for triplets). However, if the two sites of the pair are
    // equivalent (which, for example, is always the case in systems with one atom
    // in the primitive cell, such as FCC) [0, 1] and [1, 0] are considered equivalent
    // and we will only get one of them.
    for (const auto &cvElement : primitiveOrbits.clusterVectorElements())
    {
      vector<double> clusterVectorElements(numElements, 0.0);

      for (size_t elementIndex = 0; elementIndex < numElements; ++elementIndex)
      {
        const auto &counts = clusterCounts[elementIndex];

        // Loop over all the counts for this orbit
        for (const auto &cluster_count : counts)
        {
          // Loop over all permutations belonging to this multi-component vector
          for (const auto &permutation : cvElement.sitePermutations)
          {
            clusterVectorElements[elementIndex] += cluster_count.second * clusterSpace_.evaluateClusterProduct(
                                                                              cvElement.multiComponentVector,
                                                                              allowedOccupations,
                                                                              cluster_count.first,
                                                                              indicesOfRepresentativeSites,
                                                                              permutation);
          }
        }
      }

      // Normalization: Usually we could have counted multiplicity by simply adding the number of
      // clusters in the orbit (clusterCount.second), but in the case of
      // local cluster vectors or changes in cluster vectors, we have only counted
      // a subset of the clusters. We therefore use the pre-computed multiplicity.
      const double normalizationFactor = (1.0 / static_cast<double>(cvElement.multiplicity)) *
                                         (static_cast<double>(primitiveOrbitList_->structure().size()) /
                                          static_cast<double>(config.GetNumLattices()));

      for (size_t elementIndex = 0; elementIndex < numElements; ++elementIndex)
      {
        multiElementClusterVectors[elementIndex].push_back(clusterVectorElements[elementIndex] * normalizationFactor);
      }
    }
  }

  return multiElementClusterVectors;
}

void SymmetricCE::TestSymmetricCE(
    Config &config,
    const vector<double> &eciVector,
    const pair<size_t, size_t> &latticeIdPair)
{
  auto localOrbitEncoding = GetLocalOrbitsEncoding();

  auto canonicalSortedSitesLatticeId1 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdPair.first,
      1);
  canonicalSortedSitesLatticeId1.emplace_back(latticeIdPair.first);

  auto canonicalSortedSitesLatticeId2 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdPair.second,
      1);
  canonicalSortedSitesLatticeId2.emplace_back(latticeIdPair.second);

  auto nAtoms = config.GetNumAtoms();

  // Lambda to compute formation energy
  auto computeFormationEnergy = [](const vector<double> &icetEci,
                                   const vector<double> &ceClusterVector,
                                   const int nAtoms)
  {
    // Will return total formation energy
    return nAtoms * inner_product(icetEci.begin(), icetEci.end(),
                                  ceClusterVector.begin(), 0.0);
  };

  // Chemical potentials (example: Mo, Ta)
  unordered_map<string, double> chemicalPotentials{
      {"Mo", -10.93308598}, // Mo (Z=42)
      {"Ta", -11.81233671}  // Ta (Z=73)
  };

  // Lambda to compute total energy = formation + sum(N_i * mu_i)
  auto computeTotalEnergy = [&](double formationEnergy,
                                const unordered_map<string, int> &elementCountMap)
  {
    double muContribution = 0.0;
    for (const auto &[element, count] : elementCountMap)
    {
      auto it = chemicalPotentials.find(element);
      if (it != chemicalPotentials.end())
      {
        muContribution += count * it->second;
      }
    }
    return formationEnergy + muContribution;
  };

  // -------------------------------
  // Before Jump
  // -------------------------------
  auto structureBeforeJump = ConvertConfigToStructure(config);

  auto cvSupercellBeforeJump = clusterSpace_.getClusterVector(
      structureBeforeJump, fractionalPositionTolerance_);

  auto clusterVector1 = GetLocalClusterVector(config,
                                              canonicalSortedSitesLatticeId1,
                                              localOrbitEncoding);

  auto clusterVector2 = GetLocalClusterVector(config,
                                              canonicalSortedSitesLatticeId2,
                                              localOrbitEncoding);

  double formationEnergyBefore = computeFormationEnergy(eciVector,
                                                        cvSupercellBeforeJump,
                                                        nAtoms);

  double formationEnergyLocalBefore = computeFormationEnergy(eciVector,
                                                             clusterVector1,
                                                             nAtoms) +
                                      computeFormationEnergy(eciVector,
                                                             clusterVector2,
                                                             nAtoms);

  unordered_map<string, int> elementCountMapBefore;
  for (auto ele : config.GetAtomVector())
  {
    elementCountMapBefore[ele.GetElementString()]++;
  }

  double totalEnergyBefore = computeTotalEnergy(formationEnergyBefore,
                                                elementCountMapBefore);

  cout << "=== Before Jump ===\n";
  cout << "Supercell Formation Energy: " << formationEnergyBefore << "\n";
  cout << "Local Formation Energy (site1 + site2): " << formationEnergyLocalBefore << "\n";
  cout << "Supercell Total Energy: " << totalEnergyBefore << "\n";

  // -------------------------------
  // After Jump
  // -------------------------------
  config.LatticeJump(latticeIdPair);

  auto structureAfterJump = ConvertConfigToStructure(config);

  auto cvSupercellAfterJump = clusterSpace_.getClusterVector(
      structureAfterJump, fractionalPositionTolerance_);

  auto clusterVector1After = GetLocalClusterVector(config,
                                                   canonicalSortedSitesLatticeId1,
                                                   localOrbitEncoding);

  auto clusterVector2After = GetLocalClusterVector(config,
                                                   canonicalSortedSitesLatticeId2,
                                                   localOrbitEncoding);

  double formationEnergyAfter = computeFormationEnergy(eciVector,
                                                       cvSupercellAfterJump,
                                                       nAtoms);

  double formationEnergyLocalAfter = computeFormationEnergy(eciVector,
                                                            clusterVector1After,
                                                            nAtoms) +
                                     computeFormationEnergy(eciVector,
                                                            clusterVector2After,
                                                            nAtoms);

  unordered_map<string, int> elementCountMapAfter;
  for (auto ele : config.GetAtomVector())
  {
    elementCountMapAfter[ele.GetElementString()]++;
  }

  double totalEnergyAfter = computeTotalEnergy(formationEnergyAfter,
                                               elementCountMapAfter);

  cout << "=== After Jump ===\n";
  cout << "Supercell Formation Energy: " << formationEnergyAfter << "\n";
  cout << "Local Formation Energy (site1 + site2): " << formationEnergyLocalAfter << "\n";
  cout << "Supercell Total Energy: " << totalEnergyAfter << "\n";

  // -------------------------------
  // Differences
  // -------------------------------
  cout << "=== Energy Difference (After - Before) ===\n";
  cout << "Δ Formation Energy (Supercell): "
       << (formationEnergyAfter - formationEnergyBefore) << "\n";
  cout << "Δ Formation Energy (Local): "
       << (formationEnergyLocalAfter - formationEnergyLocalBefore) << "\n";
  cout << "Δ Total Energy (Supercell): "
       << (totalEnergyAfter - totalEnergyBefore) << "\n";

  // -------------------------------
  // Restore original configuration
  // -------------------------------
  config.LatticeJump(latticeIdPair); // undo the jump

  // Recompute energies to verify restoration
  auto structureRestored = ConvertConfigToStructure(config);
  auto cvSupercellRestored = clusterSpace_.getClusterVector(
      structureRestored, fractionalPositionTolerance_);

  auto clusterVector1Restored = GetLocalClusterVector(config,
                                                      canonicalSortedSitesLatticeId1,
                                                      localOrbitEncoding);

  auto clusterVector2Restored = GetLocalClusterVector(config,
                                                      canonicalSortedSitesLatticeId2,
                                                      localOrbitEncoding);

  double formationEnergyRestored = computeFormationEnergy(eciVector,
                                                          cvSupercellRestored,
                                                          nAtoms);

  double formationEnergyLocalRestored = computeFormationEnergy(eciVector,
                                                               clusterVector1Restored,
                                                               nAtoms) +
                                        computeFormationEnergy(eciVector,
                                                               clusterVector2Restored,
                                                               nAtoms);

  unordered_map<string, int> elementCountMapRestored;
  for (auto ele : config.GetAtomVector())
  {
    elementCountMapRestored[ele.GetElementString()]++;
  }

  double totalEnergyRestored = computeTotalEnergy(formationEnergyRestored,
                                                  elementCountMapRestored);

  cout << "=== After Restoring Jump ===\n";
  cout << "Supercell Formation Energy: " << formationEnergyRestored << "\n";
  cout << "Local Formation Energy (site1 + site2): " << formationEnergyLocalRestored << "\n";
  cout << "Supercell Total Energy: " << totalEnergyRestored << "\n";

  cout << "Δ Formation Energy vs original: "
       << formationEnergyRestored - formationEnergyBefore << "\n";
  cout << "Δ Total Energy vs original: "
       << totalEnergyRestored - totalEnergyBefore << "\n";
}

// This function returns the cluster vector sum for a latticeIdPair after setting
// element to the desired elements

// Returns local cluster vector for a lattice Id for which canonicalSortedLatticeIdVector are given
pair<vector<double>, vector<double>> SymmetricCE::GetLocalClusterVectorForPairSerial(
    const Config &config, // The current configuration
    const pair<size_t, size_t> &latticeIdPair,
    const pair<Element, Element> &latticeIdPairElements,   // Elements which will be used for latticeIdPair (Order is important)
    const vector<size_t> &canonicalSortedLatticeIdVector1, // sorted lattice Ids for Id 1
    const vector<size_t> &canonicalSortedLatticeIdVector2, // sorted lattice Ids for Id 2
    const vector<vector<vector<int>>> &localOrbitEncoding) const
{
  // Check that orbit lists match
  if (primitiveOrbitList_->size() != localOrbitEncoding.size())
  {
    ostringstream msg;
    msg << "Error in `SymmetricCE::GetLocalClusterVector`: Orbit lists do not match."
        << endl
        << localOrbitEncoding.size() << " >= " << primitiveOrbitList_->size() << endl;
    throw runtime_error(msg.str());
  }

  // 1 AT THE END OF VARIABLE CORRESPONDS TO latticeIdPair.first
  // 2 AT THE END OF VARIABLE CORRESPONDS TO latticeIdPair.second

  vector<double> clusterVector1; // corresponds to latticeIdPair.first
  vector<double> clusterVector2; // corresponds to latticeIdPair.second

  // Then we will calculate a local cluster vector.
  // The zerolet (which is 1 in the full cluster vector)
  // can be considered as made up of equal contributions from
  // all sites.

  // Constant term
  clusterVector1.push_back(1.0);
  clusterVector2.push_back(1.0);

  // Define some variables before loop for performance reasons
  map<vector<int>, double> clusterCounts1;
  map<vector<int>, double> clusterCounts2;

  vector<int> allowedOccupations;
  vector<int> indicesOfRepresentativeSites;

  // Start to occupy the cluster vector orbit by orbit
  for (size_t orbitIndex = 0; orbitIndex < primitiveOrbitList_->size(); orbitIndex++)
  {

    // Will be same for both the latticeIds
    const Orbit &primitiveOrbit = primitiveOrbitList_->getOrbit(orbitIndex);
    auto encodedSupercellOrbit = localOrbitEncoding[orbitIndex];

    // Count clusters
    // for a given orbit containing the equivalent cluster count the number of atom clusters
    // but seems like they are not counting the clusters canonically but the order or
    // permutation will be determined by the primitive orbit list

    // print2DVector(supercellOrbit);

    // This is basically atom cluster count for each orbit
    clusterCounts1.clear();
    clusterCounts2.clear();

    // Here instead of storing atomCluster as vector<int>
    // One can store it as AtomClusterType ?
    // But not sure about the order though where A-A-B-A is
    // taken to be equivalent as A-A-A-B or may be not sure about this

    for (auto encodedCluster : encodedSupercellOrbit)
    {
      vector<int> atomCluster1;
      vector<int> atomCluster2;

      atomCluster1.reserve(encodedCluster.size());
      atomCluster2.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        // Atom cluster for site 1
        auto latticeId1 = canonicalSortedLatticeIdVector1[encodedIdx];

        // If latticeId is same as first site Id assign first element
        if (latticeId1 == latticeIdPair.first)
        {
          atomCluster1.emplace_back(static_cast<int>(latticeIdPairElements.first.GetAtomicIndex()));
        }
        // If latticeId is same as second site Id assign second element
        else if (latticeId1 == latticeIdPair.second)
        {
          atomCluster1.emplace_back(static_cast<int>(latticeIdPairElements.second.GetAtomicIndex()));
        }
        else
        {
          // Use the element from config object
          atomCluster1.emplace_back(static_cast<int>(config.GetElementOfLattice(latticeId1).GetAtomicIndex()));
        }

        // Atom cluster for site 2
        auto latticeId2 = canonicalSortedLatticeIdVector2[encodedIdx];

        // If latticeId is same as first site Id assign first element
        if (latticeId2 == latticeIdPair.first)
        {
          atomCluster2.emplace_back(static_cast<int>(latticeIdPairElements.first.GetAtomicIndex()));
        }
        // If latticeId is same as second site Id assign second element
        else if (latticeId2 == latticeIdPair.second)
        {
          atomCluster2.emplace_back(static_cast<int>(latticeIdPairElements.second.GetAtomicIndex()));
        }
        else
        {
          // Use the element from config object
          atomCluster2.emplace_back(static_cast<int>(config.GetElementOfLattice(latticeId2).GetAtomicIndex()));
        }
      }

      double unit = 1;
      clusterCounts1[atomCluster1] += unit;
      clusterCounts2[atomCluster2] += unit;
    }

    // Extract allowed occupations (needed to calculate point functions)
    allowedOccupations = primitiveOrbit.representativeCluster().getNumberOfAllowedSpeciesPerSite();

    // cout << "Allowed Occupations : ";
    // print1DVector(allowedOccupations);

    // Extract indices of representative sites
    // (needed to extract internal integer representation of each element)
    indicesOfRepresentativeSites.clear();
    for (const LatticeSite &site : primitiveOrbit.representativeCluster().latticeSites())
    {
      indicesOfRepresentativeSites.push_back(site.index());
    }

    // Loop over all multi-component vectors for this orbit.
    // These are vectors of integers, where the integer represents
    // a cluster function index.
    //
    // Example 1: For a binary alloy we obtain [0, 0] and [0, 0, 0]
    // for pair and triplet terms, respectively.
    //
    // Example 2: For a ternary alloy we obtain [0, 0], [0, 1], [1, 0], and [1, 1]
    // for pairs (and similarly for triplets). However, if the two sites of the pair are
    // equivalent (which, for example, is always the case in systems with one atom
    // in the primitive cell, such as FCC) [0, 1] and [1, 0] are considered equivalent
    // and we will only get one of them.

    for (auto &cvElement : primitiveOrbit.clusterVectorElements())
    {
      // For Site 1
      double clusterVectorElement1 = 0;

      /// Loop over all the counts for this orbit
      for (const auto &clusterCount : clusterCounts1)
      {
        /// Loop over all permutations belonging to this multi-component vector
        for (const auto &permutation : cvElement.sitePermutations)
        {
          clusterVectorElement1 += clusterCount.second * clusterSpace_.evaluateClusterProduct(
                                                             cvElement.multiComponentVector,
                                                             allowedOccupations,
                                                             clusterCount.first,
                                                             indicesOfRepresentativeSites,
                                                             permutation);
        }
      }
      // Usually we could have counted multiplicity by simply adding the number of
      // clusters in the orbit (clusterCount.second), but in the case of
      // local cluster vectors or changes in cluster vectors, we have only counted
      // a subset of the clusters. We therefore use the pre-computed multiplicity.
      clusterVector1.push_back(clusterVectorElement1 / (double)cvElement.multiplicity *
                               (double)primitiveOrbitList_->structure().size() /
                               (double)config.GetNumLattices());

      // For Site2

      double clusterVectorElement2 = 0;

      /// Loop over all the counts for this orbit
      for (const auto &clusterCount : clusterCounts2)
      {
        /// Loop over all permutations belonging to this multi-component vector
        for (const auto &permutation : cvElement.sitePermutations)
        {
          clusterVectorElement2 += clusterCount.second * clusterSpace_.evaluateClusterProduct(
                                                             cvElement.multiComponentVector,
                                                             allowedOccupations,
                                                             clusterCount.first,
                                                             indicesOfRepresentativeSites,
                                                             permutation);
        }
      }
      // Usually we could have counted multiplicity by simply adding the number of
      // clusters in the orbit (clusterCount.second), but in the case of
      // local cluster vectors or changes in cluster vectors, we have only counted
      // a subset of the clusters. We therefore use the pre-computed multiplicity.
      clusterVector2.push_back(clusterVectorElement2 / (double)cvElement.multiplicity *
                               (double)primitiveOrbitList_->structure().size() /
                               (double)config.GetNumLattices());
    }
  }

  return {clusterVector1, clusterVector2};
}

// Returns local cluster vector for a lattice Id for which canonicalSortedLatticeIdVector are given
// Parallel Version
pair<vector<double>, vector<double>> SymmetricCE::GetLocalClusterVectorForPair(
    const Config &config, // The current configuration
    const pair<size_t, size_t> &latticeIdPair,
    const pair<Element, Element> &latticeIdPairElements,   // Elements which will be used for latticeIdPair (Order is important)
    const vector<size_t> &canonicalSortedLatticeIdVector1, // sorted lattice Ids for Id 1
    const vector<size_t> &canonicalSortedLatticeIdVector2, // sorted lattice Ids for Id 2
    const vector<vector<vector<int>>> &localOrbitEncoding) const
{
  // Check that orbit lists match
  if (primitiveOrbitList_->size() != localOrbitEncoding.size())
  {
    ostringstream msg;
    msg << "Error in `SymmetricCE::GetLocalClusterVector`: Orbit lists do not match."
        << endl
        << localOrbitEncoding.size() << " >= " << primitiveOrbitList_->size() << endl;
    throw runtime_error(msg.str());
  }

  // 1 AT THE END OF VARIABLE CORRESPONDS TO latticeIdPair.first
  // 2 AT THE END OF VARIABLE CORRESPONDS TO latticeIdPair.second

  vector<double> clusterVector1; // corresponds to latticeIdPair.first
  vector<double> clusterVector2; // corresponds to latticeIdPair.second

  // Then we will calculate a local cluster vector.
  // The zerolet (which is 1 in the full cluster vector)
  // can be considered as made up of equal contributions from
  // all sites.

  
  // Compute total size for preallocation
  size_t total_elements = 1; // constant term
  for (size_t orbitIndex = 0; orbitIndex < primitiveOrbitList_->size(); ++orbitIndex)
  {
    const Orbit &primitiveOrbit = primitiveOrbitList_->getOrbit(orbitIndex);
    total_elements += primitiveOrbit.clusterVectorElements().size();
  }
  

  clusterVector1.resize(total_elements);
  clusterVector2.resize(total_elements);

  // Constant term
  clusterVector1[0] = 1.0;
  clusterVector2[0] = 1.0;

  // Compute offsets for each orbit's starting position in the vectors
  vector<size_t> offsets(primitiveOrbitList_->size());
  size_t current_pos = 1;
  for (size_t orbitIndex = 0; orbitIndex < primitiveOrbitList_->size(); ++orbitIndex)
  {
    offsets[orbitIndex] = current_pos;
    const Orbit &primitiveOrbit = primitiveOrbitList_->getOrbit(orbitIndex);
    current_pos += primitiveOrbit.clusterVectorElements().size();
  }

  // Parallel loop over orbits
#pragma omp parallel for schedule(dynamic)
  for (size_t orbitIndex = 0; orbitIndex < primitiveOrbitList_->size(); ++orbitIndex)
  {
    // Will be same for both the latticeIds
    const Orbit &primitiveOrbit = primitiveOrbitList_->getOrbit(orbitIndex);
    auto encodedSupercellOrbit = localOrbitEncoding[orbitIndex];

    // Count clusters
    // for a given orbit containing the equivalent cluster count the number of atom clusters
    // but seems like they are not counting the clusters canonically but the order or
    // permutation will be determined by the primitive orbit list

    // This is basically atom cluster count for each orbit
    map<vector<int>, double> clusterCounts1;
    map<vector<int>, double> clusterCounts2;

    // Here instead of storing atomCluster as vector<int>
    // One can store it as AtomClusterType ?
    // But not sure about the order though where A-A-B-A is
    // taken to be equivalent as A-A-A-B or may be not sure about this

    for (auto encodedCluster : encodedSupercellOrbit)
    {
      vector<int> atomCluster1;
      vector<int> atomCluster2;

      atomCluster1.reserve(encodedCluster.size());
      atomCluster2.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        // Atom cluster for site 1
        auto latticeId1 = canonicalSortedLatticeIdVector1[encodedIdx];

        // If latticeId is same as first site Id assign first element
        if (latticeId1 == latticeIdPair.first)
        {
          atomCluster1.emplace_back(static_cast<int>(latticeIdPairElements.first.GetAtomicIndex()));
        }
        // If latticeId is same as second site Id assign second element
        else if (latticeId1 == latticeIdPair.second)
        {
          atomCluster1.emplace_back(static_cast<int>(latticeIdPairElements.second.GetAtomicIndex()));
        }
        else
        {
          // Use the element from config object
          atomCluster1.emplace_back(static_cast<int>(config.GetElementOfLattice(latticeId1).GetAtomicIndex()));
        }

        // Atom cluster for site 2
        auto latticeId2 = canonicalSortedLatticeIdVector2[encodedIdx];

        // If latticeId is same as first site Id assign first element
        if (latticeId2 == latticeIdPair.first)
        {
          atomCluster2.emplace_back(static_cast<int>(latticeIdPairElements.first.GetAtomicIndex()));
        }
        // If latticeId is same as second site Id assign second element
        else if (latticeId2 == latticeIdPair.second)
        {
          atomCluster2.emplace_back(static_cast<int>(latticeIdPairElements.second.GetAtomicIndex()));
        }
        else
        {
          // Use the element from config object
          atomCluster2.emplace_back(static_cast<int>(config.GetElementOfLattice(latticeId2).GetAtomicIndex()));
        }
      }

      double unit = 1;
      clusterCounts1[atomCluster1] += unit;
      clusterCounts2[atomCluster2] += unit;
    }

    // Extract allowed occupations (needed to calculate point functions)
    vector<int> allowedOccupations = primitiveOrbit.representativeCluster().getNumberOfAllowedSpeciesPerSite();

    // cout << "Allowed Occupations : ";
    // print1DVector(allowedOccupations);

    // Extract indices of representative sites
    // (needed to extract internal integer representation of each element)
    vector<int> indicesOfRepresentativeSites;
    for (const LatticeSite &site : primitiveOrbit.representativeCluster().latticeSites())
    {
      indicesOfRepresentativeSites.push_back(site.index());
    }

    // Loop over all multi-component vectors for this orbit.
    // These are vectors of integers, where the integer represents
    // a cluster function index.
    //
    // Example 1: For a binary alloy we obtain [0, 0] and [0, 0, 0]
    // for pair and triplet terms, respectively.
    //
    // Example 2: For a ternary alloy we obtain [0, 0], [0, 1], [1, 0], and [1, 1]
    // for pairs (and similarly for triplets). However, if the two sites of the pair are
    // equivalent (which, for example, is always the case in systems with one atom
    // in the primitive cell, such as FCC) [0, 1] and [1, 0] are considered equivalent
    // and we will only get one of them.

    size_t pos = offsets[orbitIndex];
    for (auto &cvElement : primitiveOrbit.clusterVectorElements())
    {
      // For Site 1
      double clusterVectorElement1 = 0;

      /// Loop over all the counts for this orbit
      for (const auto &clusterCount : clusterCounts1)
      {
        /// Loop over all permutations belonging to this multi-component vector
        for (const auto &permutation : cvElement.sitePermutations)
        {
          clusterVectorElement1 += clusterCount.second * clusterSpace_.evaluateClusterProduct(
                                                             cvElement.multiComponentVector,
                                                             allowedOccupations,
                                                             clusterCount.first,
                                                             indicesOfRepresentativeSites,
                                                             permutation);
        }
      }
      // Usually we could have counted multiplicity by simply adding the number of
      // clusters in the orbit (clusterCount.second), but in the case of
      // local cluster vectors or changes in cluster vectors, we have only counted
      // a subset of the clusters. We therefore use the pre-computed multiplicity.
      clusterVector1[pos] = clusterVectorElement1 / (double)cvElement.multiplicity *
                            (double)primitiveOrbitList_->structure().size() /
                            (double)config.GetNumLattices();

      // For Site2

      double clusterVectorElement2 = 0;

      /// Loop over all the counts for this orbit
      for (const auto &clusterCount : clusterCounts2)
      {
        /// Loop over all permutations belonging to this multi-component vector
        for (const auto &permutation : cvElement.sitePermutations)
        {
          clusterVectorElement2 += clusterCount.second * clusterSpace_.evaluateClusterProduct(
                                                             cvElement.multiComponentVector,
                                                             allowedOccupations,
                                                             clusterCount.first,
                                                             indicesOfRepresentativeSites,
                                                             permutation);
        }
      }
      // Usually we could have counted multiplicity by simply adding the number of
      // clusters in the orbit (clusterCount.second), but in the case of
      // local cluster vectors or changes in cluster vectors, we have only counted
      // a subset of the clusters. We therefore use the pre-computed multiplicity.
      clusterVector2[pos] = clusterVectorElement2 / (double)cvElement.multiplicity *
                            (double)primitiveOrbitList_->structure().size() /
                            (double)config.GetNumLattices();

      pos++;
    }
  }

  return {clusterVector1, clusterVector2};
}