#include "CorrelationVector.h"

// Make use of symmetry
/*
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  vector<VectorXd> correlationFunctions;
  correlationFunctions.reserve(equivalentEncodedClusters.size());

  for (const auto &encodedOrbitPair : equivalentEncodedClusters)
  {
    
    if (encodedOrbitPair.first.size() == 1)
    {
      if (encodedOrbitPair.first[0].empty())
      {
        VectorXd orbitCorrelationFunction(1);
        orbitCorrelationFunction(0) = 1.0; // φ₀ = 1

        correlationFunctions.emplace_back(orbitCorrelationFunction);

        continue;
      }
    }

    VectorXd orbitCorrelationFunction = GetCorrelationFunction(
        config,
        atomicBasis,
        canonicalSortedLatticeIds,
        encodedOrbitPair.first);

    
    if (orbitCorrelationFunction.size() > 0)
      correlationFunctions.emplace_back(orbitCorrelationFunction);
  }

  // Concat the correlationFunctions and return a VectorXd

  size_t totalSize = 0;
  for (const auto &v : correlationFunctions)
    totalSize += v.size();

  VectorXd correlationVector(totalSize);

  size_t offset = 0;
  for (const auto &v : correlationFunctions)
  {
    correlationVector.segment(offset, v.size()) = v;
    offset += v.size();
  }

  return correlationVector;
<<<<<<< Updated upstream
=======
}
*/
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

// GetCorrelationVector with a element at specific site
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        targetLatticeId,
        elementToAssign,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}
*/
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

// GetCorrelationVector with a element at specific site
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        targetLatticeId,
        elementToAssign,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}
*/
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

// GetCorrelationVector with a element at specific site
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        targetLatticeId,
        elementToAssign,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}
*/
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

// GetCorrelationVector with a element at specific site
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        targetLatticeId,
        elementToAssign,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}
*/
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

// GetCorrelationVector with a element at specific site
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        targetLatticeId,
        elementToAssign,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentOrbitVector)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentOrbitVector.size());

  for (size_t idx = 0; idx < equivalentOrbitVector.size(); ++idx)
  {
    const auto &orbitPair = equivalentOrbitVector[idx];

    // Handle empty cluster
    if (orbitPair.first.size() == 1 && orbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitPair.first);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

double GetOrbitCorrelationFunction(const Config &config,
                                   BasisSet &atomicBasis,
                                   const vector<vector<size_t>> &orbitVector)
{

  double orbitClusterFunction = 0;

  for (const auto &latticeCluster : orbitVector)
  {
    vector<Element> elementVector;
    elementVector.reserve(latticeCluster.size());

    for (const auto &latticeId : latticeCluster)
    {
      auto element = config.GetElementOfLattice(latticeId);
      elementVector.emplace_back(element);
    }

    AtomClusterType atomClusterType(elementVector);

    double clusterFunction = atomicBasis.GetBasisProduct(atomClusterType);
    orbitClusterFunction += clusterFunction;
  }

  orbitClusterFunction = orbitClusterFunction / static_cast<double>(orbitVector.size());

  return orbitClusterFunction;
}

// Assiging element at a given targetLatticeId
double GetOrbitCorrelationFunction(const Config &config,
                                   BasisSet &atomicBasis,
                                   const size_t &targetLatticeId,
                                   const Element &elementToAssign,
                                   const vector<vector<size_t>> &orbitVector)
{

  double orbitClusterFunction = 0;

  for (const auto &latticeCluster : orbitVector)
  {
    vector<Element> elementVector;
    elementVector.reserve(latticeCluster.size());

    for (const auto &latticeId : latticeCluster)
    {
      // Assigning the element to a specific site
      if (latticeId == targetLatticeId)
      {
        elementVector.emplace_back(elementToAssign);
      }
      else
      {
        auto element = config.GetElementOfLattice(latticeId);
        elementVector.emplace_back(element);
      }
    }

    AtomClusterType atomClusterType(elementVector);

    double clusterFunction = atomicBasis.GetBasisProduct(atomClusterType);
    orbitClusterFunction += clusterFunction;
  }

  orbitClusterFunction = orbitClusterFunction / static_cast<double>(orbitVector.size());

  return orbitClusterFunction;
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
}