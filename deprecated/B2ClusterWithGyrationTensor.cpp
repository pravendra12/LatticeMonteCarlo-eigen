#include "B2Cluster.h"

B2Cluster::B2Cluster(const Config &config) : config_(config)
{
  BuildB2Clusters();
}

/*
void B2Cluster::WriteB2ClusterConfig(const string &filename)
{
  size_t numSites = config_.GetNumLattices();

  // Auxiliary lists
  Config::VectorVariant clusterIdVec = vector<int>(numSites, -1);
  Config::VectorVariant clusterSizeVec = vector<size_t>(numSites, 0);
  Config::VectorVariant isSharedVec = vector<int>(numSites, 0); // 0 = not shared, 1 = shared

  unordered_map<size_t, size_t> atomClusterCount;

  // First pass: count how many clusters each atom belongs to
  for (size_t clusterId = 0; clusterId < b2ClusterVector_.size(); ++clusterId)
  {
    for (size_t latticeId : b2ClusterVector_[clusterId])
    {
      auto atomId = config_.GetAtomIdOfLattice(latticeId);
      atomClusterCount[atomId]++;
    }
  }

  // Second pass: fill auxiliary lists
  for (size_t clusterId = 0; clusterId < b2ClusterVector_.size(); ++clusterId)
  {
    size_t clusterSize = b2ClusterVector_[clusterId].size();
    for (size_t latticeId : b2ClusterVector_[clusterId])
    {
      auto atomId = config_.GetAtomIdOfLattice(latticeId);
      get<vector<int>>(clusterIdVec)[atomId] = int(clusterId);
      get<vector<size_t>>(clusterSizeVec)[atomId] = clusterSize;

      if (atomClusterCount[atomId] > 1)
        get<vector<int>>(isSharedVec)[atomId] = 1; // mark as shared
    }
  }

  map<string, Config::VectorVariant> auxiliaryLists;
  auxiliaryLists["clusterId"] = clusterIdVec;
  auxiliaryLists["clusterSize"] = clusterSizeVec;
  auxiliaryLists["sharedAtom"] = isSharedVec; // new shared flag

  map<string, Config::ValueVariant> globalList; // empty

  Config::WriteXyzExtended(filename, config_, auxiliaryLists, globalList);
}
*/


void B2Cluster::WriteB2ClusterConfig(const string &filename)
{
  size_t numSites = config_.GetNumLattices();

  // Auxiliary lists
  Config::VectorVariant clusterIdVec = vector<int>(numSites, -1);
  Config::VectorVariant clusterSizeVec = vector<size_t>(numSites, 0);
  Config::VectorVariant isSharedVec = vector<int>(numSites, 0); // 0 = not shared, 1 = shared

  // Gyration tensor components (symmetric, so only the unique entries are stored)
  Config::VectorVariant gyrationS11Vec = vector<double>(numSites, 0.0);
  Config::VectorVariant gyrationS12Vec = vector<double>(numSites, 0.0);
  Config::VectorVariant gyrationS13Vec = vector<double>(numSites, 0.0);
  Config::VectorVariant gyrationS22Vec = vector<double>(numSites, 0.0);
  Config::VectorVariant gyrationS23Vec = vector<double>(numSites, 0.0);
  Config::VectorVariant gyrationS33Vec = vector<double>(numSites, 0.0);

  unordered_map<size_t, size_t> atomClusterCount;

  // First pass: count how many clusters each atom belongs to
  for (size_t clusterId = 0; clusterId < b2ClusterVector_.size(); ++clusterId)
  {
    for (size_t latticeId : b2ClusterVector_[clusterId])
    {
      auto atomId = config_.GetAtomIdOfLattice(latticeId);
      atomClusterCount[atomId]++;
    }
  }

  // Second pass: fill auxiliary lists
  for (size_t clusterId = 0; clusterId < b2ClusterVector_.size(); ++clusterId)
  {
    size_t clusterSize = b2ClusterVector_[clusterId].size();
    const Matrix3d &gyrationTensor = b2ClusterGyrationTensorVector_[clusterId];

    for (size_t latticeId : b2ClusterVector_[clusterId])
    {
      auto atomId = config_.GetAtomIdOfLattice(latticeId);
      get<vector<int>>(clusterIdVec)[atomId] = int(clusterId);
      get<vector<size_t>>(clusterSizeVec)[atomId] = clusterSize;

      if (atomClusterCount[atomId] > 1)
        get<vector<int>>(isSharedVec)[atomId] = 1; // mark as shared

      // same gyration tensor for every atom in this cluster
      get<vector<double>>(gyrationS11Vec)[atomId] = gyrationTensor(0, 0);
      get<vector<double>>(gyrationS12Vec)[atomId] = gyrationTensor(0, 1);
      get<vector<double>>(gyrationS13Vec)[atomId] = gyrationTensor(0, 2);
      get<vector<double>>(gyrationS22Vec)[atomId] = gyrationTensor(1, 1);
      get<vector<double>>(gyrationS23Vec)[atomId] = gyrationTensor(1, 2);
      get<vector<double>>(gyrationS33Vec)[atomId] = gyrationTensor(2, 2);
    }
  }

  map<string, Config::VectorVariant> auxiliaryLists;
  auxiliaryLists["clusterId"] = clusterIdVec;
  auxiliaryLists["clusterSize"] = clusterSizeVec;
  auxiliaryLists["sharedAtom"] = isSharedVec; // new shared flag

  // gyration tensor components (unique entries only, since Sij == Sji)
  auxiliaryLists["S11"] = gyrationS11Vec;
  auxiliaryLists["S12"] = gyrationS12Vec;
  auxiliaryLists["S13"] = gyrationS13Vec;
  auxiliaryLists["S22"] = gyrationS22Vec;
  auxiliaryLists["S23"] = gyrationS23Vec;
  auxiliaryLists["S33"] = gyrationS33Vec;

  map<string, Config::ValueVariant> globalList; // empty

  Config::WriteXyzExtended(filename, config_, auxiliaryLists, globalList);
}


vector<unordered_set<size_t>> B2Cluster::GetB2Clusters()
{
  return b2ClusterVector_;
}

vector<Matrix3d> B2Cluster::GetGyrationTensors()
{
  return b2ClusterGyrationTensorVector_;
}

void B2Cluster::BuildB2Clusters()
{
  const size_t numSites = config_.GetNumLattices();

  // Convert relative displacement vectors to Cartesian displacement vectors
  const Matrix3d basisTranspose =
      config_.GetBasis().transpose();

  // Clear previously calculated cluster information
  b2ClusterVector_.clear();
  b2ClusterGyrationTensorVector_.clear();

  unordered_set<size_t> visitedB2; // only B2 sites marked visited

  for (size_t latticeId = 0; latticeId < numSites; ++latticeId)
  {
    // skip if site already processed as B2 center
    if (visitedB2.count(latticeId))
      continue;

    // check if this site is a B2 center
    if (!isB2(latticeId))
      continue;

    // start new cluster
    unordered_set<size_t> currentCluster;
    queue<size_t> bfsQueue;

    /*
     * Local data used only to calculate the gyration tensor.
     *
     * currentPositionIndex maps each lattice ID to the corresponding
     * position in currentUnwrappedCartesianPositions.
     */
    unordered_map<size_t, size_t> currentPositionIndex;
    vector<Vector3d> currentUnwrappedCartesianPositions;

    // use the first B2 center as the unwrapping anchor
    currentPositionIndex.emplace(latticeId, 0);
    currentUnwrappedCartesianPositions.push_back(
        config_.GetCartesianPositionOfLattice(latticeId));

    bfsQueue.push(latticeId);
    visitedB2.insert(latticeId);

    while (!bfsQueue.empty())
    {
      size_t center = bfsQueue.front();
      bfsQueue.pop();

      // add center to cluster
      currentCluster.insert(center);

      // get 1NN + 2NN neighbors
      auto nnUpToSecond = config_.GetNeighborLatticeIdsUpToOrder(center, 2);

      const size_t centerPositionIndex =
          currentPositionIndex.at(center);

      /*
       * Make a copy rather than a reference because pushing new positions
       * into the vector can reallocate its storage.
       */
      const Vector3d centerUnwrappedPosition =
          currentUnwrappedCartesianPositions[centerPositionIndex];

      /*
       * Single pass over the neighbor list: add every neighbor (B2 or not)
       * to the cluster, expand the BFS through unvisited B2 sites, and
       * unwrap every neighbor's Cartesian position.
       */
      for (size_t nnId : nnUpToSecond)
      {
        // add all NN (B2 or not) into the cluster
        currentCluster.insert(nnId);

        // BFS expansion: only through real, not-yet-queued B2 sites
        if (isB2(nnId) && !visitedB2.count(nnId))
        {
          visitedB2.insert(nnId);
          bfsQueue.push(nnId);
        }

        // minimum-image displacement from center to nnId
        const Vector3d relativeDisplacement =
            config_.GetRelativeDistanceVectorLattice(center, nnId);

        // convert relative displacement to Cartesian displacement
        const Vector3d cartesianDisplacement =
            basisTranspose * relativeDisplacement;

        const Vector3d candidateUnwrappedPosition =
            centerUnwrappedPosition + cartesianDisplacement;

        /*
         * Insert the position only when this lattice site is encountered
         * for the first time in the current cluster.
         */
        const auto [positionIterator, positionInserted] =
            currentPositionIndex.emplace(
                nnId, currentUnwrappedCartesianPositions.size());

        if (positionInserted)
        {
          currentUnwrappedCartesianPositions.push_back(
              candidateUnwrappedPosition);
        }
        else
        {
          /*
           * A neighboring site can be reached from more than one B2 center.
           * Both paths should produce the same unwrapped position.
           */
          const Vector3d &existingUnwrappedPosition =
              currentUnwrappedCartesianPositions[positionIterator->second];

          const double positionDifference =
              (existingUnwrappedPosition - candidateUnwrappedPosition).norm();

          if (positionDifference > constants::kEpsilon)
          {
            throw runtime_error(
                "B2 cluster cannot be unwrapped consistently. "
                "The cluster may percolate through the "
                "periodic simulation cell.");
          }
        }
      }
    }

    // calculate the geometry using the temporary unwrapped positions
    const ClusterGeometry clusterGeometry(currentUnwrappedCartesianPositions);

    // store the original B2 cluster
    b2ClusterVector_.push_back(std::move(currentCluster));

    // store only the gyration tensor
    b2ClusterGyrationTensorVector_.push_back(
        clusterGeometry.GetGyrationTensor());
  }
}

/*
void B2Cluster::BuildB2Clusters()
{
  const size_t numSites = config_.GetNumLattices();
  unordered_set<size_t> visitedB2; // only B2 sites marked visited

  for (size_t latticeId = 0; latticeId < numSites; ++latticeId)
  {
    // skip if site already processed as B2 center
    if (visitedB2.count(latticeId))
      continue;

    // check if this site is a B2 center
    if (!isB2(latticeId))
      continue;

    // start new cluster
    unordered_set<size_t> currentCluster;
    queue<size_t> bfsQueue;

    bfsQueue.push(latticeId);
    visitedB2.insert(latticeId);

    while (!bfsQueue.empty())
    {
      size_t center = bfsQueue.front();
      bfsQueue.pop();

      // add center to cluster
      currentCluster.insert(center);

      // get 1NN + 2NN neighbors
      auto nnUpToSecond = config_.GetNeighborLatticeIdsUpToOrder(center, 2);

      // add all NN (B2 or not) into the cluster
      for (size_t nnId : nnUpToSecond)
      {
        currentCluster.insert(nnId);
      }

      // BFS expansion: only through real B2 sites
      for (size_t nnId : nnUpToSecond)
      {
        if (visitedB2.count(nnId))
          continue;

        if (isB2(nnId))
        {
          visitedB2.insert(nnId);
          bfsQueue.push(nnId);
        }
      }
    }

    b2ClusterVector_.push_back(std::move(currentCluster));
  }
}
*/

bool B2Cluster::isB2(const size_t latticeId)
{
  const Element centralElement = config_.GetElementOfLattice(latticeId);

  if (centralElement == Element("X"))
    return false;

  const auto firstNN = config_.GetNeighborLatticeIdVectorOfLattice(latticeId, 1);

  // The 1NN element type must be uniform and opposite to central
  const Element firstNNElement = config_.GetElementOfLattice(firstNN[0]);

  if (firstNNElement == centralElement)
    return false;

  for (size_t nnId : firstNN)
  {
    if (config_.GetElementOfLattice(nnId) != firstNNElement)
      return false;
  }

  const auto secondNN = config_.GetNeighborLatticeIdVectorOfLattice(latticeId, 2);

  for (size_t nnId : secondNN)
  {
    if (config_.GetElementOfLattice(nnId) != centralElement)
      return false;
  }
  return true;
}
