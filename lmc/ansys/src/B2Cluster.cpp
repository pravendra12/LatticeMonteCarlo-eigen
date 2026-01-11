#include "B2Cluster.h"

B2Cluster::B2Cluster(const Config &config) : config_(config)
{
  BuildB2Clusters();
}

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

vector<unordered_set<size_t>> B2Cluster::GetB2Clusters()
{
  return b2ClusterVector_;
}

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
