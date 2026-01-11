#include "B2Cluster.h"

B2Cluster::B2Cluster(const Config &config) : config_(config),
                                                         clusters_{},
                                                         visited_{}
{
  detectB2Clusters();
}

vector<unordered_set<size_t>> B2Cluster::GetB2Clusters()
{
  return clusters_;
}

void B2Cluster::detectB2Clusters()
{
  size_t numAtoms = config_.GetNumAtoms();
  for (size_t i = 0; i < numAtoms; ++i)
  {
    if (visited_.count(i))
      continue;

    unordered_set<size_t> cluster;
    if (growB2Cluster(i, cluster) && !cluster.empty())
    {
      clusters_.push_back(move(cluster));
    }
  }

  mergeAllClusters();
  
  /*
  size_t vacancyId = static_cast<size_t>(-1);
  
  try
  {
    vacancyId = config_.GetVacancyAtomId();
  }
  catch (const std::exception &e)
  {
    std::cerr << "Warning: No vacancy found, setting vacancyId = -1.\n";
  }

  vector<int> atomClusterTypeVector(numAtoms, -1);
  vector<int> atomClusterSizeVector(numAtoms, -1);

  int clusterId = 1;

  for (const auto &cluster : clusters_)
  {
    int clusterSize = cluster.size();

    // If vacancy is part of the cluster, we exclude it from the count
    if (cluster.find(vacancyId) != cluster.end())
    {
      clusterSize -= 1;
    }

    for (auto atomId : cluster)
    {
      atomClusterTypeVector[atomId] = clusterId;
      atomClusterSizeVector[atomId] = clusterSize;
    }

    clusterSizeVector.emplace_back(clusterSize);

    clusterId++;
  }

  auxiliaryList["clusterType"] = atomClusterTypeVector;
  auxiliaryList["clusterSize"] = atomClusterSizeVector;
  */
}



bool B2Cluster::growB2Cluster(size_t atomId,
                                    unordered_set<size_t> &cluster)
{
  if (!visited_.emplace(atomId).second)
    return false;

  bool isB2Center = B2OrderParameter::isB2Ordered(config_, atomId);

  if (!isB2Center)
    return false;

  cluster.insert(atomId);
  // First nearest neighbours
  const auto neighbors = config_.GetNeighborAtomIdVectorOfAtom(atomId, 1);

  for (size_t nId : neighbors)
  {
    cluster.insert(nId); // Include all neighbors regardless of order
    if (B2OrderParameter::isB2Ordered(config_, nId))
    {
      growB2Cluster(nId, cluster); // Recurse on ordered neighbors only
    }
  }
  return true;
}

unordered_set<size_t> B2Cluster::mergeIfIntersect(
    const unordered_set<size_t> &set1,
    const unordered_set<size_t> &set2)
{

  for (const size_t &val : set1)
  {
    if (set2.count(val))
    {
      unordered_set<size_t> result = set1;
      result.insert(set2.begin(), set2.end());
      return result;
    }
  }
  return {};
}

void B2Cluster::mergeAllClusters()
{
  for (size_t i = 0; i < clusters_.size(); ++i)
  {
    for (size_t j = i + 1; j < clusters_.size();)
    {
      auto combined = mergeIfIntersect(clusters_[i], clusters_[j]);
      if (!combined.empty())
      {
        clusters_.erase(clusters_.begin() + j);
        clusters_.erase(clusters_.begin() + i);
        clusters_.push_back(combined);
        i = static_cast<size_t>(-1);
        break;
      }
      else
      {
        ++j;
      }
    }
  }
}
