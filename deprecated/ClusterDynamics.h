#ifndef LMC_ANSYS_INCLUDE_CLUSTERDYNAMICS_H_
#define LMC_ANSYS_INCLUDE_CLUSTERDYNAMICS_H_

#include "Config.h"
#include "B2OrderParameter.h"

using namespace std;

class ClusterDynamics
{
public:
  ClusterDynamics(const Config &config);

  void detectB2Clusters(map<string, Config::VectorVariant> &auxiliaryList,
                        vector<int> &clusterSizeVector);

private:
  bool growB2Cluster(size_t latticeId,
                     unordered_set<size_t> &cluster);

  unordered_set<size_t> mergeIfIntersect(const unordered_set<size_t> &set1,
                                         const unordered_set<size_t> &set2);
  void mergeAllClusters();

  const Config &config_;
  vector<unordered_set<size_t>> clusters_;
  unordered_set<size_t> visited_;
};

#endif // LMC_ANSYS_INCLUDE_CLUSTERDYNAMICS_H_
