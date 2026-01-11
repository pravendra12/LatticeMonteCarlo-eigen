#ifndef LMC_ANSYS_INCLUDE_B2CLUSTER_H_
#define LMC_ANSYS_INCLUDE_B2CLUSTER_H_

#include "Config.h"
#include "B2OrderParameter.h"

using namespace std;

class B2Cluster
{
public:
  B2Cluster(const Config &config);

  
  // Returns the atom Ids of all the b2 cluster present in the configuration
  vector<unordered_set<size_t>> GetB2Clusters();
  
  
  
  private:
  void detectB2Clusters();
  bool growB2Cluster(size_t latticeId,
                     unordered_set<size_t> &cluster);

  unordered_set<size_t> mergeIfIntersect(const unordered_set<size_t> &set1,
                                         const unordered_set<size_t> &set2);
  void mergeAllClusters();

  const Config &config_;
  vector<unordered_set<size_t>> clusters_;
  unordered_set<size_t> visited_;
};

#endif // LMC_ANSYS_INCLUDE_B2CLUSTER_H_
