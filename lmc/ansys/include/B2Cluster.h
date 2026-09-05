#ifndef LMC_ANSYS_INCLUDE_B2CLUSTER_H_
#define LMC_ANSYS_INCLUDE_B2CLUSTER_H_

#include "Config.h"
#include <queue>
#include <unordered_set>

using namespace std;

class B2Cluster
{
public:
  B2Cluster(const Config &config);

  void WriteB2ClusterConfig(const string &filename);
  vector<unordered_set<size_t>> GetB2Clusters();

  /**
   * @brief Writes the atom IDs belonging to each B2 cluster.
   *
   * Each line has the format:
   *
   * clusterId: atomId atomId atomId ...
   *
   * An atom ID may appear in more than one cluster if the atom is shared
   * between overlapping B2 clusters.
   *
   * The output is compressed when the filename ends in .gz or .bz2.
   *
   * @param filename Name of the output file.
   */
  void WriteB2ClusterAtomIds(const string &filename) const;

private:
  void BuildB2Clusters();
  bool isB2(const size_t latticeId);

  const Config &config_;
  vector<unordered_set<size_t>> b2ClusterVector_{};
};

#endif // LMC_ANSYS_INCLUDE_B2CLUSTER_H_
