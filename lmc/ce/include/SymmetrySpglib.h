#ifndef LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_
#define LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_

#include "Config.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include <spglib.h>
#include "SymmetrySpglib.h"
#include "ClusterExpansion.h"
#include "PrintUtility.h"
#include <queue>
#include "UnionFind.h"

using namespace std;
using namespace Eigen;

// Returns space group symmetry operation based on the configuration
vector<pair<Matrix3d, Vector3d>> GetSymmetryOperations(const Config &config);

/*


Space group: R-3m (#166)
Point group: -3m
0 : {55, 243, 51, 198, 194, 6}
1 : {50, 1, 244, 5, 199, 248}
2 : {54, 70, 205, 179, 195, 228, 21, 201, 240, 9, 44, 48}
4 : {25, 224}
5 : {29, 45, 204, 24, 220, 225}
7 : {49, 200, 229, 20, 4, 245}

Encoded equivalent groups:
0 : {0, 1, 11, 26, 36, 37}
1 : {2, 12, 14, 23, 25, 35}
2 : {3, 4, 9, 10, 13, 16, 21, 24, 27, 28, 33, 34}
4 : {5, 32}
5 : {6, 7, 18, 19, 30, 31}
7 : {8, 15, 17, 20, 22, 29}
*/

vector<pair<Matrix3d, Vector3d>> GetSymmetryOperations(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const bool debug = false,
    const map<size_t, size_t> &latticeIdToIndexMap = map<size_t, size_t>{},
    const double symprec = 1e-5);

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
vector<set<vector<size_t>>> GetEquivalentClusters(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec = 1e-5,
    const bool debug = false);


=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
unordered_map<size_t, Eigen::RowVector3d> GetCenteredNeighborsAlongJumpDirection(
    const Config &config,
    const size_t maxBondOrder,
    const Vector3d &jumpDirection);

vector<size_t> GetCanonicalSortedSitesForPair(
    const Config &config, const pair<size_t, size_t> &latticeIdJumpPair,
    const size_t &maxBondOrder,
    const unordered_map<size_t, Eigen::RowVector3d> &referenceLatticeIdHashmap,
    const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations);

unordered_map<size_t, Eigen::RowVector3d> GetCenteredNeighboursSite(
    const Config &config,
    const size_t latticeId,
    const size_t maxBondOrder);

vector<size_t> GetCanonicalSortedSitesForSite(
    const Config &config,
    const size_t latticeId,
    const size_t &maxBondOrder);



/*
vector<set<vector<size_t>>> GetEquivalentClusters(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec = 1e-5,
    const bool debug = false);
*/
/*
// This is specifically for the KRA
// Where one need to provide the reference map about which the orbits encoding
// will be defined
vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEquivalentClustersEncoding(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxBondOrderOfCluster,
    const size_t &maxClusterSize,
    const unordered_map<size_t, Eigen::RowVector3d> &canonicalReferenceMap,
    const bool debug = false);

vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEquivalentClustersEncoding(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const size_t &maxBondOrderOfCluster,
    const bool debug = false,
    const double symprec = 1e-5);


vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEquivalentClustersEncoding(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const bool debug = false,
    const double symprec = 1e-5);
*/

#endif // LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_
