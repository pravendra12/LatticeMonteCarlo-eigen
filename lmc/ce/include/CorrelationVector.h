#ifndef LMC_CE_INCLUDE_CORRELATIONVECTOR_H_
#define LMC_CE_INCLUDE_CORRELATIONVECTOR_H_

#include <Eigen/Dense>
#include "Config.h"
#include "BasisSet.h"
#include "LatticeClusterType.hpp"
#include "CorrelationFunction.h"

using namespace std;
using namespace Eigen;

VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

// Return correlation vector with a given element assigned to a lattice Id
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

VectorXd GetCorrelationVector(
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
  const Config &config, 
  BasisSet &atomicBasis, 
  const vector<size_t> &canonicalSortedLatticeIds, 
  const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);
=======
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
    const Config &config,
    BasisSet &atomicBasis,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentOrbitVector);

double GetOrbitCorrelationFunction(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<vector<size_t>> &orbitVector);

double GetOrbitCorrelationFunction(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<vector<size_t>> &orbitVector);
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

#endif // LMC_CE_INCLUDE_CORRELATIONVECTOR_H_
