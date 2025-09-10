#include "SymmetrySpglib.h"

// Returns space group operation for a passed lattice type
vector<pair<Matrix3d, Vector3d>> GetSymmetryOperations(
    const Config &config)
{
  // First NN
  auto nnIdsVector = config.GetNeighborLatticeIdVectorOfLattice(0, 1);

  unordered_set<size_t> nnIdsSet(nnIdsVector.begin(), nnIdsVector.end());

  return GetSymmetryOperations(config, nnIdsSet, true);
}

vector<pair<Matrix3d, Vector3d>> GetSymmetryOperations(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const bool debug,
    const map<size_t, size_t> &latticeIdToIndexMap,
    const double symprec)
{
  double lattice[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  int numAtoms = latticeIdSet.size();

  double (*positions)[3] = new double[numAtoms][3];
  int *types = new int[numAtoms];

  int i = 0;
  for (const auto &latticeId : latticeIdSet)
  {
    VectorXd relativePosition = config.GetRelativePositionOfLattice(latticeId);

    positions[i][0] = relativePosition(0);
    positions[i][1] = relativePosition(1);
    positions[i][2] = relativePosition(2);

    types[i] = 1;
    i += 1;
  }

  SpglibDataset *dataset =
      spg_get_dataset(lattice, positions, types, numAtoms, symprec);

  if (!dataset)
  {
    cerr << "Error in `GetSymmetryOperations`: Failed to get spglib dataset."
         << endl;
    exit(1);
  }

  if (debug)
  {
    cout << "Space group: " << dataset->international_symbol << " (#"
         << dataset->spacegroup_number << ")" << endl;

    cout << "Point group: " << dataset->pointgroup_symbol << endl;
  }

  // Symmetry Operations
  vector<pair<Matrix3d, Vector3d>> symmetryOperations;
  for (int j = 0; j < dataset->n_operations; ++j)
  {
    Matrix3d rot;
    Vector3d trans;

    for (int r = 0; r < 3; ++r)
    {
      trans(r) = dataset->translations[j][r];
      for (int c = 0; c < 3; ++c)
      {
        rot(r, c) = dataset->rotations[j][r][c];
      }
    }
    pair<Matrix3d, Vector3d> symOpPair = {rot, trans};
    symmetryOperations.emplace_back(symOpPair);
  }

  // Print the equivalent group under the symmetry operations
  if (debug)
  {
    cout << "\n Equivalent Groups: " << endl;
    map<size_t, vector<size_t>> equivalentSitesMap;

    i = 0;
    for (const auto &latticeId : latticeIdSet)
    {
      size_t key = dataset->equivalent_atoms[i];

      if (equivalentSitesMap.find(key) != equivalentSitesMap.end())
      {
        equivalentSitesMap[key].push_back(latticeId);
      }
      else
      {
        equivalentSitesMap[key] = {latticeId};
      }
      i += 1;
    }

    for (const auto &group : equivalentSitesMap)
    {
      cout << group.first << " : ";
      print1DVector(group.second);
    }

    if (!latticeIdToIndexMap.empty())
    {
      cout << "Encoded equivalent groups: " << endl;
      for (const auto &group : equivalentSitesMap)
      {
        vector<size_t> indices;
        indices.reserve(group.second.size());
        for (const auto &id : group.second)
          indices.push_back(latticeIdToIndexMap.at(id));

        sort(indices.begin(), indices.end());

        cout << group.first << " : ";
        print1DVector(indices);
      }
    }
  }

  spg_free_dataset(dataset);
  delete[] positions;

  return symmetryOperations;
}

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
// Helper function to get the equivalent clusters
static vector<set<vector<size_t>>> GetEquivalentGroups(
    const map<vector<size_t>, set<vector<size_t>>> &equivalentMap,
    map<vector<size_t>, int> &clustersToGroupMap)
{
  set<vector<size_t>> visited;
  vector<set<vector<size_t>>> groups;

  for (const auto &kv : equivalentMap)
  {
    const vector<size_t> &start = kv.first;
    if (visited.find(start) != visited.end())
      continue;

    set<vector<size_t>> cls;
    queue<vector<size_t>> q;
    q.push(start);

    while (!q.empty())
    {
      vector<size_t> cur = q.front();
      q.pop();
      if (visited.find(cur) != visited.end())
        continue;

      visited.insert(cur);
      cls.insert(cur);

      // forward neighbors
      auto it = equivalentMap.find(cur);
      if (it != equivalentMap.end())
      {
        for (const auto &nbr : it->second)
          if (visited.find(nbr) == visited.end())
            q.push(nbr);
      }

      // reverse neighbors: other keys that list cur as equivalent
      for (const auto &other_kv : equivalentMap)
      {
        const vector<size_t> &other_key = other_kv.first;
        if (visited.find(other_key) != visited.end())
          continue;
        const auto &other_set = other_kv.second;
        if (other_set.find(cur) != other_set.end())
          q.push(other_key);
      }
    }

    groups.push_back(cls);
  }

  // produce mapping triplet -> classId
  clustersToGroupMap.clear();
  for (size_t cid = 0; cid < groups.size(); ++cid)
  {
    for (const auto &t : groups[cid])
    {
      clustersToGroupMap[t] = static_cast<int>(cid);
    }
  }

  return groups;
}
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
/*
vector<set<vector<size_t>>> GetEquivalentClusters(
    const Config &config, const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>>
        &latticeClusterSet,
    const double symprec, const bool debug)
{

  auto symmetryOperations =
      GetSymmetryOperations(config, latticeIdSet, debug, {}, symprec);
  // Map latticeId to Position
  unordered_map<size_t, Vector3d> latticeIdToPositionMap;

  // Map Position to latticeId
  // String key inplace or positions
  unordered_map<Vector3d, size_t, Vector3dHash2, Vector3dEqual>
      positionToLatticeIdMap;

  for (auto latticeId : latticeIdSet)
  {
    Vector3d position = config.GetRelativePositionOfLattice(latticeId);
    latticeIdToPositionMap[latticeId] = position;
    positionToLatticeIdMap[position] = latticeId;
  }

  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (const auto &cluster : latticeClusterSet)
  {
    set<vector<size_t>> equivalentSet;
    set<LatticeCluster> equivalentLatticeClusters;

    auto originalCluster = cluster.GetLatticeIdVector();

    // if (originalCluster.size() <= 1)
    // {
    //   // Skip the singlets and empty cluster for now;
    //   continue;
    // }

    // Iterate over the symmetry operations to get the equivalent clusters

    for (const auto &operation : symmetryOperations)
    {
      vector<size_t> transformedCluster;
      transformedCluster.reserve(originalCluster.size());
      bool allFound = true;

      // Tranform the cluster
      for (auto latticeId : originalCluster)
      {
        Vector3d position = latticeIdToPositionMap[latticeId];

        // operation = [rotation, translation]
        Vector3d transformedPosition =
            operation.first * position + operation.second;

        // wrap to [0, 1)
        for (int k = 0; k < 3; ++k)
        {
          transformedPosition[k] =
              transformedPosition[k] - floor(transformedPosition[k]);
          if (transformedPosition[k] < 0)
            transformedPosition[k] += 1.0;
        }

        // Get the lattice Id corresponding to the tranformed position
        auto it = positionToLatticeIdMap.find(transformedPosition);
        if (it != positionToLatticeIdMap.end())
        {
          transformedCluster.push_back(it->second);
        }
        // May be need to add other conditions
        else
        {
          allFound = false;
        }
      }

      if (!allFound)
        continue;

      auto latticeClusterType =
          IdentifyLatticeClusterType(config, transformedCluster);
      LatticeCluster tranformedLatticeCluster(latticeClusterType,
                                              transformedCluster);

      equivalentSet.insert(tranformedLatticeCluster.GetLatticeIdVector());
    }

    equivalentMap[originalCluster] = equivalentSet;
  }

  map<vector<size_t>, int> clustersToGroupMap;

  auto equivalentGroups =
      GetEquivalentGroups(equivalentMap, clustersToGroupMap);

  if (debug)
  {

    cout << "Found " << equivalentGroups.size() << " equivalence classes.\n\n";

    for (size_t cid = 0; cid < equivalentGroups.size(); ++cid)
    {
      cout << "Equivalence Class " << cid << ":\n";
      for (const auto &t : equivalentGroups[cid])
      {
        cout << " Cluster: ";
        for (auto id : t)
          cout << id << " ";

        auto clusterType = IdentifyLatticeClusterType(config, t);
        cout << clusterType << "\n";
      }
      cout << "\n";
    }
  }

  return equivalentGroups;
}
  */
/*
vector<set<vector<size_t>>> GetEquivalentClusters(
    const Config &config, const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>>
        &latticeClusterSet,
    const double symprec, const bool debug)
{

  auto symmetryOperations =
      GetSymmetryOperations(config, latticeIdSet);

  // Map latticeId to Position (fractional)
  unordered_map<size_t, Vector3d> latticeIdToPositionMap;

  // No need for positionToLatticeIdMap anymore, since we'll search id_to_pos

  for (auto latticeId : latticeIdSet)
  {
    Vector3d position = config.GetRelativePositionOfLattice(latticeId);
    latticeIdToPositionMap[latticeId] = position;
  }

  // Helper lambda to find the closest lattice ID to a given fractional coordinate (targetFractionalCoordinate)
  // within a specified symmetry tolerance (symprec), accounting for periodic boundary conditions.
  //
  // Steps:
  // 1. Iterate over all lattice sites in latticeIdToPositionMap (fractional coordinates).
  // 2. Compute fractional difference between targetFractionalCoordinate and the candidate lattice site.
  // 3. Apply periodic wrapping by subtracting the nearest integer (round) so the difference
  //    lies within [-0.5, 0.5] for each dimension (handles supercell periodicity).
  // 4. Convert the wrapped fractional difference into Cartesian coordinates using the basis matrix.
  // 5. Compute squared Cartesian distance to avoid unnecessary sqrt until the end.
  // 6. Keep track of the lattice ID with the smallest distance.
  // 7. Return true if the minimum distance is less than the symmetry tolerance (symprec),
  //    and update matchedLatticeId with the closest lattice ID.

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId) -> bool
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
=======

  auto FindClosestLatticeId = [&](const Vector3d &targetFractionalCoordinate, size_t &matchedLatticeId, double tol = 1e-5) -> bool
>>>>>>> Stashed changes
  {
    const Matrix3d basis = config.GetBasis();
    double minDistance = numeric_limits<double>::max();
    matchedLatticeId = static_cast<size_t>(-1);

    for (const auto &pair : latticeIdToPositionMap)
    {
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
      Vector3d diffFractional = pair.second - targetFractionalCoordinate;
      for (int k = 0; k < 3; ++k)
      {
        diffFractional[k] -= round(diffFractional[k]);
      }
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
=======
>>>>>>> Stashed changes
      // Fractional difference
      Vector3d diffFractional = pair.second - targetFractionalCoordinate;

      // Apply periodic wrapping [-0.5, 0.5)
      for (int k = 0; k < 3; ++k)
        diffFractional[k] -= round(diffFractional[k]);

      // Cartesian distance squared
>>>>>>> Stashed changes
      Vector3d diffCartesian = basis * diffFractional;
      double dist_sq = diffCartesian.squaredNorm();

      if (dist_sq < minDistance)
      {
        minDistance = dist_sq;
        matchedLatticeId = pair.first;
      }
    }

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
    return (sqrt(minDistance) < symprec);
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
=======
    // If minDistance is below tolerance (in Å^2), consider it a match
    return (sqrt(minDistance) < tol);
>>>>>>> Stashed changes
  };

  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (const auto &cluster : latticeClusterSet)
  {
    set<vector<size_t>> equivalentSet;
    set<LatticeCluster> equivalentLatticeClusters;

    auto originalCluster = cluster.GetLatticeIdVector();

    // if (originalCluster.size() <= 1)
    // {
    //   // Skip the singlets and empty cluster for now;
    //   continue;
    // }

    // Iterate over the symmetry operations to get the equivalent clusters

    for (const auto &operation : symmetryOperations)
    {
      vector<size_t> transformedCluster;
      transformedCluster.reserve(originalCluster.size());
      bool allFound = true;

      // Transform the cluster
      for (auto latticeId : originalCluster)
      {
        Vector3d position = latticeIdToPositionMap[latticeId];

        // operation = [rotation, translation]
        Vector3d transformedPosition =
            operation.first * position + operation.second;

        // wrap to [0, 1)
        for (int k = 0; k < 3; ++k)
        {
          transformedPosition[k] =
              transformedPosition[k] - floor(transformedPosition[k]);
          if (transformedPosition[k] < 0)
            transformedPosition[k] += 1.0;
          // Additional boundary snap for stability (optional, tune epsilon if needed)
          const double eps = 1e-10;
          if (transformedPosition[k] > 1.0 - eps)
            transformedPosition[k] = 0.0;
        }

        size_t matchedId;
        if (FindClosestLatticeId(transformedPosition, matchedId))
        {
          transformedCluster.push_back(matchedId);
        }
        else
        {
          allFound = false;
          break; // Early exit if any site fails
        }
      }

      if (!allFound)
        continue;

      auto latticeClusterType =
          IdentifyLatticeClusterType(config, transformedCluster);
      LatticeCluster tranformedLatticeCluster(latticeClusterType,
                                              transformedCluster);

      equivalentSet.insert(tranformedLatticeCluster.GetLatticeIdVector());
    }

    equivalentMap[originalCluster] = equivalentSet;
  }

  map<vector<size_t>, int> clustersToGroupMap;

  auto equivalentGroups =
      GetEquivalentGroups(equivalentMap, clustersToGroupMap);

  if (debug)
  {
    cout << "Found " << equivalentGroups.size() << " equivalence classes.\n\n";

    for (size_t cid = 0; cid < equivalentGroups.size(); ++cid)
    {
      cout << "Equivalence Class " << cid << ":\n";
      for (const auto &t : equivalentGroups[cid])
      {
        cout << " Cluster: ";
        for (auto id : t)
          cout << id << " ";

        auto clusterType = IdentifyLatticeClusterType(config, t);
        cout << clusterType << "\n";
      }
      cout << "\n";
    }
  }

  return equivalentGroups;
}
*/

// Canonical LCE
// Choose one axis as representative
// Then for other axis apply symmetry operations to get the sorted atoms
// relative to the canonical representation

static bool areVectorsClose(const Vector3d &a, const Vector3d &b,
                            double tol = 1e-6)
{
  return (a - b).norm() < tol;
}

static bool ComparePositionMaps(const unordered_map<size_t, RowVector3d> &map1,
                                const unordered_map<size_t, RowVector3d> &map2,
                                double tol = 1e-6)
{
  if (map1.size() != map2.size())
    return false;

  // Extract positions into vectors
  vector<RowVector3d> vec1, vec2;
  for (const auto &[id, pos] : map1)
    vec1.push_back(pos);
  for (const auto &[id, pos] : map2)
    vec2.push_back(pos);

  // Custom comparator for sorting Vector3d lex order
  auto vectorLess = [](const RowVector3d &a, const RowVector3d &b)
  {
    if (abs(a.x() - b.x()) > 1e-8)
      return a.x() < b.x();
    if (abs(a.y() - b.y()) > 1e-8)
      return a.y() < b.y();
    return a.z() < b.z();
  };

  sort(vec1.begin(), vec1.end(), vectorLess);
  sort(vec2.begin(), vec2.end(), vectorLess);

  // Compare element-wise with tolerance
  for (size_t i = 0; i < vec1.size(); ++i)
  {
    if (!areVectorsClose(vec1[i], vec2[i], tol))
      return false;
  }

  return true;
}

unordered_map<size_t, Eigen::RowVector3d>
GetCenteredNeighours(const Config &config,
                     const pair<size_t, size_t> &latticeIdJumpPair,
                     const size_t maxBondOrder)
{
  auto neighbouringLatticeIds =
      config.GetNeighboringLatticeIdSetOfPair(latticeIdJumpPair, maxBondOrder);

  size_t numSites = neighbouringLatticeIds.size();
  Eigen::RowVector3d moveToCenter =
      Eigen::RowVector3d(0.5, 0.5, 0.5) -
      config.GetLatticePairCenter(latticeIdJumpPair);

  unordered_map<size_t, Eigen::RowVector3d> latticeIdHashmap;
  latticeIdHashmap.reserve(numSites);

  // Move lattice IDs to center
  for (const auto id : neighbouringLatticeIds)
  {
    Eigen::RowVector3d relativePosition =
        config.GetRelativePositionOfLattice(id).transpose();
    relativePosition += moveToCenter;
    relativePosition -=
        relativePosition.unaryExpr([](double x)
                                   { return floor(x); });
    latticeIdHashmap.emplace(id, relativePosition);
  }

  return latticeIdHashmap;
}

unordered_map<size_t, Eigen::RowVector3d>
GetCenteredNeighborsAlongJumpDirection(const Config &config,
                                       const size_t maxBondOrder,
                                       const Vector3d &jumpDirection)
{
  Vector3d normJumpDir = jumpDirection.normalized();

  pair<size_t, size_t> latticeIdJumpPair;

  bool isDirectionVectorFound = false;
  for (auto id : config.GetNeighborLatticeIdVectorOfLattice(0, 1))
  {
    Vector3d directionVector =
        config.GetRelativeDistanceVectorLattice(0, id).normalized();

    if (areVectorsClose(directionVector, normJumpDir))
    {
      latticeIdJumpPair = {0, id};
      isDirectionVectorFound = true;
      cout << "Along jump direction <" << jumpDirection.transpose() << ">, "
           << "the lattice ID pair is {" << latticeIdJumpPair.first << ", "
           << latticeIdJumpPair.second << "}" << endl;
      break;
    }
  }

  if (!isDirectionVectorFound)
  {
    cerr << "Error: No neighbor lattice vector found matching jump direction "
         << jumpDirection.transpose() << endl;
    return {};
  }

  return GetCenteredNeighours(config, latticeIdJumpPair, maxBondOrder);
}

inline bool PositionCompareState(const pair<size_t, Eigen::RowVector3d> &lhs,
                                 const pair<size_t, Eigen::RowVector3d> &rhs)
{
  const auto &relative_position_lhs = lhs.second;
  const auto &relative_position_rhs = rhs.second;

  // Compare individual components (x, y, z)
  const double diff_x = relative_position_lhs[0] - relative_position_rhs[0];
  if (diff_x < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_x > constants::kEpsilon)
  {
    return false;
  }

  const double diff_y = relative_position_lhs[1] - relative_position_rhs[1];
  if (diff_y < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_y > constants::kEpsilon)
  {
    return false;
  }

  const double diff_z = relative_position_lhs[2] - relative_position_rhs[2];
  if (diff_z < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_z > constants::kEpsilon)
  {
    return false;
  }

  return false;
}

// symmetryOperations : including the central site
// 0 and its neighours
vector<size_t> GetCanonicalSortedSitesForPair(
    const Config &config, const pair<size_t, size_t> &latticeIdJumpPair,
    const size_t &maxBondOrder,
    const unordered_map<size_t, Eigen::RowVector3d> &referenceLatticeIdHashmap,
    const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations)
{

  auto centeredLatticeIdHashmap =
      GetCenteredNeighours(config, latticeIdJumpPair, maxBondOrder);

  unordered_map<size_t, Eigen::RowVector3d> matchedTransformedPositionMap;

  for (auto &operation : symmetryOperations)
  {
    unordered_map<size_t, Eigen::RowVector3d> transformedPositionMap;

    for (auto entry : centeredLatticeIdHashmap)
    {
      Vector3d position = entry.second;

      // operation = [rotation, translation]
      Vector3d transformedPosition =
          operation.first * position + operation.second;

      // wrap to [0, 1)
      for (int k = 0; k < 3; ++k)
      {
        transformedPosition[k] =
            transformedPosition[k] - floor(transformedPosition[k]);
        if (transformedPosition[k] < 0)
          transformedPosition[k] += 1.0;
      }

      transformedPositionMap[entry.first] = transformedPosition;
    }
    auto areSame =
        ComparePositionMaps(transformedPositionMap, referenceLatticeIdHashmap);
    if (areSame)
    {
      matchedTransformedPositionMap = transformedPositionMap;
      break;
    }
  }

  if (matchedTransformedPositionMap.empty())
  {
    cerr << "Error in `GetCanonicalSortedSitesForPair`: No equivalent lattice "
            "ID hashmap exists for the given symmetry operations."
         << endl;
    return {};
  }

  // Sort the transformed positions
  // Convert unordered_map to vector for sorting
  vector<pair<size_t, Eigen::RowVector3d>> latticeIdVector(
      matchedTransformedPositionMap.begin(),
      matchedTransformedPositionMap.end());

  // Sort the lattice vector based on PositionCompare
  sort(latticeIdVector.begin(), latticeIdVector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  vector<size_t> sortedLatticeIdVector;
  sortedLatticeIdVector.reserve(latticeIdVector.size());
  for (const auto &pair : latticeIdVector)
  {
    sortedLatticeIdVector.push_back(pair.first);
  }

  return sortedLatticeIdVector;
}

unordered_map<size_t, Eigen::RowVector3d>
GetCenteredNeighboursSite(const Config &config,
                          size_t latticeId,
                          size_t maxBondOrder)
{
  auto neighbourIds = config.GetNeighborLatticeIdsUpToOrder(latticeId, maxBondOrder);
  // cout << neighbourIds.size() << endl;

  unordered_map<size_t, Eigen::RowVector3d> centeredPositions;
  centeredPositions.reserve(neighbourIds.size());

  // Reference position of the central site
  Eigen::RowVector3d centerPos = config.GetRelativePositionOfLattice(latticeId).transpose();

  for (size_t id : neighbourIds)
  {
    Eigen::RowVector3d relPos = config.GetRelativePositionOfLattice(id).transpose();

    // Translate so central site is at origin
    relPos -= centerPos;

    // Wrap into [-0.5, 0.5) for minimal image (symmetric around 0)
    relPos = relPos.unaryExpr([](double x)
                              {
                                      double wrapped = x - std::floor(x + 0.5);  // Directly computes [-0.5, 0.5)
                                      return wrapped; });

    // Shift to place central site at (0.5, 0.5, 0.5) and neighbors in [0,1)
    relPos += Eigen::RowVector3d(0.5, 0.5, 0.5);

    centeredPositions.emplace(id, relPos);
  }

  // Optionally add the central site itself at (0.5, 0.5, 0.5) if needed
  // centeredPositions.emplace(latticeId, Eigen::RowVector3d(0.5, 0.5, 0.5));

  return centeredPositions;
}

vector<size_t> GetCanonicalSortedSitesForSite(
    const Config &config,
    const size_t latticeId,
    const size_t &maxBondOrder)
{
  auto centeredLatticeIdHashmap = GetCenteredNeighboursSite(config,
                                                            latticeId,
                                                            maxBondOrder);

  // Sort the transformed positions
  // Convert unordered_map to vector for sorting
  vector<pair<size_t, Eigen::RowVector3d>> latticeIdVector(
      centeredLatticeIdHashmap.begin(),
      centeredLatticeIdHashmap.end());

  // Sort the lattice vector based on PositionCompare
  sort(latticeIdVector.begin(), latticeIdVector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  vector<size_t> sortedLatticeIdVector;
  sortedLatticeIdVector.reserve(latticeIdVector.size());
  for (const auto &pair : latticeIdVector)
  {
    sortedLatticeIdVector.push_back(pair.first);
  }

  return sortedLatticeIdVector;
}

/*
// For KRA
vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEquivalentClustersEncoding(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxBondOrderOfCluster,
    const size_t &maxClusterSize,
    const unordered_map<size_t, Eigen::RowVector3d> &canonicalReferenceMap,
    const bool debug)
{
  pair<size_t, size_t> jumpPair = {0, config.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]};

  auto nnLatticeSites = config.GetNeighboringLatticeIdSetOfPair(jumpPair, maxBondOrder);

  auto symmetryOperation = GetSymmetryOperations(config, nnLatticeSites);

  auto canonicalSortedLatticeSites = GetCanonicalSortedSitesForPair(
      config,
      jumpPair,
      maxBondOrder,
      canonicalReferenceMap,
      symmetryOperation);

  auto latticeClusterHashset = FindClustersWithinAllowedSites(
      config,
      maxClusterSize,
      maxBondOrderOfCluster,
      canonicalSortedLatticeSites);

  auto equivalentClusters = GetEquivalentClusters(
      config,
      nnLatticeSites,
      latticeClusterHashset,
      debug);

  map<size_t, size_t> latticeIdToIndexMap;

  for (size_t i = 0; i < canonicalSortedLatticeSites.size(); i++)
  {
    latticeIdToIndexMap[canonicalSortedLatticeSites[i]] = i;
  }

  // Encode the clusters
  vector<pair<vector<vector<size_t>>, LatticeClusterType>> encodedEquivalentClusters;
  encodedEquivalentClusters.reserve(equivalentClusters.size());

  for (const auto &orbits : equivalentClusters)
  {
    vector<vector<size_t>> encodedOrbits;
    LatticeClusterType clusterType;
    for (const auto &latticeCluster : orbits)
    {
      clusterType = IdentifyLatticeClusterType(config, latticeCluster);
      if (latticeCluster.size() == 0)
      {
        encodedOrbits.emplace_back(latticeCluster);
      }
      else
      {
        vector<size_t> encodedCluster;
        encodedCluster.reserve(latticeCluster.size());
        for (const auto latticeId : latticeCluster)
        {
          encodedCluster.emplace_back(latticeIdToIndexMap.at(latticeId));
        }
        encodedOrbits.emplace_back(encodedCluster);
      }
    }
    encodedEquivalentClusters.emplace_back(make_pair(encodedOrbits, clusterType));
  }

  sort(encodedEquivalentClusters.begin(), encodedEquivalentClusters.end(),
       [](const auto &a, const auto &b)
       {
         return a.second < b.second;
       });

  if (debug)
  {
    cout << "Encoded Equivalent Clusters : " << endl;
    int orbitIndex = 0;
    for (const auto &pair : encodedEquivalentClusters)
    {
      const auto &encodedOrbits = pair.first;
      const auto &clusterType = pair.second;

      cout << "Cluster Type: " << clusterType << "\n";
      cout << "Orbits:" << orbitIndex << endl;

      for (const auto &orbit : encodedOrbits)
      {
        cout << "  [ ";
        for (const auto &id : orbit)
          cout << id << " ";
        cout << "]\n";
      }

      cout << "------------------------\n";

      orbitIndex++;
    }

    GetSymmetryOperations(config, nnLatticeSites, debug, latticeIdToIndexMap);
  }

  return encodedEquivalentClusters;
}
*/
/*
// For LCE
vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEquivalentClustersEncoding(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const size_t &maxBondOrderOfCluster,
    const bool debug,
    const double symprec)
{
  size_t latticeId = 0;

  // Canonical Sorting
  auto nnLatticeSites = GetCanonicalSortedSitesForSite(config, latticeId, maxBondOrder);
  // nnLatticeSites.emplace_back(latticeId);

  unordered_set<size_t> nnLatticeSitesSet(nnLatticeSites.begin(), nnLatticeSites.end());

  auto symmetryOperation = GetSymmetryOperations(config, nnLatticeSitesSet, debug, {}, symprec);

  auto latticeClusterHashset = FindClustersWithinAllowedSites(
      config,
      maxClusterSize,
      maxBondOrderOfCluster,
      nnLatticeSites);

  auto equivalentClusters = GetEquivalentClusters(
      config,
      nnLatticeSitesSet,
      latticeClusterHashset, symprec, debug);

  map<size_t, size_t> latticeIdToIndexMap;

  for (size_t i = 0; i < nnLatticeSites.size(); i++)
  {
    latticeIdToIndexMap[nnLatticeSites[i]] = i;
  }

  // Encode the clusters
  vector<pair<vector<vector<size_t>>, LatticeClusterType>> encodedEquivalentClusters;
  encodedEquivalentClusters.reserve(equivalentClusters.size());

  for (const auto &orbits : equivalentClusters)
  {
    vector<vector<size_t>> encodedOrbits;
    LatticeClusterType clusterType;
    for (const auto &latticeCluster : orbits)
    {
      clusterType = IdentifyLatticeClusterType(config, latticeCluster);
      if (latticeCluster.size() == 0)
      {
        encodedOrbits.emplace_back(latticeCluster);
      }
      else
      {
        vector<size_t> encodedCluster;
        encodedCluster.reserve(latticeCluster.size());
        for (const auto latticeId : latticeCluster)
        {
          encodedCluster.emplace_back(latticeIdToIndexMap.at(latticeId));
        }
        encodedOrbits.emplace_back(encodedCluster);
      }
    }
    encodedEquivalentClusters.emplace_back(make_pair(encodedOrbits, clusterType));
  }

  sort(encodedEquivalentClusters.begin(), encodedEquivalentClusters.end(),
       [](const auto &a, const auto &b)
       {
         return a.second < b.second;
       });

  if (debug)
  {
    cout << "Encoded Equivalent Clusters : " << endl;
    int orbitIndex = 0;
    for (const auto &pair : encodedEquivalentClusters)
    {
      const auto &encodedOrbits = pair.first;
      const auto &clusterType = pair.second;

      cout << "Cluster Type: " << clusterType << "\n";
      cout << "Orbits: " << orbitIndex << endl;

      for (const auto &orbit : encodedOrbits)
      {
        cout << "  [ ";
        for (const auto &id : orbit)
          cout << id << " ";
        cout << "]\n";
      }

      cout << "------------------------\n";
      orbitIndex++;
    }

    GetSymmetryOperations(config, nnLatticeSitesSet, debug, latticeIdToIndexMap);
  }

  return encodedEquivalentClusters;
}
*/
/*
// Not Required because this is taken care by the ICET Source code using LOLG
// For energy per site clusters
vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEquivalentClustersEncoding(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const bool debug,
    const double symprec)
{
  size_t latticeId = 0;

  // Canonical Sorting
  auto nnLatticeSites = GetCanonicalSortedSitesForSite(config, latticeId, maxBondOrder);
  nnLatticeSites.emplace_back(latticeId);

  unordered_set<size_t> nnLatticeSitesSet(nnLatticeSites.begin(), nnLatticeSites.end());

  auto symmetryOperation = GetSymmetryOperations(config, nnLatticeSitesSet, debug, {}, symprec);

  auto latticeClusterHashset = FindAllLatticeClusters(
      config,
      maxClusterSize,
      maxBondOrder,
      {latticeId});

  auto equivalentClusters = GetEquivalentClusters(
      config,
      nnLatticeSitesSet,
      latticeClusterHashset, symprec, debug);

  map<size_t, size_t> latticeIdToIndexMap;

  for (size_t i = 0; i < nnLatticeSites.size(); i++)
  {
    latticeIdToIndexMap[nnLatticeSites[i]] = i;
  }

  // Encode the clusters
  vector<pair<vector<vector<size_t>>, LatticeClusterType>> encodedEquivalentClusters;
  encodedEquivalentClusters.reserve(equivalentClusters.size());

  for (const auto &orbits : equivalentClusters)
  {
    vector<vector<size_t>> encodedOrbits;
    LatticeClusterType clusterType;
    for (const auto &latticeCluster : orbits)
    {
      clusterType = IdentifyLatticeClusterType(config, latticeCluster);
      if (latticeCluster.size() == 0)
      {
        encodedOrbits.emplace_back(latticeCluster);
      }
      else
      {
        vector<size_t> encodedCluster;
        encodedCluster.reserve(latticeCluster.size());
        for (const auto latticeId : latticeCluster)
        {
          encodedCluster.emplace_back(latticeIdToIndexMap.at(latticeId));
        }
        encodedOrbits.emplace_back(encodedCluster);
      }
    }
    encodedEquivalentClusters.emplace_back(make_pair(encodedOrbits, clusterType));
  }

  sort(encodedEquivalentClusters.begin(), encodedEquivalentClusters.end(),
       [](const auto &a, const auto &b)
       {
         return a.second < b.second;
       });

  if (debug)
  {
    cout << "Encoded Equivalent Clusters : " << endl;
    for (const auto &pair : encodedEquivalentClusters)
    {
      const auto &encodedOrbits = pair.first;
      const auto &clusterType = pair.second;

      cout << "Cluster Type: " << clusterType << "\n";
      cout << "Orbits:\n";

      for (const auto &orbit : encodedOrbits)
      {
        cout << "  [ ";
        for (const auto &id : orbit)
          cout << id << " ";
        cout << "]\n";
      }

      cout << "------------------------\n";
    }

    GetSymmetryOperations(config, nnLatticeSitesSet, debug, latticeIdToIndexMap);
  }

  return encodedEquivalentClusters;
}

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
=======
>>>>>>> Stashed changes
*/
/*
>>>>>>> Stashed changes
unordered_map<size_t, Eigen::RowVector3d>
GetCenteredNeighboursSite(const Config &config,
                          size_t latticeId,
                          size_t maxBondOrder)
{
  auto neighbourIds = config.GetNeighborLatticeIdsUpToOrder(latticeId, maxBondOrder);
  unordered_map<size_t, Eigen::RowVector3d> centeredPositions;
  centeredPositions.reserve(neighbourIds.size());

  // Reference position of the central site
  Eigen::RowVector3d centerPos = config.GetRelativePositionOfLattice(latticeId).transpose();

  for (size_t id : neighbourIds)
  {
    Eigen::RowVector3d relPos = config.GetRelativePositionOfLattice(id).transpose();

    // Translate so central site is at origin
    relPos -= centerPos;

<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
    // Wrap into [0,1) box for periodic boundary conditions
    relPos = relPos.unaryExpr([](double x)
                              {
                                      double wrapped = x - floor(x);
                                      if (wrapped >= 1.0) wrapped -= 1.0;
                                      if (wrapped < 0.0)  wrapped += 1.0;
                                      return wrapped; });
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
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

    centeredPositions.emplace(id, relPos);
  }

  return centeredPositions;
}
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream

vector<size_t> GetCanonicalSortedSitesForSite(
    const Config &config,
    const size_t latticeId,
    const size_t &maxBondOrder)
{
  auto centeredLatticeIdHashmap = GetCenteredNeighboursSite(config,
                                                            latticeId,
                                                            maxBondOrder);

  // Sort the transformed positions
  // Convert unordered_map to vector for sorting
  vector<pair<size_t, Eigen::RowVector3d>> latticeIdVector(
      centeredLatticeIdHashmap.begin(),
      centeredLatticeIdHashmap.end());

  // Sort the lattice vector based on PositionCompare
  sort(latticeIdVector.begin(), latticeIdVector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  vector<size_t> sortedLatticeIdVector;
  sortedLatticeIdVector.reserve(latticeIdVector.size());
  for (const auto &pair : latticeIdVector)
  {
    sortedLatticeIdVector.push_back(pair.first);
  }

  return sortedLatticeIdVector;
}
=======
*/
>>>>>>> Stashed changes
=======
*/
>>>>>>> Stashed changes
=======
*/
>>>>>>> Stashed changes
=======
*/
>>>>>>> Stashed changes
=======
*/
>>>>>>> Stashed changes
=======
*/
>>>>>>> Stashed changes
=======
*/
>>>>>>> Stashed changes
=======
*/
>>>>>>> Stashed changes
