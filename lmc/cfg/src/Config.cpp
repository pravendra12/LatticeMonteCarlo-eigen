/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 1/16/20 3:55 AM                                                                         *
 * @Last Modified by: pravendra12                                                                    *
 * @Last Modified time: 11/28/24 8:40 PM                                                          *
 **************************************************************************************************/

/*! \file  Config.cpp
 *  \brief File for the Config class implementation.
 */

#include "Config.h"

Config::Config() = default;

Config::Config(Eigen::Matrix3d basis,
               Eigen::Matrix3Xd relative_position_matrix,
               std::vector<Element> atom_vector)
    : basis_(std::move(basis)),
      relative_position_matrix_(std::move(relative_position_matrix)),
      atom_vector_(std::move(atom_vector))
{
  // Check that the number of atoms is the same in the basis vectors and the relative position matrix
  if (atom_vector_.size() != static_cast<size_t>(relative_position_matrix_.cols()))
  {
    throw std::runtime_error("Lattice vector and atom vector size do not match, lattice_vector.size() = " +
                             std::to_string(relative_position_matrix_.cols()) + ", atom_vector_.size() = " +
                             std::to_string(atom_vector_.size()));
  }

  cartesian_position_matrix_ = basis_ * relative_position_matrix_;
  for (size_t id = 0; id < atom_vector_.size(); ++id)
  {
    // id here is also lattice_id
    lattice_to_atom_hashmap_.emplace(id, id);
    atom_to_lattice_hashmap_.emplace(id, id);
  }
}

Config::~Config() = default;

size_t Config::GetNumAtoms() const
{
  return atom_vector_.size();
}

const std::vector<Element> &Config::GetAtomVector() const
{
  return atom_vector_;
}

size_t Config::GetNumLattices() const
{
  return static_cast<size_t>(relative_position_matrix_.cols());
}

const Eigen::Matrix3d &Config::GetBasis() const
{
  return basis_;
}

const std::vector<std::vector<std::vector<size_t>>> &Config::GetNeighborLists() const
{
  return neighbor_lists_;
}

size_t Config::GetVacancyLatticeId() const
{

  for (size_t i = 0; i < atom_vector_.size(); i++)
  {
    if (atom_vector_[i] == ElementName::X)
    {
      return atom_to_lattice_hashmap_.at(i);
    }
  }

  throw std::runtime_error("Vacancy Not Found");

  return 0;
}

const Eigen::Matrix3Xd &Config::GetRelativePositionMatrix() const
{
  return relative_position_matrix_;
}

size_t Config::GetCentralAtomLatticeId() const
{

  Eigen::Vector3d centralAtomRelativePosition{0.5, 0.5, 0.5};

  for (size_t i = 0; i < relative_position_matrix_.cols(); ++i)
  {
    if (relative_position_matrix_.col(i).isApprox(centralAtomRelativePosition, 1e-6))
    {
      return i;
    }
  }

  // Exit the program if no match is found
  std::cerr << "Error: Central atom not found in the lattice!" << std::endl;
  exit(EXIT_FAILURE); // Exit with failure status
}

size_t Config::GetVacancyAtomId() const
{

  for (size_t i = 0; i < atom_vector_.size(); i++)
  {
    if (atom_vector_[i] == ElementName::X)
    {
      return i;
    }
  }

  throw std::runtime_error("Vacancy Not Found");

  return 0;
}

std::map<Element, double> Config::GetConcentration() const
{
  std::map<Element, double> concentration;

  // Count occurrences
  for (const auto &element : atom_vector_)
  {
    concentration[element]++;
  }
  // Total number of atoms
  double total_atoms = static_cast<double>(GetNumAtoms());

  // Calculate concentration fractions
  for (auto &pair : concentration)
  {
    pair.second /= total_atoms;
  }

  return concentration;
}

std::map<Element, std::vector<size_t>>
Config::GetElementOfAtomIdVectorMap() const
{
  std::map<Element, std::vector<size_t>> element_list_map;
  for (size_t id = 0; id < atom_vector_.size(); id++)
  {
    element_list_map[atom_vector_[id]].push_back(id);
  }
  return element_list_map;
}

std::vector<size_t> Config::GetNeighborAtomIdVectorOfAtom(size_t atom_id, size_t distance_order) const
{
  auto lattice_id = atom_to_lattice_hashmap_.at(atom_id);
  const auto &neighbor_lattice_id_vector = neighbor_lists_.at(distance_order - 1).at(lattice_id);
  std::vector<size_t> neighbor_atom_id_vector;
  neighbor_atom_id_vector.reserve(neighbor_lattice_id_vector.size());

  std::transform(neighbor_lattice_id_vector.begin(), neighbor_lattice_id_vector.end(),
                 std::back_inserter(neighbor_atom_id_vector),
                 [this](const auto &neighbor_lattice_id)
                 { return lattice_to_atom_hashmap_.at(neighbor_lattice_id); });
  return neighbor_atom_id_vector;
}

std::vector<size_t> Config::GetNeighborLatticeIdVectorOfLattice(size_t lattice_id, size_t distance_order) const
{
  // Get the vector of neighbor lattice IDs for the specified lattice ID and distance order
  const auto &neighbor_lattice_id_vector = neighbor_lists_.at(distance_order - 1).at(lattice_id);

  // Return the vector directly as it already contains the neighbor lattice IDs
  return neighbor_lattice_id_vector;
}

vector<size_t> Config::GetNeighborLatticeIdsUpToOrder(size_t latticeId, size_t maxBondOrder) const
{
  // Validate maxBondOrder
  if (maxBondOrder == 0 || maxBondOrder > neighbor_lists_.size())
  {
    std::cerr << "Error: maxBondOrder (" << maxBondOrder
              << ") must be between 1 and " << neighbor_lists_.size()
              << " (the number of available neighbor shells)." << std::endl;
    throw std::out_of_range("Invalid maxBondOrder in GetNeighborLatticeIdsUpToOrder()");
  }

  std::vector<size_t> neighbors;
  for (size_t order = 1; order < maxBondOrder + 1; ++order)
  {
    const auto &ids = GetNeighborLatticeIdVectorOfLattice(latticeId, order);
    neighbors.insert(neighbors.end(), ids.begin(), ids.end());
  }

  // cout << "Number of NN upto " << maxBondOrder << " maxBondOrder: " << neighbors.size() << endl;

  return neighbors;
}

Element Config::GetElementOfLattice(size_t lattice_id) const
{
  return atom_vector_.at(lattice_to_atom_hashmap_.at(lattice_id));
}

Element Config::GetElementOfAtom(size_t atom_id) const
{
  return atom_vector_.at(atom_id);
}

size_t Config::GetAtomIdOfLattice(size_t latticeId) const
{
  return lattice_to_atom_hashmap_.at(latticeId);
}

size_t Config::GetLatticeIdOfAtom(size_t atomId) const
{
  return atom_to_lattice_hashmap_.at(atomId);
}

Eigen::Ref<const Eigen::Vector3d> Config::GetRelativePositionOfLattice(size_t lattice_id) const
{
  return relative_position_matrix_.col(static_cast<int>(lattice_id));
}

Eigen::Ref<const Eigen::Vector3d> Config::GetRelativePositionOfAtom(size_t atom_id) const
{
  return GetRelativePositionOfLattice(atom_to_lattice_hashmap_.at(atom_id));
}

Eigen::Ref<const Eigen::Vector3d> Config::GetCartesianPositionOfLattice(size_t lattice_id) const
{
  return cartesian_position_matrix_.col(static_cast<int>(lattice_id));
}

Eigen::Ref<const Eigen::Vector3d> Config::GetCartesianPositionOfAtom(size_t atom_id) const
{
  return GetCartesianPositionOfLattice(atom_to_lattice_hashmap_.at(atom_id));
}

/// Functions which are used in Symmetry Operations.

inline bool PositionCompareState(const std::pair<size_t, Eigen::RowVector3d> &lhs,
                                 const std::pair<size_t, Eigen::RowVector3d> &rhs)
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

std::vector<size_t>
Config::GetSortedLatticeVectorStateOfPair(
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    const size_t &max_bond_order) const
{

  // Take one direction as canonical direction : 111
  // Get the position for this direction
  // Then for other direction use the symmetry operation to match the cords

  auto neighboring_lattice_ids = GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair,
                                                                  max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();
  Eigen::RowVector3d move_distance =
      Eigen::RowVector3d(0.5, 0.5, 0.5) - GetLatticePairCenter(lattice_id_jump_pair);

  std::unordered_map<size_t, Eigen::RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids)
  {
    Eigen::RowVector3d relative_position = GetRelativePositionOfLattice(id).transpose();
    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return std::floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  // Rotate lattice vectors
  RotateLatticeVector(lattice_id_hashmap, GetLatticePairRotationMatrix(lattice_id_jump_pair));

  // Convert unordered_map to vector for sorting
  std::vector<std::pair<size_t, Eigen::RowVector3d>>
      lattice_id_vector(lattice_id_hashmap.begin(), lattice_id_hashmap.end());

  // Sort the lattice vector based on PositionCompare
  std::sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  std::vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto &pair : lattice_id_vector)
  {
    sorted_lattice_ids.push_back(pair.first);
    // std::cout << pair.first << " " << pair.second << std::endl;
  }

  return sorted_lattice_ids;
}

std::vector<size_t> Config::GetSortedLatticeVectorStateWithPair(
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    const size_t &max_bond_order) const
{

  auto neighboring_lattice_ids = GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair,
                                                                  max_bond_order);
  neighboring_lattice_ids.reserve(neighboring_lattice_ids.size() + 2);

  neighboring_lattice_ids.insert(lattice_id_jump_pair.first);
  neighboring_lattice_ids.insert(lattice_id_jump_pair.second);

  size_t num_sites = neighboring_lattice_ids.size();
  Eigen::RowVector3d move_distance =
      Eigen::RowVector3d(0.5, 0.5, 0.5) - GetLatticePairCenter(lattice_id_jump_pair);

  std::unordered_map<size_t, Eigen::RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids)
  {
    Eigen::RowVector3d relative_position = GetRelativePositionOfLattice(id).transpose();
    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return std::floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  // Rotate lattice vectors
  RotateLatticeVector(lattice_id_hashmap, GetLatticePairRotationMatrix(lattice_id_jump_pair));

  // Convert unordered_map to vector for sorting
  std::vector<std::pair<size_t, Eigen::RowVector3d>>
      lattice_id_vector(lattice_id_hashmap.begin(), lattice_id_hashmap.end());

  // Sort the lattice vector based on PositionCompare
  std::sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  std::vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto &pair : lattice_id_vector)
  {
    sorted_lattice_ids.push_back(pair.first);
    // std::cout << pair.first << " " << pair.second << std::endl;
  }

  return sorted_lattice_ids;
}

Eigen::RowVector3d Config::GetLatticePairCenter(
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const
{

  // Get relative positions of the lattice IDs
  Eigen::Vector3d first_relative = GetRelativePositionOfLattice(lattice_id_jump_pair.first);
  Eigen::Vector3d second_relative = GetRelativePositionOfLattice(lattice_id_jump_pair.second);

  Eigen::Vector3d center_position;
  for (int kDim = 0; kDim < 3; ++kDim)
  {
    double distance = first_relative[kDim] - second_relative[kDim];
    int period = static_cast<int>(distance / 0.5);

    // Adjust the positions once based on the period
    if (period != 0)
    {
      first_relative[kDim] -= period;
    }

    center_position[kDim] = 0.5 * (first_relative[kDim] + second_relative[kDim]);
  }

  return center_position.transpose();
}

void Config::RotateLatticeVector(
    std::unordered_map<size_t, Eigen::RowVector3d> &lattice_id_hashmap,
    const Eigen::Matrix3d &rotation_matrix)
{

  const Eigen::RowVector3d move_distance_after_rotation =
      Eigen::RowVector3d(0.5, 0.5, 0.5) - (Eigen::RowVector3d(0.5, 0.5, 0.5) * rotation_matrix);

  for (auto &lattice : lattice_id_hashmap)
  {

    auto relative_position = lattice.second;
    // rotate
    relative_position = relative_position * rotation_matrix;
    // move to new center
    relative_position += move_distance_after_rotation;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return std::floor(x); });

    lattice.second = relative_position;
  }
}

Eigen::Matrix3d Config::GetLatticePairRotationMatrix(
    const std::pair<size_t, size_t> &lattice_id_jump_pair) const
{

  // Get Pair Direction
  Eigen::Vector3d relative_distance_vector_pair = GetRelativeDistanceVectorLattice(
      lattice_id_jump_pair.first, lattice_id_jump_pair.second);
  Eigen::RowVector3d pair_direction = relative_distance_vector_pair.transpose() * basis_;
  pair_direction.normalize();

  Eigen::RowVector3d vertical_vector;

  // First nearest neighbors
  std::vector<size_t> nn_list = GetNeighborLatticeIdVectorOfLattice(lattice_id_jump_pair.first, 1);

  std::sort(nn_list.begin(), nn_list.end());

  for (const auto nn_id : nn_list)
  {
    // Get Jump Vector
    Eigen::Vector3d relative_distance_vector = GetRelativeDistanceVectorLattice(
        lattice_id_jump_pair.first, nn_id);
    Eigen::RowVector3d jump_vector = relative_distance_vector.transpose() * basis_;
    jump_vector.normalize();

    // Check for vertical direction
    const double dot_prod = pair_direction.dot(jump_vector);
    if (std::abs(dot_prod) < constants::kEpsilon)
    {
      vertical_vector = jump_vector;
      break;
    }
  }

  // Rotation Matrix
  Eigen::Matrix3d rotation_matrix;
  rotation_matrix.row(0) = pair_direction;
  rotation_matrix.row(1) = vertical_vector;
  rotation_matrix.row(2) = pair_direction.cross(vertical_vector);

  return rotation_matrix.transpose();
}

std::unordered_set<size_t> Config::GetNeighboringLatticeIdSetOfPair(
    const std::pair<size_t, size_t> &lattice_id_pair,
    const size_t &max_bond_order) const
{

  // unordered set to store the neighboring_lattice_ids
  std::unordered_set<size_t> neighboring_lattice_ids;

  // Combine neighbors for both lattice_ids in a single loop
  for (size_t bond_order = 1; bond_order <= max_bond_order; ++bond_order)
  {
    auto neighbors_vector_1 = GetNeighborLatticeIdVectorOfLattice(lattice_id_pair.first, bond_order);
    auto neighbors_vector_2 = GetNeighborLatticeIdVectorOfLattice(lattice_id_pair.second, bond_order);

    // Insert neighbors from both lattice_id1 and lattice_id2
    for (const auto &neighbor_id : neighbors_vector_1)
    {
      if (neighbor_id != lattice_id_pair.first && neighbor_id != lattice_id_pair.second)
      {
        neighboring_lattice_ids.insert(neighbor_id);
      }
    }

    for (const auto &neighbor_id : neighbors_vector_2)
    {
      if (neighbor_id != lattice_id_pair.first && neighbor_id != lattice_id_pair.second)
      {
        neighboring_lattice_ids.insert(neighbor_id);
      }
    }
  }

  return neighboring_lattice_ids;
}

void Config::SetPeriodicBoundaryCondition(const std::array<bool, 3> &periodic_boundary_condition)
{
  periodic_boundary_condition_ = periodic_boundary_condition;
}

void Config::Wrap()
{
  // periodic boundary conditions
  for (int col = 0; col < relative_position_matrix_.cols(); ++col)
  {
    Eigen::Ref<Eigen::Vector3d> related_position = relative_position_matrix_.col(col);
    // for (Eigen::Ref<Eigen::Vector3d> related_position : relative_position_matrix_.colwise()) {
    for (const int kDim : std::vector<int>{0, 1, 2})
    {
      if (periodic_boundary_condition_[static_cast<size_t>(kDim)])
      {
        while (related_position[kDim] >= 1)
        {
          related_position[kDim] -= 1;
        }
        while (related_position[kDim] < 0)
        {
          related_position[kDim] += 1;
        }
      }
    }
  }
  cartesian_position_matrix_ = basis_ * relative_position_matrix_;
}

void Config::SetElementOfAtom(size_t atom_id, Element element_type)
{
  atom_vector_.at(atom_id) = Element(element_type);

  // Get the lattice Id corresponding to a given atom id
  auto latticeId = atom_to_lattice_hashmap_.at(atom_id);
}

void Config::SetElementOfLattice(size_t lattice_id, Element element_type)
{
  auto atom_id = lattice_to_atom_hashmap_.at(lattice_id);
  atom_vector_.at(atom_id) = Element(element_type);
}

void Config::LatticeJump(const std::pair<size_t, size_t> &lattice_id_jump_pair)
{
  const auto [lattice_id_lhs, lattice_id_rhs] = lattice_id_jump_pair;

  // std::cout << "Before Swap : " << std::endl;
  // std::cout << "Lattice ID : " << "{ " << lattice_id_lhs << GetElementOfLattice(lattice_id_lhs) <<", " <<
  //                                       lattice_id_rhs << GetElementOfLattice(lattice_id_rhs) << "}" << std::endl;

  const auto atom_id_lhs = lattice_to_atom_hashmap_.at(lattice_id_lhs);
  const auto atom_id_rhs = lattice_to_atom_hashmap_.at(lattice_id_rhs);

  // std::cout << "Atom ID : " << "{ " << atom_id_lhs << GetElementOfAtom(atom_id_lhs) <<", " <<
  //                                       atom_id_rhs << GetElementOfAtom(atom_id_rhs) << "}" << std::endl;

  atom_to_lattice_hashmap_.at(atom_id_lhs) = lattice_id_rhs;
  atom_to_lattice_hashmap_.at(atom_id_rhs) = lattice_id_lhs;
  lattice_to_atom_hashmap_.at(lattice_id_lhs) = atom_id_rhs;
  lattice_to_atom_hashmap_.at(lattice_id_rhs) = atom_id_lhs;

  //  std::cout << "After Swap : " << std::endl;
  //  std::cout << "Lattice ID : " << "{ " << lattice_id_lhs << GetElementOfLattice(lattice_id_lhs) << ", " <<
  //                                        lattice_id_rhs << GetElementOfLattice(lattice_id_rhs) << "}" << std::endl;
  //
  //
  //   std::cout << "Atom ID : " << "{ " << atom_id_lhs << GetElementOfAtom(atom_id_lhs) <<", " <<
  //                                        atom_id_rhs << GetElementOfAtom(atom_id_rhs) << "}" << std::endl;
}

void Config::ReassignLattice()
{
  std::vector<std::pair<Eigen::Vector3d, size_t>> new_lattice_id_vector(GetNumLattices());
  for (size_t i = 0; i < GetNumLattices(); ++i)
  {
    new_lattice_id_vector[i] = std::make_pair(relative_position_matrix_.col(static_cast<int>(i)), i);
  }
  std::sort(new_lattice_id_vector.begin(), new_lattice_id_vector.end(),
            [](const auto &lhs, const auto &rhs) -> bool
            {
              const double x_diff = lhs.first[0] - rhs.first[0];
              if (x_diff < -constants::kEpsilon)
              {
                return true;
              }
              if (x_diff > constants::kEpsilon)
              {
                return false;
              }
              const double y_diff = lhs.first[1] - rhs.first[1];
              if (y_diff < -constants::kEpsilon)
              {
                return true;
              }
              if (y_diff > constants::kEpsilon)
              {
                return false;
              }
              return lhs.first[2] < rhs.first[2] - constants::kEpsilon;
            });

  std::unordered_map<size_t, size_t> old_lattice_id_to_new;
  for (size_t lattice_id = 0; lattice_id < GetNumAtoms(); ++lattice_id)
  {
    old_lattice_id_to_new.emplace(new_lattice_id_vector.at(lattice_id).second, lattice_id);
  }
  std::unordered_map<size_t, size_t> new_lattice_to_atom_hashmap, new_atom_to_lattice_hashmap;
  for (size_t atom_id = 0; atom_id < GetNumAtoms(); ++atom_id)
  {
    auto old_lattice_id = atom_to_lattice_hashmap_.at(atom_id);
    auto new_lattice_id = old_lattice_id_to_new.at(old_lattice_id);
    new_lattice_to_atom_hashmap.emplace(new_lattice_id, atom_id);
    new_atom_to_lattice_hashmap.emplace(atom_id, new_lattice_id);
  }
  // Copy sorted columns back into matrix
  for (size_t i = 0; i < GetNumLattices(); ++i)
  {
    relative_position_matrix_.col(static_cast<int>(i)) = new_lattice_id_vector[i].first;
  }
  cartesian_position_matrix_ = basis_ * relative_position_matrix_;
  lattice_to_atom_hashmap_ = new_lattice_to_atom_hashmap;
  atom_to_lattice_hashmap_ = new_atom_to_lattice_hashmap;
}

Eigen::Vector3d Config::GetNormalizedDirection(size_t referenceId, size_t latticeId) const
{
  auto diff = GetRelativePositionOfLattice(referenceId) -
              GetRelativePositionOfLattice(latticeId);

  double normFactor = diff.cwiseAbs().maxCoeff();

  if (normFactor == 0)
  {
    return Eigen::Vector3d::Zero();
  }

  return diff / normFactor;
}

Eigen::Vector3d Config::GetRelativeDistanceVectorLattice(size_t lattice_id1, size_t lattice_id2) const
{
  Eigen::Vector3d relative_distance_vector = relative_position_matrix_.col(static_cast<int>(lattice_id2)) - relative_position_matrix_.col(static_cast<int>(lattice_id1));
  // periodic boundary conditions
  for (const int kDim : std::vector<int>{0, 1, 2})
  {
    if (periodic_boundary_condition_[static_cast<size_t>(kDim)])
    {
      while (relative_distance_vector[kDim] >= 0.5)
      {
        relative_distance_vector[kDim] -= 1;
      }
      while (relative_distance_vector[kDim] < -0.5)
      {
        relative_distance_vector[kDim] += 1;
      }
    }
  }
  return relative_distance_vector;
}

size_t Config::GetDistanceOrder(size_t lattice_id1, size_t lattice_id2) const
{
  if (lattice_id1 == lattice_id2)
  {
    return 0;
  } // same lattice
  Eigen::Vector3d relative_distance_vector = GetRelativeDistanceVectorLattice(lattice_id1, lattice_id2);
  double cartesian_distance = (basis_ * relative_distance_vector).norm();
  auto upper = std::upper_bound(cutoffs_.begin(), cutoffs_.end(), cartesian_distance);

  // std::cout << std::endl;

  if (upper == cutoffs_.end())
  {
    return std::numeric_limits<size_t>::max(); // distance is larger than the largest cutoff
  }

  // std::cout << "Distance_Order: " << static_cast<size_t>(1 + std::distance(cutoffs_.begin(), upper)) << std::endl;

  return static_cast<size_t>(1 + std::distance(cutoffs_.begin(), upper));
}

void Config::UpdateNeighborList(std::vector<double> cutoffs)
{
  cutoffs_ = std::move(cutoffs);
  std::sort(cutoffs_.begin(), cutoffs_.end());
  std::vector<double> cutoffs_squared(cutoffs_.size());
  std::transform(cutoffs_.begin(),
                 cutoffs_.end(),
                 cutoffs_squared.begin(),
                 [](double cutoff)
                 { return cutoff * cutoff; });

  std::vector<std::vector<size_t>> neighbors_list(GetNumLattices());
  neighbor_lists_ = {cutoffs_.size(), neighbors_list};

  // // Check if the max cutoff is greater than the radius of the max sphere that can fit in the cell
  // const auto basis_inverse_tran = basis_.inverse().transpose();
  // double cutoff_pbc = std::numeric_limits<double>::infinity();
  // for (const int kDim : std::vector<int>{0, 1, 2}) {
  //   if (periodic_boundary_condition_[static_cast<size_t>(kDim)]) {
  //     cutoff_pbc = std::min(cutoff_pbc, 0.5 / basis_inverse_tran.col(kDim).norm());
  //   }
  // }
  // if (cutoff_pbc <= cutoffs_.back()) {
  //   throw std::runtime_error(
  //       "The cutoff is larger than the maximum cutoff allowed due to periodic boundary conditions, "
  //       "cutoff_pbc = " + std::to_string(cutoff_pbc) + ", cutoff_input = " + std::to_string(cutoffs_.back()));
  // }

  // Calculate the number of cells in each dimension, at least 3 and then create cells
  num_cells_ = ((basis_.colwise().norm().array()) / (0. + cutoffs_.back())).floor().cast<int>();
  num_cells_ = num_cells_.cwiseMax(Eigen::Vector3i::Constant(3));

  cells_ = std::vector<std::vector<size_t>>(static_cast<size_t>(num_cells_.prod()));
  for (size_t lattice_id = 0; lattice_id < GetNumLattices(); ++lattice_id)
  {
    const Eigen::Vector3d relative_position = relative_position_matrix_.col(static_cast<int>(lattice_id));
    Eigen::Vector3i cell_pos = (num_cells_.cast<double>().array() * relative_position.array()).floor().cast<int>();
    int cell_idx = (cell_pos(0) * num_cells_(1) + cell_pos(1)) * num_cells_(2) + cell_pos(2);
    cells_.at(static_cast<size_t>(cell_idx)).push_back(lattice_id);
  }

  // // Check if the max distance between two atoms in the neighboring cells is greater than the cutoff
  // Eigen::Matrix3d cell_basis = basis_.array().colwise() / num_cells_.cast<double>().array();
  // double cutoff_cell = 1 / (cell_basis.inverse().colwise().norm().maxCoeff());
  // if (cutoff_cell <= cutoffs_.back()) {
  //   throw std::runtime_error(
  //       "The cutoff is larger than the maximum cutoff allowed due to non-orthogonal cells "
  //       "cutoff_cell = " + std::to_string(cutoff_cell) + ", cutoff_input = " + std::to_string(cutoffs_.back()));
  // }

  static const std::vector<std::tuple<int, int, int>> OFFSET_LIST = []()
  {
    std::vector<std::tuple<int, int, int>> offsets;
    for (int x : {-1, 0, 1})
    {
      for (int y : {-1, 0, 1})
      {
        for (int z : {-1, 0, 1})
        {
          offsets.emplace_back(x, y, z);
        }
      }
    }
    return offsets;
  }();

  // Create neighbor list, iterate over each cell and find neighboring points
  for (int cell_idx = 0; cell_idx < num_cells_.prod(); ++cell_idx)
  {
    auto &cell = cells_.at(static_cast<size_t>(cell_idx));
    int i = cell_idx / (num_cells_[1] * num_cells_[2]);
    int j = (cell_idx % (num_cells_[1] * num_cells_[2])) / num_cells_[2];
    int k = cell_idx % num_cells_[2];
    // Check neighboring cells, taking into account periodic boundaries
    for (auto [di, dj, dk] : OFFSET_LIST)
    {
      int ni = (i + di + num_cells_(0)) % num_cells_(0);
      int nj = (j + dj + num_cells_(1)) % num_cells_(1);
      int nk = (k + dk + num_cells_(2)) % num_cells_(2);
      int neighbor_cell_idx = (ni * num_cells_(1) + nj) * num_cells_(2) + nk;
      auto &neighbor_cell = cells_.at(static_cast<size_t>(neighbor_cell_idx));
      // For each point in the cell, check if it's close to any point in the neighboring cell
      for (size_t lattice_id1 : cell)
      {
        for (size_t lattice_id2 : neighbor_cell)
        {
          // Make sure we're not comparing a point to itself, and don't double-count pairs within the same cell
          if (lattice_id2 >= lattice_id1)
          {
            continue;
          }
          // Calculate distance
          const double cartesian_distance_squared =
              (basis_ * GetRelativeDistanceVectorLattice(lattice_id1, lattice_id2)).squaredNorm();
          // If the distance is less than the cutoff, the points are bonded
          for (size_t cutoff_squared_id = 0; cutoff_squared_id < cutoffs_squared.size(); ++cutoff_squared_id)
          {
            if (cartesian_distance_squared < cutoffs_squared.at(cutoff_squared_id))
            {
              neighbor_lists_.at(cutoff_squared_id).at(lattice_id1).push_back(lattice_id2);
              neighbor_lists_.at(cutoff_squared_id).at(lattice_id2).push_back(lattice_id1);
              break;
            }
          }
        }
      }
    }
  }
}

Config Config::ReadMap(const std::string &lattice_filename,
                       const std::string &element_filename,
                       const std::string &map_filename)
{
  std::ifstream ifs_lattice(lattice_filename, std::ifstream::in);
  if (!ifs_lattice)
  {
    throw std::runtime_error("Cannot open " + lattice_filename);
  }
  size_t num_atoms;
  ifs_lattice >> num_atoms;
  ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  Eigen::Matrix3d basis;
  ifs_lattice >> basis(0, 0) >> basis(0, 1) >> basis(0, 2);              // lattice vector a
  ifs_lattice >> basis(1, 0) >> basis(1, 1) >> basis(1, 2);              // lattice vector b
  ifs_lattice >> basis(2, 0) >> basis(2, 1) >> basis(2, 2);              // lattice vector c
  ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // finish this line

  Eigen::Matrix3Xd relative_position_matrix(3, num_atoms);
  double position_X, position_Y, position_Z;
  for (size_t lattice_id = 0; lattice_id < num_atoms; ++lattice_id)
  {
    ifs_lattice >> position_X >> position_Y >> position_Z;
    ifs_lattice.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    relative_position_matrix.col(static_cast<int>(lattice_id)) = Eigen::Vector3d{position_X, position_Y, position_Z};
  }

  std::ifstream ifs_element(element_filename, std::ifstream::in);
  if (!ifs_element)
  {
    throw std::runtime_error("Cannot open " + element_filename);
  }

  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);
  std::string type;
  for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id)
  {
    ifs_element >> type;
    atom_vector.emplace_back(type);
    ifs_element.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  Config config(basis, relative_position_matrix, atom_vector);
  std::ifstream ifs_map(map_filename, std::ifstream::in);
  if (!ifs_map)
  {
    throw std::runtime_error("Cannot open " + map_filename);
  }
  size_t lattice_id;
  for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id)
  {
    ifs_map >> lattice_id;
    config.lattice_to_atom_hashmap_.at(lattice_id) = atom_id;
    config.atom_to_lattice_hashmap_.at(atom_id) = lattice_id;
    ifs_map.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  return config;
}

Config Config::ReadCfg(const std::string &filename)
{
  std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
  if (!ifs)
  {
    throw std::runtime_error("Could not open file: " + filename);
  }
  boost::iostreams::filtering_istream fis;
  if (boost::filesystem::path(filename).extension() == ".gz")
  {
    fis.push(boost::iostreams::gzip_decompressor());
  }
  else if (boost::filesystem::path(filename).extension() == ".bz2")
  {
    fis.push(boost::iostreams::bzip2_decompressor());
  }
  fis.push(ifs);

  // "Number of particles = %i"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  size_t num_atoms;
  fis >> num_atoms;
  // A = 1.0 Angstrom (basic length-scale)
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  double basis_xx, basis_xy, basis_xz,
      basis_yx, basis_yy, basis_yz,
      basis_zx, basis_zy, basis_zz;
  // "H0(1,1) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_xx;
  // "H0(1,2) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_xy;
  // "H0(1,3) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_xz;
  // "H0(2,1) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_yx;
  // "H0(2,2) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_yy;
  // "H0(2,3) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_yz;
  // "H0(3,1) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_zx;
  // "H0(3,2) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_zy;
  // "H0(3,3) = %lf A"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '=');
  fis >> basis_zz;
  // finish this line
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // .NO_VELOCITY.
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // "entry_count = 3"
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  auto basis = Eigen::Matrix3d{{basis_xx, basis_xy, basis_xz},
                               {basis_yx, basis_yy, basis_yz},
                               {basis_zx, basis_zy, basis_zz}};
  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);
  Eigen::Matrix3Xd relative_position_matrix(3, num_atoms);

  double mass;
  std::string type;
  Eigen::Vector3d relative_position;

  std::vector<std::vector<size_t>> first_neighbors_adjacency_list,
      second_neighbors_adjacency_list, third_neighbors_adjacency_list;

  for (size_t id = 0; id < num_atoms; ++id)
  {
    fis >> mass >> type >> relative_position(0) >> relative_position(1) >> relative_position(2);
    atom_vector.emplace_back(type);
    relative_position_matrix.col(static_cast<int>(id)) = relative_position;
  }
  Config config_in = Config{basis, relative_position_matrix, atom_vector};
  config_in.ReassignLattice();
  config_in.Wrap();
  return config_in;
}

Config Config::ReadPoscar(const std::string &filename)
{
  std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
  if (!ifs)
  {
    throw std::runtime_error("Could not open file: " + filename);
  }

  boost::iostreams::filtering_istream fis;
  if (boost::filesystem::path(filename).extension() == ".gz")
  {
    fis.push(boost::iostreams::gzip_decompressor());
  }
  else if (boost::filesystem::path(filename).extension() == ".bz2")
  {
    fis.push(boost::iostreams::bzip2_decompressor());
  }
  fis.push(ifs);

  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip comment line
  double scale;
  fis >> scale; // scale factor, typically 1.0
  Eigen::Matrix3d basis;
  fis >> basis(0, 0) >> basis(0, 1) >> basis(0, 2);              // lattice vector a
  fis >> basis(1, 0) >> basis(1, 1) >> basis(1, 2);              // lattice vector b
  fis >> basis(2, 0) >> basis(2, 1) >> basis(2, 2);              // lattice vector c
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line
  basis *= scale;
  auto inverse_basis = basis.inverse();

  // Read the elements and number of atoms
  std::string buffer;
  getline(fis, buffer); // Skip line
  std::istringstream element_iss(buffer);
  getline(fis, buffer); // Element counts
  std::istringstream count_iss(buffer);

  std::string element;
  size_t num_elems;
  size_t num_atoms = 0;
  std::vector<std::pair<std::string, size_t>> elements_counts;

  while (element_iss >> element && count_iss >> num_elems)
  {
    elements_counts.emplace_back(element, num_elems);
    num_atoms += num_elems;
  }

  getline(fis, buffer); // Check if positions are relative or Cartesian
  bool relative_option =
      buffer[0] != 'C' && buffer[0] != 'c' && buffer[0] != 'K' && buffer[0] != 'k';

  // Store atom information
  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);
  Eigen::Matrix3Xd relative_position_matrix(3, num_atoms);

  size_t id_count = 0;
  double position_X, position_Y, position_Z;

  for (const auto &[element_symbol, count] : elements_counts)
  {
    for (size_t j = 0; j < count; ++j)
    {
      fis >> position_X >> position_Y >> position_Z;
      fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      if (relative_option)
      {
        // Read fractional coordinates
        relative_position_matrix.col(static_cast<int>(id_count)) =
            Eigen::Vector3d(position_X, position_Y, position_Z);
      }
      else
      {
        // Convert Cartesian coordinates to fractional coordinates
        relative_position_matrix.col(static_cast<int>(id_count)) =
            inverse_basis * Eigen::Vector3d(position_X, position_Y, position_Z);
      }
      atom_vector.emplace_back(element_symbol);
      ++id_count;
    }
  }

  Config config_in = Config{basis, relative_position_matrix, atom_vector};
  config_in.ReassignLattice();
  config_in.Wrap();

  return config_in;
}

Config Config::ReadConfig(const std::string &filename)
{
  // Get the file extension
  std::string extension = boost::filesystem::path(filename).extension().string();
  // Determine the file type and call the appropriate function
  if (extension == ".cfg")
  {
    return Config::ReadCfg(filename);
  }
  else if (extension == ".POSCAR")
  {
    return Config::ReadPoscar(filename);
  }
  else if (extension == ".gz" || extension == ".bz2")
  {
    // Handle compressed files by checking their base name
    std::string base_extension = boost::filesystem::path(boost::filesystem::path(filename).stem()).extension().string();
    if (base_extension == ".cfg")
    {
      return Config::ReadCfg(filename);
    }
    else if (base_extension == ".POSCAR")
    {
      return Config::ReadPoscar(filename);
    }
  }

  throw std::runtime_error(
      "Unsupported file format: " + filename +
      ". Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2");
}

Config Config::ReadXyz(const std::string &filename)
{
  // Open file stream (binary mode only matters for compressed reading)
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs.is_open())
  {
    throw std::runtime_error("Could not open file: " + filename);
  }

  // Setup boost filtering stream for optional decompression
  boost::iostreams::filtering_istream fis;
  auto ext = boost::filesystem::path(filename).extension().string();
  if (ext == ".gz")
  {
    fis.push(boost::iostreams::gzip_decompressor());
  }
  else if (ext == ".bz2")
  {
    fis.push(boost::iostreams::bzip2_decompressor());
  }
  fis.push(ifs); // push the underlying file stream last

  // --- Read number of atoms ---
  size_t num_atoms;
  fis >> num_atoms;
  fis.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  std::string lattice_line;
  std::getline(fis, lattice_line);

  // Extract lattice numbers
  size_t lat_start = lattice_line.find("Lattice=\"");
  Eigen::Matrix3d basis = Eigen::Matrix3d::Identity();
  if (lat_start != std::string::npos)
  {
    lat_start += 9; // move past 'Lattice="'
    size_t lat_end = lattice_line.find('"', lat_start);
    std::string lat_str = lattice_line.substr(lat_start, lat_end - lat_start);

    std::istringstream lat_iss(lat_str);
    std::vector<double> values;
    double v;
    while (lat_iss >> v)
      values.push_back(v);

    if (values.size() == 9)
    {
      basis << values[0], values[1], values[2],
          values[3], values[4], values[5],
          values[6], values[7], values[8];
    }
    else
    {
      throw std::runtime_error("Invalid lattice line, expected 9 numbers");
    }
  }

  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);

  Eigen::Matrix3Xd cartesian_position_matrix(3, num_atoms);

  // --- Read atoms line by line ---
  std::string line;
  for (size_t i = 0; i < num_atoms; ++i)
  {
    if (!std::getline(fis, line))
    {
      throw std::runtime_error("Unexpected end of file while reading atom " + std::to_string(i + 1));
    }

    std::istringstream iss(line);
    std::string symbol;
    double x, y, z;

    if (!(iss >> symbol >> x >> y >> z))
    {
      throw std::runtime_error("Failed to parse atom line " + std::to_string(i + 1));
    }

    atom_vector.emplace_back(symbol);
    cartesian_position_matrix(0, i) = x;
    cartesian_position_matrix(1, i) = y;
    cartesian_position_matrix(2, i) = z;
  }

  Eigen::Matrix3d inv_basis = basis.inverse();

  Eigen::Matrix3Xd relative_positions_matrix(3, cartesian_position_matrix.cols());
  relative_positions_matrix = inv_basis * cartesian_position_matrix;

  Config config_in = Config{basis, relative_positions_matrix, atom_vector};

  config_in.ReassignLattice();
  config_in.Wrap();

  return config_in;
}

Config Config::GenerateSupercell(
    size_t supercell_size,
    double lattice_param,
    const std::string &element_symbol,
    const std::string &structure_type)
{
  size_t num_atoms = (structure_type == "BCC") ? supercell_size * supercell_size * supercell_size * 2
                                               : supercell_size * supercell_size * supercell_size * 4;

  std::vector<Element> atom_vector;
  atom_vector.reserve(num_atoms);
  Eigen::Matrix3Xd relative_position_matrix(3, num_atoms);

  double lattice_param_supercell = lattice_param * supercell_size;
  // Define the basis for BCC or FCC
  Eigen::Matrix3d basis = lattice_param_supercell * Eigen::Matrix3d::Identity();

  // Define lattice points based on structure type
  std::vector<Eigen::Vector3d> lattice_points;
  if (structure_type == "BCC")
  {
    lattice_points = {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0.5, 0.5, 0.5)};
  }
  else if (structure_type == "FCC")
  {
    lattice_points = {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0.5, 0.5, 0),
                      Eigen::Vector3d(0.5, 0, 0.5), Eigen::Vector3d(0, 0.5, 0.5)};
  }

  size_t id_count = 0;
  for (int x = 0; x < supercell_size; ++x)
  {
    for (int y = 0; y < supercell_size; ++y)
    {
      for (int z = 0; z < supercell_size; ++z)
      {
        for (const auto &point : lattice_points)
        {
          Eigen::Vector3d pos = (Eigen::Vector3d(x, y, z) + point) / supercell_size;
          relative_position_matrix.col(id_count) = pos;
          atom_vector.emplace_back(element_symbol);
          ++id_count;
        }
      }
    }
  }

  Config supercell = Config{basis, relative_position_matrix, atom_vector};

  supercell.ReassignLattice();
  supercell.Wrap();

  return supercell;
}

Config Config::GenerateAlloySupercell(
    size_t supercell_size,
    double lattice_param,
    std::string structure_type,
    const std::vector<std::string> &element_vector,
    const std::vector<double> &composition_vector,
    unsigned seed)
{

  string baseElement;

  for (int i = 0; i < composition_vector.size(); i++)
  {
    if (composition_vector[i] != 0)
    {
      baseElement = element_vector[i];
      break;
    }
  }

  // First element in the element vector is the base element
  auto base_supercell = GenerateSupercell(supercell_size,
                                          lattice_param,
                                          baseElement,
                                          structure_type);

  auto atom_vector = base_supercell.GetAtomVector();
  auto num_atoms = base_supercell.GetNumAtoms();

  int start_index = 0;

  for (int el = 1; el < composition_vector.size(); ++el)
  {
    int num_atoms_el = composition_vector[el] * num_atoms / 100;

    // std::cout << num_atoms_el << std::endl;

    for (int i = start_index; i < start_index + num_atoms_el; ++i)
    {
      atom_vector[i] = Element(element_vector[el]);
    }
    start_index += num_atoms_el;
  }

  std::mt19937 mersenne_twister_rng_(seed);

  std::shuffle(atom_vector.begin(), atom_vector.end(), mersenne_twister_rng_);

  Config supercell = Config{base_supercell.GetBasis(),
                            base_supercell.GetRelativePositionMatrix(),
                            atom_vector};

  supercell.ReassignLattice();
  supercell.Wrap();

  return supercell;
}

void Config::WriteConfig(const std::string &filename, const Config &config_out)
{
  WriteConfigExtended(filename, config_out, {});
}

void Config::WriteConfigExtended(
    const std::string &filename,
    const Config &config_out,
    const std::map<std::string, std::vector<double>> &auxiliary_lists)
{
  std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fos;
  if (boost::filesystem::path(filename).extension() == ".gz")
  {
    fos.push(boost::iostreams::gzip_compressor());
  }
  else if (boost::filesystem::path(filename).extension() == ".bz2")
  {
    fos.push(boost::iostreams::bzip2_compressor());
  }
  fos.push(ofs);
  fos.precision(8);
  fos << "Number of particles = " << config_out.GetNumAtoms() << '\n';
  fos << "A = 1.0 Angstrom (basic length-scale)\n";
  fos << "H0(1,1) = " << config_out.basis_(0, 0) << " A\n";
  fos << "H0(1,2) = " << config_out.basis_(0, 1) << " A\n";
  fos << "H0(1,3) = " << config_out.basis_(0, 2) << " A\n";
  fos << "H0(2,1) = " << config_out.basis_(1, 0) << " A\n";
  fos << "H0(2,2) = " << config_out.basis_(1, 1) << " A\n";
  fos << "H0(2,3) = " << config_out.basis_(1, 2) << " A\n";
  fos << "H0(3,1) = " << config_out.basis_(2, 0) << " A\n";
  fos << "H0(3,2) = " << config_out.basis_(2, 1) << " A\n";
  fos << "H0(3,3) = " << config_out.basis_(2, 2) << " A\n";
  fos << ".NO_VELOCITY.\n";
  fos << "entry_count = " << 3 + auxiliary_lists.size() << "\n";

  size_t auxiliary_index = 0;
  for (const auto &auxiliary_list : auxiliary_lists)
  {
    fos << "auxiliary[" << auxiliary_index << "] = " << auxiliary_list.first << " [reduced unit]\n";
    ++auxiliary_index;
  }
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n", "", "", "", "");
  for (size_t it = 0; it < config_out.atom_vector_.size(); ++it)
  {
    const auto &atom = config_out.atom_vector_[it];
    const auto &relative_position = config_out.relative_position_matrix_.col(static_cast<int>(config_out.atom_to_lattice_hashmap_.at(it)));
    fos << atom.GetMass() << '\n'
        << atom << '\n'
        << relative_position.transpose().format(fmt);
    for (const auto &[key, auxiliary_list] : auxiliary_lists)
    {
      fos << ' ' << auxiliary_list.at(it);
    }
    fos << std::endl;
  }
}

void Config::WriteXyzExtended(const std::string &filename,
                              const Config &config_out,
                              const std::map<std::string, VectorVariant> &auxiliary_lists,
                              const std::map<std::string, ValueVariant> &global_list)
{
  std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream fos;
  if (boost::filesystem::path(filename).extension() == ".gz")
  {
    fos.push(boost::iostreams::gzip_compressor());
  }
  else if (boost::filesystem::path(filename).extension() == ".bz2")
  {
    fos.push(boost::iostreams::bzip2_compressor());
  }
  fos.push(ofs);
  fos.precision(8);
  fos << config_out.GetNumAtoms() << '\n';
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");
  fos << "Lattice=\"" << config_out.basis_.format(fmt) << "\" ";
  fos << "pbc=\"" << config_out.periodic_boundary_condition_[0] << " "
      << config_out.periodic_boundary_condition_[1] << " "
      << config_out.periodic_boundary_condition_[2] << "\" ";
  for (const auto &[key, value] : global_list)
  {
    fos << key << "=";
    std::visit([&fos](const auto &val)
               { fos << val << ' '; }, value);
  }
  fos << "Properties=species:S:1:pos:R:3";
  for (const auto &[key, auxiliary_list] : auxiliary_lists)
  {
    fos << ":" << key << ":";
    std::any ret;
    std::visit([&ret](const auto &vec)
               { ret = vec.at(0); }, auxiliary_list);
    if (ret.type() == typeid(int) || ret.type() == typeid(size_t))
    {
      fos << "I:1";
    }
    else if (ret.type() == typeid(double))
    {
      fos << "R:1";
    }
    else if (ret.type() == typeid(std::string))
    {
      fos << "S:1";
    }
    else if (ret.type() == typeid(Eigen::Vector3d) && std::any_cast<Eigen::Vector3d>(ret).size() == 3)
    {
      fos << "R:3";
    }
    else
    {
      throw std::runtime_error("Unsupported type");
    }
  }
  fos << '\n';

  for (size_t it = 0; it < config_out.atom_vector_.size(); ++it)
  {
    const auto &atom = config_out.atom_vector_[it];
    const auto &cartesian_position =
        config_out.cartesian_position_matrix_.col(static_cast<int>(config_out.atom_to_lattice_hashmap_.at(it)));
    fos << atom << ' '
        << cartesian_position.transpose().format(fmt) << ' ';

    for (const auto &[key, auxiliary_list] : auxiliary_lists)
    {
      std::any ret;
      std::visit([&ret, it](const auto &vec)
                 { ret = vec.at(it); }, auxiliary_list);
      if (ret.type() == typeid(int))
      {
        fos << std::any_cast<int>(ret) << ' ';
      }
      else if (ret.type() == typeid(size_t))
      {
        fos << std::any_cast<size_t>(ret) << ' ';
      }
      else if (ret.type() == typeid(double))
      {
        fos << std::any_cast<double>(ret) << ' ';
      }
      else if (ret.type() == typeid(std::string))
      {
        fos << std::any_cast<std::string>(ret) << ' ';
      }
      else if (ret.type() == typeid(Eigen::Vector3d))
      {
        fos << std::any_cast<Eigen::Vector3d>(ret).format(fmt) << ' ';
      }
      else
      {
        throw std::runtime_error("Unsupported type");
      }
    }
    fos << std::endl;
  }
}

//
//
void Config::WriteLattice(const std::string &filename, size_t &max_bond_order) const
{
  boost::iostreams::stream_buffer<boost::iostreams::file_sink> buffer(filename);
  std::ostream ofs(&buffer);
  ofs.precision(16);

  ofs << GetNumAtoms() << " positions in total" << '\n';
  ofs << basis_ << std::endl;

  Eigen::IOFormat pos_format(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", " ", "", "", "", "");

  for (size_t i = 0; i < GetNumLattices(); ++i)
  {
    auto relative_position_lattice = GetRelativePositionOfLattice(i);
    ofs << relative_position_lattice.transpose().format(pos_format) << " # ";

    // Write neighbors for 1st, 2nd, .. upto max_bond_order
    for (size_t order = 1; order <= max_bond_order; ++order)
    {
      const auto &neighbors = GetNeighborLatticeIdVectorOfLattice(i, order);
      for (const auto neighbor : neighbors)
      {
        ofs << neighbor << ' ';
      }
    }
    ofs << '\n';
  }
  buffer.close(); // Ensures all data is written at once
}

static void ExpandShell(
    const Config &config,
    const std::unordered_set<size_t> &currentShell,
    std::unordered_set<size_t> &visitedSet,
    std::unordered_set<size_t> &nextShell)
{
  for (const auto latticeId : currentShell)
  {
    const auto &neighbors = config.GetNeighborLatticeIdVectorOfLattice(latticeId, 1);
    for (const auto nId : neighbors)
    {
      if (visitedSet.insert(nId).second) // Insert only if not already present
      {
        nextShell.insert(nId);
      }
    }
  }
}

Config Config::ExtractLocalSupercell(
    const std::pair<size_t, size_t> &latticeJumpPair,
    size_t supercellSize,
    double latticeParam) const
{
  if (supercellSize < 2)
  {
    throw std::invalid_argument("Supercell size should be at least 2");
  }

  int maxShell = static_cast<int>(supercellSize - 1);

  std::unordered_set<size_t> latticeIdSet{latticeJumpPair.first, latticeJumpPair.second};
  std::unordered_set<size_t> latticeIdSetNN{latticeJumpPair.second};

  std::unordered_set<size_t> currentShell{latticeJumpPair.first};
  std::unordered_set<size_t> currentShellNN{latticeJumpPair.second};

  for (int shell = 0; shell < maxShell; ++shell)
  {
    std::unordered_set<size_t> nextShell;
    std::unordered_set<size_t> nextShellNN;

    ExpandShell(*this, currentShell, latticeIdSet, nextShell);
    ExpandShell(*this, currentShellNN, latticeIdSetNN, nextShellNN);

    currentShell = std::move(nextShell);
    currentShellNN = std::move(nextShellNN);
  }

  // Merge the two lattice ID sets
  latticeIdSet.insert(latticeIdSetNN.begin(), latticeIdSetNN.end());

  // Construct the output configuration
  size_t numSites = latticeIdSet.size();

  Eigen::Matrix3Xd relativePositionMatrix(3, numSites);
  std::vector<Element> atomVector;
  atomVector.reserve(numSites);

  size_t i = 0;
  for (const auto latticeId : latticeIdSet)
  {
    Eigen::Vector3d relativePosition = GetRelativePositionOfLattice(latticeId);

    relativePositionMatrix.col(i++) = relativePosition;
    atomVector.emplace_back(GetElementOfLattice(latticeId));
  }

  // Eigen::Matrix3d newBasis = Eigen::Matrix3d::Identity() * supercellSize * latticeParam;
  Eigen::Matrix3d newBasis = GetBasis();
  return Config(newBasis, relativePositionMatrix, atomVector);
}
