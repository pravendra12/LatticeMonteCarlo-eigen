#include "Config.h"
#include "PotentialEnergyEstimator.h"
#include "VacancyMigrationPredictor.h"
#include "SymmetrySpglib.h"
#include "EnergyPredictor.h"
#include "ClusterExpansionParameters.h"
#include "BasisSet.h"
#include "AtomClusterType.hpp"
#include "CorrelationVector.h"
#include "ConfigEncoding.h"
#include <filesystem>
#include <nlohmann/json.hpp>
#include <random>
#include <vector>
#include <string>

#include "JsonUtility.h"

#include "GetEquivalentClusters.h"
namespace fs = std::filesystem;
using json = nlohmann::json;

#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <Eigen/Dense> // Assuming you use Eigen for vectors/matrices

using namespace std;
using namespace Eigen;

// Small helper for minimum image distance
double MinimumImageDistance(const Vector3d &r1, const Vector3d &r2, const Matrix3d &cell)
{
  Vector3d dr = r1 - r2;
  // Convert to fractional coordinates
  Vector3d frac = cell.inverse() * dr;
  // Apply PBC
  for (int i = 0; i < 3; ++i)
  {
    frac[i] -= round(frac[i]);
  }
  // Back to Cartesian
  Vector3d d_cart = cell * frac;
  return d_cart.norm();
}
// Compute NN cutoffs automatically based on lattice IDs
vector<double> ComputeBCCNNCutoffs(const Config &cfg, size_t maxNN = 3, double eps = 1e-2)
{
  size_t numAtoms = cfg.GetNumAtoms();
  Matrix3d cell = cfg.GetBasis();

  set<double> uniqueDistances;

  for (size_t i = 0; i < numAtoms; ++i)
  {
    Vector3d ri = cfg.GetCartesianPositionOfLattice(i);

    for (size_t j = i + 1; j < numAtoms; ++j)
    {
      Vector3d rj = cfg.GetCartesianPositionOfLattice(j);

      double d = MinimumImageDistance(ri, rj, cell);

      // Round to avoid tiny floating point differences
      d = round(d / eps) * eps;

      uniqueDistances.insert(d);
    }
  }

  // Convert to sorted vector
  vector<double> sortedDistances(uniqueDistances.begin(), uniqueDistances.end());
  sort(sortedDistances.begin(), sortedDistances.end());

  // Take first `maxNN` distances as NN cutoffs
  vector<double> nnCutoffs;
  for (size_t i = 0; i < min(maxNN, sortedDistances.size()); ++i)
  {
    nnCutoffs.push_back(sortedDistances[i] + eps); // slightly bigger to include neighbors
  }

  return nnCutoffs;
}

VectorXd GetKRAEncoding(const Config &config,
                        const pair<size_t, size_t> &latticeIdPair,
                        BasisSet &atomicBasis,
                        const size_t &maxBondOrder,
                        const size_t &maxClusterSize,
                        const size_t &maxBondOrderOfCluster,
                        const unordered_map<size_t, Eigen::RowVector3d> &referenceLatticeIdHashmap,
                        const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations)
{
  auto sortedLatticeIds = GetCanonicalSortedSitesForPair(
      config,
      latticeIdPair,
      maxBondOrder,
      referenceLatticeIdHashmap,
      symmetryOperations);

  unordered_set<size_t> sortedLatticeIdsSet(sortedLatticeIds.begin(), sortedLatticeIds.end());

  auto allLatticeClustersSet = FindClustersWithinAllowedSites(
      config,
      maxClusterSize,
      maxBondOrderOfCluster,
      sortedLatticeIds);

  auto orbitVector = GetOrbitsNew(
      config,
      sortedLatticeIdsSet,
      allLatticeClustersSet,
      symmetryOperations,
      1e-5, true);

  auto correlationVector = GetCorrelationVector(config,
                                                atomicBasis,
                                                orbitVector);

  return correlationVector;
}

VectorXd GetLCEEncoding(const Config &config,
                        const size_t centralLatticeId,
                        BasisSet &atomicBasis,
                        const size_t &maxBondOrder,
                        const size_t &maxClusterSize,
                        const size_t &maxBondOrderOfCluster,
                        const unordered_map<size_t, Eigen::RowVector3d> &referenceLatticeIdHashmap,
                        const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations)
{
  auto sortedLatticeIds = GetCanonicalSortedSitesForSite(
      config,
      centralLatticeId,
      maxBondOrder);

  cout << "Size of soretedLatticeIds : " << sortedLatticeIds.size() << endl;

  unordered_set<size_t> sortedLatticeIdsSet(sortedLatticeIds.begin(), sortedLatticeIds.end());

  auto allLatticeClustersSet = FindClustersWithinAllowedSites(
      config,
      maxClusterSize,
      maxBondOrderOfCluster,
      sortedLatticeIds);

  auto orbitVector = GetOrbitsNew(
      config,
      sortedLatticeIdsSet,
      allLatticeClustersSet,
      symmetryOperations,
      1e-5, true);

  auto correlationVector = GetCorrelationVector(config,
                                                atomicBasis,
                                                orbitVector);

  return correlationVector;
}

VectorXd GetConfigEncoding(const Config &config,
                           BasisSet &atomicBasis,
                           const size_t maxBondOrder,
                           const size_t maxClusterSize,
                           const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations)
{

  unordered_set<size_t> latticeIdsSet;

  for (size_t id = 0; id < config.GetNumLattices(); id++)
  {
    latticeIdsSet.insert(id);
  }

  auto allLatticeClustersSet = FindAllLatticeClusters(
      config,
      maxClusterSize,
      maxBondOrder, {});

  auto orbitVector = GetOrbitsNew(
      config,
      latticeIdsSet,
      allLatticeClustersSet,
      symmetryOperations,
      1e-5, true);

  auto correlationVector = GetCorrelationVector(config,
                                                atomicBasis,
                                                orbitVector);

  return correlationVector;
}

void GetEncodeVectorCE(
    const Config &config,
    const size_t maxClusterSize,
    const size_t maxBondOrder,
    const unordered_map<ClusterType, size_t, boost::hash<ClusterType>> &clusterTypeCountHashMapReference,
    set<ClusterType> &initializedClusterHashSet)
{
  auto allLatticeClustersSet = FindAllLatticeClusters(config, maxClusterSize, maxBondOrder, {});

  auto clusterTypeCountHashMap = clusterTypeCountHashMapReference;

  for (const auto &latticeCluster : allLatticeClustersSet)
  {
    auto atomClusterType = IdentifyAtomClusterType(config, latticeCluster.GetLatticeIdVector());
    cout << latticeCluster.GetClusterType() << " : " << atomClusterType << endl;
    clusterTypeCountHashMap[ClusterType(atomClusterType, latticeCluster.GetClusterType())]++;
  }

  for (auto entry : clusterTypeCountHashMap)
  {
    cout << entry.first << " : " << entry.second << endl;
  }

  VectorXd encodeVector(initializedClusterHashSet.size());
  int idx = 0;
  for (const auto &clusterType : initializedClusterHashSet)
  {
    auto clusterTypeCount = clusterTypeCountHashMap.at(clusterType);

    // Count of Cluster Types for given configuration
    auto countCluster = static_cast<double>(clusterTypeCount);
    // Count of Cluster Types for normalization
    /// auto totalLatticeCluster = static_cast<double>(latticeClusterCountHashMap_.at(clusterType.lattice_cluster_type_));

    // encodeVector(idx) = countCluster / totalLatticeCluster;

    std::cout << clusterType << " : " << countCluster << endl;
    // " : " << totalLatticeCluster << std::endl;
    ++idx;
  }
  // return encodeVector;
}

unordered_map<ClusterType, size_t, boost::hash<ClusterType>> ConvertSetToHashMapMain(
    const set<ClusterType> &clusterTypeSet)
{
  unordered_map<ClusterType, size_t, boost::hash<ClusterType>> clusterTypeCount;
  for (const auto &clusterType : clusterTypeSet)
  {
    clusterTypeCount[clusterType] = 0;
  }
  return clusterTypeCount;
}
/*
int main()
{
  size_t maxClusterSizeCE = 3;
  size_t maxBondOrderCE = 3;
  size_t maxBondOrderOfCluster = 3;

  size_t maxBondOrderKRA = 2;
  size_t maxClusterSizeKRA = 3;

  /////////////////// Reference Config ///////////////////////////////////
  auto cfgRef = Config::ReadPoscar("//media/sf_Phd/primitiveBCC.POSCAR");
  // cfgRef.UpdateNeighborList({3, 4.7, 5.6});
  vector<double> nnCutoffsRef = ComputeBCCNNCutoffs(cfgRef, 3, 0.01);
  cout << "Automatic NN cutoffs: ";
  for (auto c : nnCutoffsRef)
    cout << c << " ";

  cout << endl;

  cfgRef.UpdateNeighborList(nnCutoffsRef);

  cout << cfgRef.GetNeighborLatticeIdVectorOfLattice(0, 1).size() << endl;
  cout << cfgRef.GetNeighborLatticeIdVectorOfLattice(0, 2).size() << endl;
  cout << cfgRef.GetNeighborLatticeIdVectorOfLattice(0, 3).size() << endl;

  set<Element> elementSet = {Element("Mo"), Element("Ta")};

  auto initializedClusterHashSet = InitializeClusterTypeSet(
      cfgRef,
      elementSet,
      maxClusterSizeCE,
      maxBondOrderCE);

  auto clusterTypeCountHashMap = ConvertSetToHashMapMain(initializedClusterHashSet);

  vector<pair<ClusterType, size_t>> entries(clusterTypeCountHashMap.begin(), clusterTypeCountHashMap.end());

  // Sort by ClusterType (you may need a custom comparator)
  sort(entries.begin(), entries.end(), [](const auto &a, const auto &b)
       {
         return a.first < b.first; // requires operator< for ClusterType
       });

  // Print
  for (const auto &entry : entries)
  {
    cout << entry.first << " : " << entry.second << endl;
  }

  string basisType = "Chebyshev";

  BasisSet atomicBasis(elementSet, basisType, true);

  //////////////////////////////////////////////////////////////////////

  auto cfg2 = Config::ReadPoscar("//media/sf_Phd/8.poscar");


  cfg2.UpdateNeighborList({2.85, 3.3, 4.7});

  cout << cfg2.GetNumAtoms() << endl;

  cout << cfg2.GetNeighborLatticeIdVectorOfLattice(0, 1).size() << endl;
  cout << cfg2.GetNeighborLatticeIdVectorOfLattice(0, 2).size() << endl;
  cout << cfg2.GetNeighborLatticeIdVectorOfLattice(0, 3).size() << endl;

  GetEncodeVectorCE(cfg2,
                    maxClusterSizeCE,
                    maxBondOrderCE,
                    clusterTypeCountHashMap,
                    initializedClusterHashSet);

  cout << cfg2.GetBasis() << endl;
}
*/
#include "Cluster.hpp"
#include "Structure.hpp"
#include "Config.h"
#include <pybind11/numpy.h>
#include <pybind11/embed.h> // Added for embedding Python interpreter
#include "InterfaceUtility.h"

#include "Geometry.h"

Config structure_to_config(const Structure &structure)
{
  // 1. Get cell (basis)
  Eigen::Matrix3d basis = structure.getCell();

  Eigen::MatrixXd positions = structure.getPositions();         // Nx3
  Eigen::Matrix3Xd relativePositionMatrix(3, positions.rows()); // 3 x N

  for (Eigen::Index i = 0; i < positions.rows(); ++i)
  {
    Eigen::Vector3d pos = positions.row(i);
    relativePositionMatrix.col(i) = pos;
  }

  // cout << relativePositionMatrix << endl;

  // Eigen::MatrixXd relativePositionMatrix = get_scaled_positions(
  //     positions_cart, basis, true, {true, true, true}); // wrap positions

  // 3. Get atomic numbers
  auto atomic_numbers_py = structure.getAtomicNumbers(); // py::array_t<int>
  auto buf = atomic_numbers_py.request();                // buffer info
  int *ptr = static_cast<int *>(buf.ptr);
  size_t n_atoms = buf.shape[0];

  std::vector<Element> atomVector(n_atoms);
  for (size_t i = 0; i < n_atoms; ++i)
  {
    atomVector[i] = Element(ptr[i]); // assumes Element has constructor from atomic number
  }

  // 4. Build Config
  return Config(basis, relativePositionMatrix, atomVector);
}

#include "Python.h"

/*
int main()
{
  py::scoped_interpreter guard{}; // <-- initializes Python

  auto cfg = Config::GenerateSupercell(4, 3.4, "Mo", "BCC");

  auto structure = ConvertConfigToStructure(cfg);

  cout << structure.getCell() << endl;

  auto prim = get_primitive_structure(structure, 0, 1, 1e-5);

  cout << prim.getCell() << endl;

  auto newCfg = structure_to_config(prim);

  Config::WriteConfig("/home/pravendra3/Desktop/test.cfg.gz", newCfg);

  cout << prim.size() << endl;


  vector<vector<string>> chemicalSymbols = {{"Mo", "Ta"}};



  auto decorated_prim_pair  = get_occupied_primitive_structure(
    prim, chemicalSymbols, 1e-5
  );


  auto newDecoratedStructure = decorated_prim_pair.first;
  auto primitiveChemicalSymbols = decorated_prim_pair.second;


  auto decoratedCfg = structure_to_config(newDecoratedStructure);

  Config::WriteConfig("/home/pravendra3/Desktop/test2.cfg.gz", decoratedCfg);


  print2DStringVector(primitiveChemicalSymbols);




}
  */
#include <iostream>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include "Structure.hpp"

namespace py = pybind11;

#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

namespace py = pybind11;
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For automatic conversion of Python lists/tuples to C++ vectors/tuples
#include <vector>
#include <tuple>
#include <iostream>
#include "Structure.hpp"
namespace py = pybind11;
#include "OrbitList.hpp"

#include "ClusterSpace.hpp"

/*
void printClusterSpace(const ClusterSpace& cs) {
    const auto& orbitList = cs.getOrbitList();
    std::cout << "Cluster Space\n";
    std::cout << "Total number of orbits: " << orbitList.size() << "\n";

    // Print basic info about each orbit
    for (size_t i = 0; i < orbitList.size(); ++i) {
        const auto& orbit = orbitList.getOrbit(i);
        std::cout << "Orbit " << std::setw(3) << i
                  << " | Order: " << orbit.getOrder()
                  << " | Radius: " << std::fixed << std::setprecision(4) << orbit.getRadius()
                  << " | Multiplicity: " << orbit.getMultiplicity()
                  << "\n";

        // Optional: print cluster sites
        std::cout << "  Sites: ";
        for (const auto& site : orbit.getRepresentativeCluster().getSites()) {
            std::cout << site.index << " ";
        }
        std::cout << "\n";
    }
}
*/

#include "ClusterExpansionCalculator.hpp"
#include "chrono"
#include <pybind11/pybind11.h>
#include <iostream>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "Structure.hpp"
#include "ClusterSpace.hpp"
#include "OrbitList.hpp"
#include "LocalOrbitListGenerator.hpp"
#include "VectorOperations.hpp"

#include <iostream>

double total_energy_from_formation(double E_form_per_atom, int x, int y, double muA, double muB)
{
  return E_form_per_atom * (x + y) + x * muA + y * muB;
}
#include <malloc.h>

void printMemoryUsage()
{
  malloc_stats(); // prints memory usage to stdout
}

void process_orbit_list_inputs(
    py::object standalone,
    py::object atoms_obj,
    py::list py_allowed,
    py::list py_cutoffs, double symprec)
{
  // Call GetInputsForOrbitList

  // py::module_ ase_atoms = py::module_::import("ase.atoms");
  // py::object Atoms = ase_atoms.attr("Atoms");

  auto prim_structure = standalone.attr("GetInputsForOrbitListPrimStructre")(
      atoms_obj, py_allowed, py_cutoffs, symprec);

  // const std::vector<std::vector<LatticeSite>> &matrixOfEquivalentSites,
  auto pm_lattice_sites_py = standalone.attr("GetInputsForOrbitListPMLatticeSites")(
      atoms_obj, py_allowed, py_cutoffs, symprec);

  // py::print(pm_lattice_sites_py);

  ///

  std::vector<std::vector<LatticeSite>> matrixOfEquivalentSites;

  for (auto row_obj : pm_lattice_sites_py)
  {
    std::vector<LatticeSite> row_vec;

    for (auto item : row_obj)
    {
      auto tuple_item = item.cast<py::tuple>(); // (int, tuple(int,int,int))
      int id = tuple_item[0].cast<int>();

      auto vec_tuple = tuple_item[1].cast<py::tuple>();
      Eigen::Vector3i vec3(vec_tuple[0].cast<int>(),
                           vec_tuple[1].cast<int>(),
                           vec_tuple[2].cast<int>());

      row_vec.emplace_back(id, vec3);
    }

    matrixOfEquivalentSites.push_back(std::move(row_vec));
  }

  // 1)), (0, (2, 0, -1)), (0, (1, 1, 3)), (0, (-1, -1, -3)), (0, (-1, -3, -1)), (0, (1, 3, 1)), (0, (-2, -2, -3)), (0, (2, 2, 3)), (0, (1, 0, -2)), (0, (-1, 0, 2)), (0, (2, -1, 0)), (0, (-2, 1, 0)), (0, (2, 2, 3)), (0, (-2, -2, -3)), (0, (0, -1, 2)), (0, (0, 1, -2)), (0, (1, -2, 0)), (0, (-1, 2, 0)), (0, (3, 1, 1)), (0, (-3, -1, -1)), (0, (-1, -1, -3)), (0, (1, 1, 3)), (0, (-3, -2, -2)), (0, (3, 2, 2)), (0, (-2, 1, 0)), (0, (2, -1, 0)), (0, (0, 2, -1)), (0, (0, -2, 1))]]

  // cout << matrixOfEquivalentSites << endl;

  // const std::vector<std::vector<std::vector<LatticeSite>>> &neighborLists
  auto neighbor_lists_py = standalone.attr("GetInputsForOrbitListNeighbourList")(
      atoms_obj, py_allowed, py_cutoffs, symprec);

  // py::print(neighbor_lists_py);
  // Conversion

  std::vector<std::vector<std::vector<LatticeSite>>> neighborLists;

  for (auto row_obj : neighbor_lists_py)
  {
    py::list row = row_obj.cast<py::list>();
    std::vector<std::vector<LatticeSite>> row_vec;
    for (auto col_obj : row)
    {
      py::list col = col_obj.cast<py::list>();
      std::vector<LatticeSite> col_vec;
      for (auto item : col)
      {
        auto tuple_item = item.cast<py::tuple>(); // (int, tuple(int,int,int))
        int id = tuple_item[0].cast<int>();

        auto vec_tuple = tuple_item[1].cast<py::tuple>();
        Eigen::Vector3i vec3(vec_tuple[0].cast<int>(),
                             vec_tuple[1].cast<int>(),
                             vec_tuple[2].cast<int>());

        col_vec.emplace_back(id, vec3);
      }
      row_vec.push_back(col_vec);
    }
    neighborLists.push_back(row_vec);
  }

  auto position_tolerance = standalone.attr("GetInputsForOrbitListTolerance")(
      atoms_obj, py_allowed, py_cutoffs, symprec);

  auto fractional_position_tolerance = standalone.attr("GetInputsForOrbitListFractionalTolerance")(
      atoms_obj, py_allowed, py_cutoffs, symprec);

  py::print(fractional_position_tolerance);
  py::print(position_tolerance);

  py::print(prim_structure);

  py::array_t<double> cell_py = prim_structure.attr("get_cell")();
  auto cell_buf = cell_py.unchecked<2>(); // 2D array

  Eigen::Matrix3d cell;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      cell(i, j) = cell_buf(i, j);

  std::cout << "Cell:\n"
            << cell << std::endl;

  py::array_t<double> pos_py = prim_structure.attr("get_positions")();
  auto pos_buf = pos_py.unchecked<2>();

  size_t n_atoms = pos_py.shape(0);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> positions(n_atoms, 3);

  for (size_t i = 0; i < n_atoms; i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      positions(i, j) = pos_buf(i, j);
    }
  }

  cout << positions << endl;

  py::list pbc_list = prim_structure.attr("get_pbc")();
  std::vector<bool> pbc(3);
  for (int i = 0; i < 3; i++)
    pbc[i] = pbc_list[i].cast<bool>();

  py::array_t<int> atomic_number_list = prim_structure.attr("get_atomic_numbers")();

  // py::print(atomic_number_list);

  Structure prim_structure_orbit_list(
      positions,
      atomic_number_list,
      cell,
      pbc);

  vector<Element> allowedElement = {Element("Mo"), Element("Ta")};

  std::vector<std::vector<int>> allowedAtomicNumbers;

  vector<int> atomicNumberVector;
  for (const auto &ele : allowedElement)
  {
    atomicNumberVector.emplace_back(ele.GetAtomicIndex());
  }

  allowedAtomicNumbers = {atomicNumberVector};

  prim_structure_orbit_list.setAllowedAtomicNumbers(allowedAtomicNumbers);

  OrbitList orbitList(
      prim_structure_orbit_list,
      matrixOfEquivalentSites,
      neighborLists,
      position_tolerance.cast<double>());

  orbitList.removeOrbitsWithInactiveSites();

  auto orbitListPtr = std::make_shared<OrbitList>(orbitList);
  ClusterSpace cs(orbitListPtr, position_tolerance.cast<double>(), fractional_position_tolerance.cast<double>());

  auto cfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");
  cfg.UpdateNeighborList({3, 4, 5});
  // cfg.SetElementOfLattice(0, Element("Ta"));

  /*
  for (size_t i = 0; i < 10; i ++ )
  {
    auto newCfg = cfg;
    newCfg.UpdateNeighborList({13});

    for (auto id : cfg.GetNeighborLatticeIdVectorOfLattice(i,1))
    {
      newCfg.SetElementOfLattice(id, Element("O"));
    }

    for (auto id : cfg.GetNeighborLatticeIdVectorOfLattice(i,2))
    {
      newCfg.SetElementOfLattice(id, Element("N"));
    }

    for (auto id : cfg.GetNeighborLatticeIdVectorOfLattice(i,3))
    {
      newCfg.SetElementOfLattice(id, Element("C"));
    }

    // newCfg.SetElementOfLattice(i, Element("O"));

    auto reducedMap = GetCenteredNeighboursSite(newCfg, i, 1);

    /*
    vector<Element> atomVector;

    Eigen::Matrix3Xd relativePositionMatrix(3, reducedMap.size());
    int k = 0;
    for (auto entry : reducedMap)
    {
      atomVector.emplace_back(newCfg.GetElementOfLattice(entry.first));

      relativePositionMatrix.col(k) = entry.second.transpose();

      cout << entry.second << endl;

      k++;
    }

    auto shiftedCfg = Config(
      newCfg.GetBasis(),
      relativePositionMatrix,
      atomVector
    );

    Config::WriteConfig("/home/pravendra3/Documents/icetClusters/shiftingToOrigin/shiftedEnv" + to_string(i) + ".cfg.gz", shiftedCfg);
  }
    */

  // exit(1);

  auto structureFromCfg = ConvertConfigToStructure(cfg);

  auto clusterVector = cs.getClusterVector(structureFromCfg, fractional_position_tolerance.cast<double>());

  print1DVector(clusterVector);

  cout << clusterVector.size() << endl;

  // i hope that it works

  auto supercell = std::make_shared<Structure>(structureFromCfg);

  LocalOrbitListGenerator LOLG = LocalOrbitListGenerator(cs.getPrimitiveOrbitList(), supercell, fractional_position_tolerance.cast<double>());
  // Create a full orbit list, used to calculate full cluster vectors.
  // _fullOrbitList = LOLG.getFullOrbitList();

  // Map indices in the supercell to offsets in the primitive cell,
  // and for each unique offset, precompute all possible
  // self-contained local orbit lists for the supercell.
  // The latter can be used to calculate local cluster vectors
  // and cluster vector differences.

  auto centralIdx = cfg.GetCentralAtomLatticeId();

  auto fractionalPositionTolerance = fractional_position_tolerance.cast<double>();

  // Find offset of this site in terms of the primitive structure
  Vector3d position = supercell->positionByIndex(centralIdx);
  Vector3i offset = cs.primitiveStructure().findLatticeSiteByPosition(position, fractionalPositionTolerance).unitcellOffset();

  // Create map from atom index to offset
  // _indexToOffset[i] = offset;

  cfg.UpdateNeighborList({13});

  /*
  auto nnOfCentralId = cfg.GetNeighborLatticeIdVectorOfLattice(centralIdx, 1);

  cout << nnOfCentralId.size() << endl;

  print1DVector(nnOfCentralId);

  auto newCfg = cfg;

  for (auto id : nnOfCentralId)
  {
    newCfg.SetElementOfLattice(id, Element("O"));
  }

  newCfg.SetElementOfLattice(centralIdx, Element("X"));
  Config::WriteConfig("/home/pravendra3/Documents/icetClusters/NearestNeighbours.cfg.gz", newCfg);
  */

  auto localOrbits = LOLG.getLocalOrbitList(offset, true);

  auto allOrbits = cs.getClusterVectorFromOrbitListCheck(localOrbits, supercell);

  VectorXd icetEci = ReadParametersFromJson("/home/pravendra3/Documents/LatticeMonteCarlo-eigen/icetECIs.json", "icet", "eci");

  // cluster Vector // Pure Mo

  double formationEnergy = 0;
  int i = 0;
  for (auto coef : icetEci)
  {
    formationEnergy += coef * clusterVector[i];
    i++;
  }
  cout << formationEnergy << endl;

  double mu_Mo = -10.93308598;
  double mu_Ta = -11.81233671;

  double E_total = total_energy_from_formation(formationEnergy, structureFromCfg.size(), 0, mu_Mo, mu_Ta);
  std::cout << "Formation Energy: " << formationEnergy << endl;
  std::cout << "Total energy of compound: " << E_total << " eV\n";

  cout << "Ta1Mo1999" << endl;

  auto centralLatticeId = cfg.GetCentralAtomLatticeId();
  cfg.SetElementOfLattice(centralLatticeId, Element("Ta"));

  auto structureFromCfg2 = ConvertConfigToStructure(cfg);

  auto clusterVector2 = cs.getClusterVector(structureFromCfg2, fractional_position_tolerance.cast<double>());

  double formationEnergy2 = 0;
  i = 0;
  for (auto coef : icetEci)
  {
    formationEnergy2 += coef * clusterVector2[i];
    i++;
  }
  cout << formationEnergy2 << endl;

  mu_Mo = -10.93308598;
  mu_Ta = -11.81233671;

  double E_total2 = total_energy_from_formation(formationEnergy2, structureFromCfg2.size() - 1, 1, mu_Mo, mu_Ta);
  std::cout << "Formation Energy: " << formationEnergy2 << endl;
  std::cout << "Total energy of compound: " << E_total2 << " eV\n";

  cout << "Adding Ta - pure Mo : " << E_total2 - E_total << endl;

  cout << "Local Cluster Expansion" << endl;

  auto localClusterVector1 = cs.getClusterVectorFromOrbitList(localOrbits, supercell);

  auto supercell2 = std::make_shared<Structure>(structureFromCfg2);

  auto localClusterVector2 = cs.getClusterVectorFromOrbitList(localOrbits, supercell2);

  double formationEnergyContribution1 = 0;
  double formationEnergyContribution2 = 0;

  i = 0;
  for (auto coef : icetEci)
  {
    formationEnergyContribution1 += coef * localClusterVector1[i];
    formationEnergyContribution2 += coef * localClusterVector2[i];
    i++;
  }

  cout << "Formation Energy for central site locally for Mo : " << formationEnergyContribution1 << endl;
  cout << "Formation Energy for central site locally for Ta : " << formationEnergyContribution2 << endl;

  cout << "dE formation energy : " << formationEnergyContribution2 - formationEnergyContribution1 << endl;

  double dELocal = 2000 * (formationEnergyContribution2 - formationEnergyContribution1) + (1999 - 2000) * mu_Mo + (1 - 0) * mu_Ta;

  cout << "dE Local : " << dELocal << endl;

  cout << "Local cluster vector check: " << endl;

  // lets do it for second

  cout << "From original function : " << endl;
  print1DVector(localClusterVector2);

  auto localClusterVectorCheck = cs.getLocalClusterVectorCheck(allOrbits, cfg);

  cout << "From Modified function : " << endl;
  print1DVector(localClusterVectorCheck);

  // Lets see what will be effect on the dE change during monte carlo simulaiton

  // cfg has Ta at center
  // migrate it to first nn simple

  // Define Cluster Space
  // Get the local orbit list for central site
  // Encode the clusters after canonically sorting the lattice Id vector which eventually can be stored
  // Then either one can use the icet's cs logic or since its bcc
  // One can directly have another function or class something like LocalClusterSpace
  // so this class can have a similar formulation as icet does or better is to
  // just have our own function i think both would work efficiently but
  // with lru the function stack can be reduced drastically
  // One more thing this LocalClusterSpace function will take the same argument as
  // ClusterSpace which is the primitive orbit just for the namesake
  // the one funciton which would be an exact copy of how ClusterSpace returns the localClusterVector
  // One funciton specific to cubic system or sysmter where the order of element cluster does not matter
  // due to space group symmetry which will basically use lru cache to reduce the redundant
  // function calls

  /*
  int j = 0;
  for (auto orbit : allOrbits)
  {
    // cout << "Orbit " << j << endl;

    // print2DVector(orbit);

    // auto newCfgLocal = newCfg;

    for (auto cluster : orbit)
    {
      for (auto id : cluster)
      {
        cout << " " << endl;
        // newCfgLocal.SetElementOfLattice(id, Element("N"));
      }
    }

    // Config::WriteConfig("/home/pravendra3/Documents/icetClusters/LocalOrbit" + to_string(j) + ".cfg.gz", newCfgLocal);

    j++;
  }
  */

  exit(1);

  ClusterExpansionCalculator ceCalculator(cs, structureFromCfg, fractional_position_tolerance.cast<double>());
  printMemoryUsage();

  /*
  double formationEnergy = 0;
  int i = 0;
  for (auto coef : icetEci)
  {
    formationEnergy += coef * clusterVector[i];
    i++;
  }
  cout << formationEnergy << endl;

  Config::WriteConfig("/home/pravendra3/Documents/LatticeMonteCarlo-eigen/testFromEnergy.cfg", cfg);

  double mu_Mo = -10.93308598;
  double mu_Ta = -11.81233671;

  double E_total = total_energy_from_formation(formationEnergy, 124, 4, mu_Mo, mu_Ta);
  std::cout << "Total energy of compound: " << E_total << " eV\n";

  auto occupationVector = structureFromCfg.getAtomicNumbers();

  py::print(occupationVector);
  /*
  auto ceClusterVector = ceCalculator.getClusterVector(occupationVector);

  formationEnergy = 0;
  i = 0;
  for (auto coef : icetEci)
  {
    formationEnergy += coef * ceClusterVector[i];
    i++;
  }
  cout << formationEnergy << endl;
  */

  auto structureCfg2 = ConvertConfigToStructure(cfg);

  auto occupationVector2 = structureCfg2.getAtomicNumbers();

  py::print(occupationVector2);
  /*
  // auto ceClusterVector2 = ceCalculator.getClusterVector(occupationVector2);
  formationEnergy = 0;
  i = 0;
  for (auto coef : icetEci)
  {
    formationEnergy += coef * ceClusterVector2[i];
    i++;
  }
  cout << formationEnergy << endl;

  cout << cfg.GetElementOfLattice(0).GetElementString() << endl;
  */

  /*
  auto testingCfg = Config::GenerateSupercell(3, 3.2, "Mo", "BCC");

  for (size_t i = 0; i < testingCfg.GetNumAtoms(); i++)
  {
    testingCfg.SetElementOfLattice(i, Element("Ta"));
    cout << i << endl;

    auto strTest = ConvertConfigToStructure(testingCfg);

    auto occVector = strTest.getAtomicNumbers();
    py::print(occVector);

    cout << "element at " << i << " : " << testingCfg.GetElementOfLattice(i) << endl;

    cout << "from structure : " << strTest.positionByIndex(i).transpose() << endl;
    cout << "from config : " << testingCfg.GetRelativePositionOfLattice(i).transpose() << endl;

    cout << "\n\n\n"
         << endl;
  }
         */
  auto start = std::chrono::high_resolution_clock::now();

  auto localCE = ceCalculator.getLocalClusterVector(structureFromCfg.getAtomicNumbers(), 0);

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

  std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;

  print1DVector(localCE);
  // cout << cs.to_string(50, 10) << endl;

  auto centralSiteIndex = cfg.GetCentralAtomLatticeId();

  cout << "Central Lattice Id: " << centralSiteIndex << endl;

  // ceCalculator.print(centralSiteIndex);

  vector<vector<size_t>> largest2DCluster = {
      {
          160,
          137,
      },
      {
          152,
          137,
      },
      {
          115,
          137,
      },
      {
          103,
          137,
      },
      {
          12,
          137,
      },
      {
          137,
          110,
      },
      {
          137,
          102,
      },
      {
          137,
          165,
      },
      {
          137,
          153,
      },
      {
          137,
          12,
      },
      {
          137,
          17,
      },
      {
          137,
          13,
      },
      {
          137,
          115,
      },
      {
          137,
          103,
      },
      {
          137,
          160,
      },
      {
          137,
          152,
      },
      {
          137,
          18,
      },
      {
          165,
          137,
      },
      {
          153,
          137,
      },
      {
          110,
          137,
      },
      {
          102,
          137,
      },
      {
          18,
          137,
      },
      {
          13,
          137,
      },
      {
          17,
          137,
      },
  };

  unordered_set<size_t> reqLatticeIds;

  for (auto cluster : largest2DCluster)
  {
    for (auto id : cluster)
    {
      reqLatticeIds.insert(id);
    }
  }

  auto cfgCopy = cfg;
  for (auto id : reqLatticeIds)
  {
    cfgCopy.SetElementOfLattice(static_cast<size_t>(id), Element("X"));
  }

  Config::WriteConfig("/home/pravendra3/Desktop/largestCluster.cfg", cfgCopy);

  /*
  OrbitList::OrbitList(const Structure &structure,
                       const std::vector<std::vector<LatticeSite>> &matrixOfEquivalentSites,
                       const std::vector<std::vector<std::vector<LatticeSite>>> &neighborLists,
                       const double positionTolerance) */
}

#include "SymmetricCEInterface.h"

  #include "EnergyPredictor.h"
  #include "ClusterExpansionParameters.h"
#include "SymmetricCE.h"

int main()
{
  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");
  vector<string> allowedElements = {"Mo", "Ta"};
  vector<double> clusterCutoffs = {13, 8, 6};

  auto supercellCfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");

  // Number of Ta atoms you want to add
  int numTa = 1000; // example

  // Get total number of sites
  int numSites = supercellCfg.GetNumAtoms();

  // Random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, numSites - 1);

  // Keep track of which sites are already replaced
  std::vector<size_t> replacedSites;

  // Randomly replace Mo with Ta
  for (int i = 0; i < numTa; ++i)
  {
    size_t site;
    do
    {
      site = dis(gen);
    } while (std::find(replacedSites.begin(), replacedSites.end(), site) != replacedSites.end());

    replacedSites.push_back(site);
    supercellCfg.SetElementOfLattice(site, Element("Ta")); // assuming such a function exists
  }

  vector<double> cutoffs = {3};
  supercellCfg.UpdateNeighborList(cutoffs);

  pair<size_t, size_t> latticeIdPair = {20, supercellCfg.GetNeighborLatticeIdVectorOfLattice(20, 1)[0]};
  // Get elements at the two sites
  string elem1 = supercellCfg.GetElementOfLattice(latticeIdPair.first).GetElementString();
  string elem2 = supercellCfg.GetElementOfLattice(latticeIdPair.second).GetElementString();

  // Ensure they are different
  if (elem1 == elem2)
  {
    if (elem1 == "Mo")
    {
      supercellCfg.SetElementOfLattice(latticeIdPair.second, Element("Ta"));
    }
    else
    {
      supercellCfg.SetElementOfLattice(latticeIdPair.second, Element("Mo"));
    }
  }

  // Print to check
  cout << supercellCfg.GetElementOfLattice(latticeIdPair.first).GetElementString() << endl;
  cout << supercellCfg.GetElementOfLattice(latticeIdPair.second).GetElementString() << endl;

  /*
  // Around Site 20
  for (auto id : supercellCfg.GetNeighborLatticeIdVectorOfLattice(latticeIdPair.first, 1))
  {
    supercellCfg.SetElementOfLattice(id, Element("Ta"));
  }

  // Around Site 100
  int i = 0;
  for (auto id : supercellCfg.GetNeighborLatticeIdVectorOfLattice(latticeIdPair.second, 1))
  {
    if (i < 4)
    {
      supercellCfg.SetElementOfLattice(id, Element("Ta"));
    }
    i += 1;
  }

  supercellCfg.SetElementOfLattice(latticeIdPair.second, Element("Ta"));
  */


  vector<double> maxCutoffVector = {
      *max_element(clusterCutoffs.begin(), clusterCutoffs.end())};
  supercellCfg.UpdateNeighborList(maxCutoffVector);


  VectorXd eciVector = ReadParametersFromJson("/home/pravendra3/Documents/LatticeMonteCarlo-eigen/icetECIs.json", "ce", "eci");

  std::vector<double> eciStdVector;
  eciStdVector.reserve(eciVector.size());
  for (int i = 0; i < eciVector.size(); ++i)
  {
    eciStdVector.push_back(eciVector[i]);
  }



  // Lets Choose Lattice Ids
  /*
  {SymmetricCE symCE(
      supercellCfg,
      primCfg,
      allowedElements,
      clusterCutoffs);


  symCE.TestSymmetricCE(supercellCfg,
                        eciStdVector,
                        latticeIdPair);

  print1DVector(symCE.GetClusterVectorForConfig(supercellCfg));}
  */
  
  string predictorFile = "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/icetECIs.json";
  SymmetricCEPredictor energyPredictor(
    predictorFile, 
    supercellCfg, 
    primCfg, 
    allowedElements, 
    clusterCutoffs
  );

  cout << "Energy of Config : " << energyPredictor.ComputeEnergyOfConfig(supercellCfg) << endl;

  auto dE = energyPredictor.GetDeSwap(supercellCfg, 
  latticeIdPair);

  cout << "Energy change : " << dE << endl;

  exit(1);

  // print1DVector(clusterVector);

  py::scoped_interpreter guard{}; // start Python interpreter

  try
  {
    // Add Python utility path
    py::module_ sys = py::module_::import("sys");
    sys.attr("path").attr("append")("/home/pravendra3/Documents/LatticeMonteCarlo-eigen/lmc/utility/src/");

    // Import Python modules
    py::module_ standalone = py::module_::import("standalone");
    py::module_ ase_atoms = py::module_::import("ase.atoms");

    // Create C++ structure
    auto cfg = Config::GenerateSupercell(4, 3.2, "Mo", "BCC");
    auto structure = ConvertConfigToStructure(cfg);

    // Convert Eigen matrices to NumPy arrays
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> pos = structure.getPositions();
    py::array_t<double> positions({pos.rows(), pos.cols()}, pos.data());

    Eigen::Matrix3d c = structure.getCell();
    py::array_t<double> cell({3, 3}, c.data());

    py::array_t<int> numbers = structure.getAtomicNumbers();

    // Convert PBC
    py::list pbc;
    for (bool b : structure.getPBC())
      pbc.append(b);

    // Create ASE Atoms object
    py::object Atoms = ase_atoms.attr("Atoms");
    py::object atoms_obj = Atoms(
        py::arg("positions") = positions,
        py::arg("numbers") = numbers,
        py::arg("cell") = cell,
        py::arg("pbc") = pbc);

    // Prepare allowed_species as py::list of py::list
    std::vector<std::string> allowed = {"Mo", "Ta"}; // example
    py::list py_allowed;
    for (auto &sub : allowed)
    {
      py_allowed.append(py::str(sub));
    }

    // Prepare cutoffs
    vector<double> cutoffs = {13, 8, 5};

    py::list py_cutoffs;

    for (auto cutoff : cutoffs)
    {
      py_cutoffs.append(cutoff);
    }

    double symprec = 1e-5;

    // Call the Python function
    py::list result1 = standalone.attr("get_chemical_symbols")(
        py_allowed, atoms_obj);

    // Convert to C++ vector<vector<string>>
    std::vector<std::vector<std::string>> chemical_symbols;
    for (auto row_obj : result1)
    {
      py::list row = row_obj.cast<py::list>();
      std::vector<std::string> row_vec;
      for (auto item : row)
      {
        row_vec.push_back(item.cast<std::string>());
      }
      chemical_symbols.push_back(row_vec);
    }

    // print2DStringVector(chemical_symbols);

    // Call GetInputsForOrbitList
    process_orbit_list_inputs(standalone, atoms_obj, py_allowed, py_cutoffs, symprec);
  }
  catch (const py::error_already_set &e)
  {
    std::cerr << "Python error: " << e.what() << std::endl;
    return 1;
  }
  catch (const std::exception &e)
  {
    std::cerr << "C++ error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
