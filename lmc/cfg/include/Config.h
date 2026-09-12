/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 1/16/20 3:55 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/30/23 3:10 PM                                                          *
 **************************************************************************************************/

/*! \file  Config.h
 *  \brief File for the Config class definition.
 */

#ifndef LMC_CFG_INCLUDE_CONFIG_H_
#define LMC_CFG_INCLUDE_CONFIG_H_

#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <utility>
#include <variant>
#include <any>
#include "Eigen/Dense"
#include "Constants.hpp"
#include "Element.hpp"
#include <utility>

#include <stdexcept>
#include <optional>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
// #include "Atom.hpp"

using namespace std;

/*! \brief Class for defining a configuration of atoms and their positions.
 */
class Config
{
public:
  /*! \brief Default constructor for Config.
   */
  Config();

  /*! \brief Constructor for setting up the configuration of atoms and their positions.
   *  \param basis                     The basis vectors (3x3) of the configuration.
   *  \param relative_position_matrix  The relative position matrix (3，n) of the configuration.
   *  \param atom_vector               The atom vector of the configuration，n atoms in total.
   */
  Config(Eigen::Matrix3d basis, Eigen::Matrix3Xd relative_position_matrix, std::vector<Element> atom_vector);

  /*! \brief Default destructor for Config.
   */
  virtual ~Config();

  /*! \brief Query for the number of atoms in the configuration.
   *  \return  The number of atoms in the configuration.
   */
  [[nodiscard]] size_t GetNumAtoms() const;

  /*! \brief Query for the vector of atoms in the configuration.
   *  \return  The vector of atoms in the configuration.
   */
  [[nodiscard]] const std::vector<Element> &GetAtomVector() const;

  /*! \brief Query for the number of lattice sites in the configuration.
   *  \return  The number of lattice sites in the configuration.
   */
  [[nodiscard]] size_t GetNumLattices() const;

  /*! \brief Query for the basis vectors of the configuration.
   *  \return  The basis vectors of the configuration.
   */
  [[nodiscard]] const Eigen::Matrix3d &GetBasis() const;

  /*! \brief Query for the lists of neighbor sites of each site in the configuration.
   *  \return  The lists of neighbor sites of each site in the configuration.
   */
  [[nodiscard]] const std::vector<std::vector<std::vector<size_t>>> &GetNeighborLists() const;

  /*! \brief Query for the lattice ID of vacancy for the configuration.
   *  \return  The vacancy lattice ID site for the configuration.
   */
  [[nodiscard]] size_t GetVacancyLatticeId() const;

  /*! \brief Returns the relative position matrix for the configuration
   */
  [[nodiscard]] const Eigen::Matrix3Xd &GetRelativePositionMatrix() const;

  [[nodiscard]] size_t GetCentralAtomLatticeId() const;

  /*! \brief Query for the Atom ID of vacancy for the configuration.
   *  \return  The vacancy Atom ID site for the configuration.
   */
  [[nodiscard]] size_t GetVacancyAtomId() const;

  /*! \brief Returns the concentration corresponding to each element present in
   *         in the configuration.
   *  \return  The concentration of elements.
   */
  [[nodiscard]] std::map<Element, double> GetConcentration() const;

  /*! \brief Returns the element and their atom ID at which they are present.
   *  \return  Element and Atom ID vectors at which those element are present.
   */
  [[nodiscard]] std::map<Element, std::vector<size_t>>
  GetElementOfAtomIdVectorMap() const;

  /*! \brief Occupancy vector contains the atomic number of element present at each lattice Id
      \return Returns the occupancy vector sorted by the lattice Id
  */
  vector<size_t> GetOccupationVector() const;

  /*! \brief Query for the set of neighbor atom id of an atom.
   *  \param atom_id         The atom id of the atom.
   *  \param distance_order  The order of distance between the two lattice.
   *  \return                The vector of neighbor atom IDs.
   */
  [[nodiscard]] std::vector<size_t> GetNeighborAtomIdVectorOfAtom(size_t atom_id, size_t distance_order) const;

  /*! \brief Query for the set of neighbor lattice id of a lattice.
   *  \param lattice_id      The lattice id of the lattice.
   *  \param distance_order  The order of distance between the two lattice.
   *  \return                The vector of neighbor lattice IDs .
   */
  [[nodiscard]] std::vector<size_t> GetNeighborLatticeIdVectorOfLattice(size_t lattice_id, size_t distance_order) const;

  /*! \brief Query for the set of neighbor lattice id of a lattice.
   *  \param lattice_id      The lattice id of the lattice.
   *  \param distance_order  The order of distance between the two lattice.
   *  \return                The vector of neighbor lattice IDs .
   */
  [[nodiscard]] vector<size_t> GetNeighborLatticeIdsUpToOrder(size_t latticeId, size_t maxBondOrder) const;

  /*! \brief Query for the order of distance between two lattice.
   *  \param lattice_id1  The lattice id of the first lattice.
   *  \param lattice_id2  The lattice id of the second lattice.
   *  \return             The order of distance between the two lattice. 0 means the same lattice or not neighbors.
   */
  [[nodiscard]] size_t GetDistanceOrder(size_t lattice_id1, size_t lattice_id2) const;

  /*! \brief Query for the element of a lattice site.
   *  \param lattice_id  The lattice id of the lattice site.
   *  \return            The element of the lattice site.
   */
  [[nodiscard]] Element GetElementOfLattice(size_t lattice_id) const;

  /*! \brief Query for the element of an atom.
   *  \param atom_id  The atom id of the atom.
   *  \return         The element of the atom.
   */
  [[nodiscard]] Element GetElementOfAtom(size_t atom_id) const;

  /*! \brief Query for atom id of a lattice.
   *  \param latticeId  The lattice id of the lattice.
   *  \return         The atom id for the atom which is at that lattice.
   */
  [[nodiscard]] size_t GetAtomIdOfLattice(size_t latticeId) const;

  /*! \brief Query for lattice id of the atom.
   *  \param atomId  The atom id of the atom.
   *  \return         The lattice id for the atom.
   */
  [[nodiscard]] size_t GetLatticeIdOfAtom(size_t atomId) const;

  /*! \brief Query for the cartesian position of a lattice site. //prav: relative position of lattice
   *  \param lattice_id  The lattice id of the lattice site.
   *  \return            The cartesian position of the lattice site.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetRelativePositionOfLattice(size_t lattice_id) const;

  /*! \brief Query for the relative position of an atom.
   *  \param atom_id  The atom id of the atom.
   *  \return         The relative position of the atom.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetRelativePositionOfAtom(size_t atom_id) const;

  /*! \brief Query for the cartesian position of an lattice site.
   *  \param lattice_id  The lattice id of the lattice site.
   *  \return            The cartesian position of the lattice site.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetCartesianPositionOfLattice(size_t lattice_id) const;

  /*! \brief Query for the cartesian position of an atom.
   *  \param atom_id  The atom id of the atom.
   *  \return         The cartesian position of the atom.
   */
  [[nodiscard]] Eigen::Ref<const Eigen::Vector3d> GetCartesianPositionOfAtom(size_t atom_id) const;

  /// Functions which are used in Symmetry Operations.

  std::vector<size_t> GetSortedLatticeVectorStateOfPair(
      const std::pair<size_t, size_t> &lattice_id_pair, const size_t &max_bond_order) const;

  // Returns the sorted lattice vector including the lattice pair
  std::vector<size_t> GetSortedLatticeVectorStateWithPair(
      const std::pair<size_t, size_t> &lattice_id_jump_pair,
      const size_t &max_bond_order) const;

  /*! \brief Computes the center position of a lattice pair while accounting for
               periodic boundary conditions.
     *
     * This function calculates the geometric center of two lattice points specified
     * by their IDs. It adjusts the relative positions of the lattice points to
     * ensure that the computed distance between them is within the range (0, 0.5)
     * in each dimension, considering periodic boundaries.
     *
     *  \param config               A reference to the Config object containing
     *                              lattice configurations and relative positions.
     *  \param lattice_id_jump_pair A pair of lattice IDs representing the two
     *                              lattice points.
     *  \return                     Eigen::Vector3d The computed center position of
     *                              the lattice pair in three dimensions.
     */
  [[nodiscard]] Eigen::RowVector3d GetLatticePairCenter(
      const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

  static void RotateLatticeVector(std::unordered_map<size_t, Eigen::RowVector3d> &lattice_id_hashmap,
                                  const Eigen::Matrix3d &rotation_matrix);

  [[nodiscard]] Eigen::Matrix3d GetLatticePairRotationMatrix(const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

  /*! \brief Retrieves the neighboring lattice IDs for a specified lattice ID pair
   *         up to the maximum bond order.
   *
   *  \param lattice_id_pair  The pair of lattice IDs whose neighbors are to be
   *                          found.
   *  \param max_bond_order   The maximum bond order to consider for neighbor
   *                          identification.
   *  \return                 A set containing the neighboring lattice IDs.
   */
  [[nodiscard]] std::unordered_set<size_t> GetNeighboringLatticeIdSetOfPair(
      const std::pair<size_t, size_t> &lattice_id_pair, const size_t &max_bond_order) const;

  [[nodiscard]] Eigen::Vector3d GetNormalizedDirection(size_t referenceId, size_t latticeId) const;

  /*! \brief Query for the relative distance vector between two lattice.
   *  \param lattice_id1  The lattice id of the first lattice.
   *  \param lattice_id2  The lattice id of the second lattice.
   *  \return             The relative distance vector between the two lattice.
   */
  [[nodiscard]] Eigen::Vector3d GetRelativeDistanceVectorLattice(size_t lattice_id1, size_t lattice_id2) const;

  /*! \brief Set the periodic boundary condition of the configuration.
   *  \param periodic_boundary_condition  The periodic boundary condition of the configuration.
   */
  void SetPeriodicBoundaryCondition(const std::array<bool, 3> &periodic_boundary_condition);

  /*! \brief Wrap atoms outside the unit cell into the cell.
   */
  void Wrap();

  /*! \brief Set the element type at the atom with given atom id.
   *  \param atom_id       The atom id of the atom.
   *  \param element_type  The new element type.
   */

  void SetElementOfAtom(size_t atom_id, Element element_type);

  void SetElementOfLattice(size_t lattice_id, Element element_type);

  /*! \brief Modify the atom configuration.
   *  \param lattice_id_jump_pair  The pair of lattice ids to modify the configuration.
   */
  void LatticeJump(const std::pair<size_t, size_t> &lattice_id_jump_pair);

  /*! \brief Update the neighbor list of the configuration with the given cutoffs.
   *  \param cutoffs  The cutoffs to update the neighbor list.
   */
  void UpdateNeighborList(std::vector<double> cutoffs);

  /*! \brief Read the configuration from a lattice file, element file and map file.
   *  \param lattice_filename  The name of the lattice file.
   *  \param element_filename  The name of the element file.
   *  \param map_filename      The name of the map file.
   *  \return                  The configuration read from the file.
   */
  static Config ReadMap(const std::string &lattice_filename,
                        const std::string &element_filename,
                        const std::string &map_filename);

  /*! \brief Read the configuration from a CFG file.
   *  \param filename  The name of the CFG file.
   *  \return          The configuration read from the file.
   */
  static Config ReadCfg(const std::string &filename);

  /*! \brief Read the configuration from a POSCAR file.
   *  \param filename  The name of the POSCAR file.
   *  \return          The configuration read from the file.
   */
  static Config ReadPoscar(const std::string &filename);

  /*! \brief Read the configuration from a POSCAR or CFG file.
   *  \param filename  The name of the configuration file.
                        Supported formats are: .cfg, .POSCAR, .cfg.gz,
                        .cfg.bz2, .POSCAR.gz, .POSCAR.bz2.
   *  \return          The configuration read from the file.
   */
  static Config ReadConfig(const std::string &filename);

  /*! \brief Read the configuration from a XYZ file.
   *  \param filename  The name of the configuration file.
   *  \return          The configuration read from the file.
   */
  static Config ReadXyz(const std::string &filename);

  /*! \brief Generate the supercell configuration.
   *  \param supercell_size  Size of the Supercell.
   *  \param lattice_param   Lattice Parameter of the unit cell.
   *  \param element_symbol  Element Symbol eg. "Al".
   *  \param structure_type  Type of Structure  "BCC" or "FCC".
   *  \return                The configuration of the generated structure.
   */
  static Config GenerateSupercell(size_t supercell_size, double lattice_param, const std::string &element_symbol, const std::string &structure_type);

  static Config GenerateAlloySupercell(
      size_t supercell_size,
      double lattice_param,
      std::string structure_type,
      const std::vector<std::string> &element_vector,
      const std::vector<double> &composition_vector,
      unsigned seed);

  /*! \brief Write the configuration to a file.
   *  \param filename  The name of the file to write the configuration to.
   */
  static void WriteConfig(const std::string &filename,
                          const Config &config_out);

  /*! \brief Write the extended configuration to a file. Check the extended CFG format
   *         at <http://li.mit.edu/Archive/Graphics/A/#extended_CFG>.
   *  \param filename         The name of the file to write the extended configuration to.
   *  \param auxiliary_lists  The auxiliary lists to write to the file.
   */
  static void WriteConfigExtended(
      const std::string &filename,
      const Config &config_out,
      const std::map<std::string, std::vector<double>> &auxiliary_lists);

  using VectorVariant = std::variant<std::vector<int>,
                                     std::vector<size_t>,
                                     std::vector<double>,
                                     std::vector<std::string>,
                                     std::vector<Eigen::Vector3d>>;
  using ValueVariant = std::variant<int, double, unsigned long long, std::string>;

  /*! \brief Write the extended XXY to a file. Check the extended XYZ format
   *         at <https://web.archive.org/web/20190811094343/https://libatoms.github.io/QUIP/io.html#extendedxyz>
   *         and <https://www.ovito.org/docs/current/reference/file_formats/input/xyz.html#file-formats-input-xyz>
   *  \param filename         The name of the file to write the extended configuration to.
   *  \param auxiliary_lists  The lists of atom information .
   *  \param global_list      The list of global information.
   */
  static void WriteXyzExtended(
      const std::string &filename,
      const Config &config_out,
      const std::map<std::string, VectorVariant> &auxiliary_lists,
      const std::map<std::string, ValueVariant> &global_list);

  // /*! \brief Write the lattice configuration to a file.
  //  *         Outputs a lattice configuration with bonds up to the specified maximum order.
  //  *         This function writes the lattice in an easily readable format, storing
  //  *         bond orders and other lattice properties for further analysis or visualization.
  //  *  \param filename       The name of the file to write the lattice configuration to.
  //  *  \param max_bond_order: The maximum bond order to include in the written configuration.
  //  */
  void WriteLattice(const std::string &filename, size_t &max_bond_order) const;

  Config ExtractLocalSupercell(const std::pair<size_t, size_t> &latticeJumpPair, size_t supercellSize, double latticeParam) const;
  Config ExtractLocalSupercellGrok(const std::pair<size_t, size_t> &latticeJumpPair, size_t supercellSize, double latticeParam) const;

private:
  /*! \brief Sort lattice sites by the positions (x, y, z)
   */
  void ReassignLattice();

  /// The periodic boundary condition status in all three directions.
  std::array<bool, 3> periodic_boundary_condition_{true, true, true};

  /// The basis matrix of the configuration.
  Eigen::Matrix3d basis_{};

  /// The relative position matrix of the configuration.
  Eigen::Matrix3Xd relative_position_matrix_{};

  /// The cartesian position matrix of the configuration.
  Eigen::Matrix3Xd cartesian_position_matrix_{};

  /// The vector of atoms in the configuration.
  std::vector<Element> atom_vector_{};

  /// Mapping from lattice points to atom ids.
  std::unordered_map<size_t, size_t> lattice_to_atom_hashmap_{};

  /// Mapping from atom ids to lattice points.
  std::unordered_map<size_t, size_t> atom_to_lattice_hashmap_{};

  /// Nearest neighbor lists. The first key is the cutoff distance between the two lattice sites in Angstrom.
  std::vector<std::vector<std::vector<size_t>>> neighbor_lists_{};

  /// Cutoffs for the neighbor lists.
  std::vector<double> cutoffs_{};

  /// Number of cells in each direction.
  Eigen::Vector3i num_cells_{};

  /// The cells of the configuration. Each cell is a vector of lattice ids.
  std::vector<std::vector<size_t>> cells_{};
};

#endif // LMC_CFG_INCLUDE_CONFIG_H_
