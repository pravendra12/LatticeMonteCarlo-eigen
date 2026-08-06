#include "LocalOrbitListGenerator.hpp"

/**
@param primitiveOrbitList orbit list for the underlying primitive cell
@param supercell supercell structure for which to set up the local orbit list generation
@param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates
*/
LocalOrbitListGenerator::LocalOrbitListGenerator(const OrbitList &primitiveOrbitList,
                                                 std::shared_ptr<Structure> supercell,
                                                 const double fractionalPositionTolerance)
    : _primitiveOrbitList(primitiveOrbitList),
      _supercell(supercell),
      _fractionalPositionTolerance(fractionalPositionTolerance)
{
  _primitiveStructure = std::make_shared<Structure>(primitiveOrbitList.structure());

  // Extract unique primitive cell offsets, i.e., all offsets in units of the
  // primitive lattice vectors that together tile the supercell
  for (size_t i = 0; i < _supercell->size(); i++)
  {
    // Find the lattice site of the supercell site in terms of the primitive cell
    Vector3d position_i = _supercell->positionByIndex(i);

    LatticeSite primitiveSite = _primitiveStructure->findLatticeSiteByPosition(position_i, _fractionalPositionTolerance);

    // We only add the cell offsets for atoms with primitive index 0.
    // Atoms with non-zero index may sit in cells with offsets that
    // are not covered here. Still, such offsets are not added here,
    // because here we only want to make sure that we have sufficient
    // offsets to tile the supercell with primitive cells. This means that
    // getLocalOrbitList may be called with an offset that is not in
    // uniqueCellOffets, but that is not a problem.
    if (primitiveSite.index() == 0)
    {
      _uniquePrimitiveCellOffsets.push_back(primitiveSite.unitcellOffset());
    }
  }
  if (_uniquePrimitiveCellOffsets.size() != _supercell->size() / _primitiveStructure->size())
  {
    std::ostringstream msg;
    msg << "Wrong number of unitcell offsets found (LocalOrbitListGenerator constructor)." << std::endl;
    msg << "Expected: " << _supercell->size() / _primitiveStructure->size() << std::endl;
    msg << "Found:    " << _uniquePrimitiveCellOffsets.size();
    throw std::runtime_error(msg.str());
  }
  std::sort(_uniquePrimitiveCellOffsets.begin(), _uniquePrimitiveCellOffsets.end(), Vector3iCompare());
}

/**
@details Generates and returns a local orbit for the supplied offset.
@param offset Cell offset in multiples of primitive lattice vectors
@param selfContained
    If this orbit list is used on its own to calculate local cluster vectors or
    differences in cluster vector, this parameter needs to be true (if false, not all
    clusters involving this offset will be included).
*/
OrbitList LocalOrbitListGenerator::getLocalOrbitList(const Vector3i &offset, bool selfContained = false)
{
  OrbitList localOrbitList = OrbitList();
  for (const Orbit &orbit : _primitiveOrbitList.orbits())
  {
    // Copy the orbit.
    Orbit supercellOrbit = orbit;

    // If this orbit list will be used standalone for calculating
    // local cluster vectors or cluster vector differences,
    // we need to add clusters that include the present cell offset,
    // but would otherwise belong to the local orbit list of
    // another cell offset.
    if (selfContained)
    {
      // We will loop over the clusters in the orbit and add more
      // clusters inside the loop, so we first extract the original
      // clusters to avoid modifying the list we are looping over.
      std::vector<Cluster> clusters = supercellOrbit.clusters();
      for (auto cluster : clusters)
      {

        // Extract all versions of the clusters for which the
        // original cluster has been translated such that one
        // of the sites sits in the {0, 0, 0} cell offset.
        std::vector<std::vector<LatticeSite>> translatedSiteGroups = _primitiveOrbitList.getSitesTranslatedToUnitcell(cluster.latticeSites(), false);
        for (auto translatedSites : translatedSiteGroups)
        {
          // Only add clusters that are not duplicates of previus clusters.
          if (!supercellOrbit.contains(translatedSites))
          {
            supercellOrbit.addCluster(Cluster(translatedSites, _primitiveStructure));
          }
        }
      }
    }

    // Translate all clusters of the new orbit.
    supercellOrbit.translate(offset);

    // Technically we should use the fractional position tolerance
    // corresponding to the cell metric of the supercell structure.
    // This is, however, not uniquely defined. Moreover, the difference
    // would only matter for very large supercells. We (@angqvist,
    // @erikfransson, @erhart) therefore decide to defer this issue
    // until someone encounters the problem in a practical situation.
    // In principle, one should not handle coordinates (floats) at this
    // level anymore. Rather one should transform any (supercell)
    // structure into an effective representation in terms of lattice
    // sites before any further operations.
    supercellOrbit.transformToSupercell(_supercell, _primitiveToSupercellMap, _fractionalPositionTolerance);
    localOrbitList.addOrbit(supercellOrbit);
  }
  return localOrbitList;
}

/// Generates the complete orbit list (the sum of all local orbit lists).
OrbitList LocalOrbitListGenerator::getFullOrbitList()
{
  OrbitList orbitList = OrbitList();

  // Generate the full orbit list by adding together the local orbit lists for each unique offset
  for (auto &offset : _uniquePrimitiveCellOffsets)
  {
    orbitList += getLocalOrbitList(offset);
  }

  if (orbitList.size() != _primitiveOrbitList.size())
  {
    std::ostringstream msg;
    msg << "Full orbitlist size is not the same as local orbitlist size (LocalOrbitListGenerator::getFullOrbitList)" << std::endl;
    msg << " full orbitlist size: " << orbitList.size() << std::endl;
    msg << " local orbitlist size: " << _primitiveOrbitList.size() << std::endl;
    throw std::runtime_error(msg.str());
  }
  return orbitList;
}
