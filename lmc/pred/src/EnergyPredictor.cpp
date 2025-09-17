#include "EnergyPredictor.h"

EnergyPredictor::EnergyPredictor(
    SymmetricCEPredictor &symCEEnergyPredictor,
    LVFEPredictor &lvfePredictor) : symCEEnergyPredictor_(symCEEnergyPredictor),
                                    lvfePredictor_(lvfePredictor),
                                    allowedElements_(
                                        symCEEnergyPredictor_.GetAllowedElements())
{
}

double EnergyPredictor::GetEnergyChange(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  auto firstElement = config.GetElementOfLattice(latticeIdJumpPair.first);
  auto secondElement = config.GetElementOfLattice(latticeIdJumpPair.second);

  if (firstElement == secondElement)
  {
    return 0;
  }

  bool firstIsVacancy = (firstElement == Element("X"));
  bool secondIsVacancy = (secondElement == Element("X"));

  double dE = 0;

  if (firstIsVacancy || secondIsVacancy)
  {
    // Any of the two site is vacancy
    dE = GetEnergyChangeWithVacancy(
        config,
        latticeIdJumpPair);
  }
  else
  {
    // Both are atoms
    dE = symCEEnergyPredictor_.GetDeSwapConst(
        config,
        latticeIdJumpPair);
  }

  return dE;
}

double EnergyPredictor::GetEnergyChangeWithVacancy(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  /*
    Local Vacancy Formation Energy + Symmetric CE contribution for a pair of sites:

     Compute total energy change due to a vacancy jump:

     dEv = dELVFE + (1/N) * sum_{s ∈ species} dEs

     where:
       dEv    : total energy change for the vacancy jump
       dELVFE : change in local vacancy formation energy
       dEs    : energy change contribution from species s (from symCE)
       N      : total number of allowed species (size of allowedElements_)

     Note:
       - If s is the migrating species, its contribution dEs = 0
       - This generalizes the binary case where we had 1/2 * (dEMo + dETa)
       - For three species, this becomes 1/3 * sum(dEs for all species)


    For binary system applicable to multicomponent system, averaging need to be done
    for all the elements except vacancy.
    Let:
      Ev^1  = Energy of config with vacancy at site 1
               Ev^1 = ELVFE^1 + 1/2 *(EMo^1 + ETa^1)
      Ev^2  = Energy of config with vacancy at site 2
               Ev^2 = ELVFE^2 + 1/2 *(EMo^2 + ETa^2)

    Here:
      Ev   : Energy of configuration with vacancy
      ELVFE: Local Vacancy Formation Energy (from fitting)
      EMo  : Energy of configuration when Mo is assigned to the vacancy site
      ETa  : Energy of configuration when Ta is assigned to the vacancy site
      Superscript 1 or 2 indicates site index.

    Energy difference between configurations:
      dEv = Ev^2 - Ev^1
           = (ELVFE^2 - ELVFE^1) + 1/2 * ((EMo^2 - EMo^1) + (ETa^2 - ETa^1))

    Explanation for EMo (same applies to ETa):
      - Suppose site 1 has a vacancy and site 2 has the migrating atom.
      - If the migrating atom is Mo:
          * Before jump:
              - EMo is computed for site 1 (vacancy site) for whole config
          * After jump:
              - Site 1 now has Mo (migrated atom)
              - EMo for site 2 (previously Mo) is identical to EMo before
              - Therefore, the change in EMo contribution is zero
      - This logic ensures that if the migrating element is the same as the atom at the other site,
        the energy contribution for that species does not artificially change.

    General rule:
      - Compute dEv as the sum of:
          1. Change in local vacancy formation energy (dELVFE)
          2. Half the sum of changes in Mo and Ta energies (dEMo + dETa)
      - Contributions for a given element are zero if the migrating atom is of the same type,
        because the occupancy swap does not change the total energy for that element.

    Benefit:
      - Accurately accounts for local energy changes due to vacancy jumps
      - Avoids double counting the energy of migrating atoms already present
  */

  size_t vacancyLatticeId;
  size_t migratingAtomLatticeId;

  bool firstIsVacancy = (config.GetElementOfLattice(latticeIdJumpPair.first) == Element("X"));
  bool secondIsVacancy = (config.GetElementOfLattice(latticeIdJumpPair.second) == Element("X"));

  if (firstIsVacancy && !secondIsVacancy)
  {
    vacancyLatticeId = latticeIdJumpPair.first;
    migratingAtomLatticeId = latticeIdJumpPair.second;
  }
  else if (!firstIsVacancy && secondIsVacancy)
  {
    vacancyLatticeId = latticeIdJumpPair.second;
    migratingAtomLatticeId = latticeIdJumpPair.first;
  }
  else
  {
    throw std::runtime_error("Error in `EnergyPredictor::GetEnergyChange`: One of the lattice sites must be a vacancy.");
  }

  auto migratingElement = config.GetElementOfLattice(migratingAtomLatticeId);

  // This give the dELVFE term
  // This applies for both migration and swap because
  // before jump the config is used to get the elements
  // After jump the migrating element will be specifically set
  // for the other site.
  double dElvfe = lvfePredictor_.GetDeForVacancyMigration(
      config,
      latticeIdJumpPair);

  // Iterate over allowedElements_ vector and compute the dE using
  // symCEEnergyPredictor_ for elements other than migrating elements

  double dESum = 0;

  for (const auto &eleString : allowedElements_)
  {
    if (eleString == migratingElement.GetElementString())
    {
      // dE_migrating_element will be 0
      /*
      - Suppose site 1 has a vacancy and site 2 has the migrating atom.
      - If the migrating atom is Mo:
          * Before jump:
              - EMo is computed for site 1 (vacancy site) for whole config
          * After jump:
              - Site 1 now has Mo (migrated atom)
              - EMo for site 2 (previously Mo) is identical to EMo before
              - Therefore, the change in EMo contribution is zero
      - This logic ensures that if the migrating element is the same as the atom at the other site,
        the energy contribution for that species does not artificially change.
      */
      continue;
    }

    // From other element compute the dE
    double dEELement = symCEEnergyPredictor_.GetDeSwapConst(
        config,
        make_pair(vacancyLatticeId, migratingAtomLatticeId),
        make_pair(Element(eleString), migratingElement));

    dESum += dEELement;
  }

  double dEAverage = dESum / static_cast<double>(allowedElements_.size());

  double dE = dElvfe + dEAverage;

  return dE;
}
