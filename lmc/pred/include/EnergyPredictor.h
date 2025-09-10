#ifndef LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

#include "Config.h"
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
#include "ClusterExpansion.h"
#include "SymmetrySpglib.h"
#include "JsonUtility.h"
#include <string>
#include "ClusterExpansionParameters.h"
#include "CorrelationVector.h"
#include "BasisSet.h"
=======
#include "SymmetricCE.h"
#include "PrintUtility.h"
#include "ClusterExpansionParameters.h"
>>>>>>> Stashed changes
=======
#include "SymmetricCE.h"
#include "PrintUtility.h"
#include "ClusterExpansionParameters.h"
>>>>>>> Stashed changes
=======
#include "SymmetricCE.h"
#include "PrintUtility.h"
#include "ClusterExpansionParameters.h"
>>>>>>> Stashed changes

using namespace std;

class EnergyPredictor
{
public:
<<<<<<< Updated upstream
<<<<<<< Updated upstream
<<<<<<< Updated upstream
  EnergyPredictor(
    const ClusterExpansionParameters &ceParams,
      const Config &config);

  // Returns the energy of a site
  double ComputeEnergyOfSite(
      const Config &config,
      const size_t &latticeId);
=======
    EnergyPredictor(
        const ClusterExpansionParameters &ceParams,
        const Config &supercellConfig,
        const Config &primitiveConfig);

    double ComputeLocalFormationEnergyOfSite(
        const Config &config,
        const size_t &latticeId);
>>>>>>> Stashed changes

    // Returns total energy of config using Symmetric site centered cluster expansion
    // E = sum_i * Ei
    double ComputeEnergyOfConfig(const Config &config);

<<<<<<< Updated upstream
  // Returns the dE due to vacancy migration
  double GetDeMigration(
    const Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair);

  // Returns the dE due to atom swap 
  double GetDeSwap(
    Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair);
=======
    // Returns the dE due to atom swap
    double GetDeSwap(
        Config &config,
        const pair<size_t, size_t> &latticeIdJumpPair);

>>>>>>> Stashed changes

  private:
    // Returns total formation energy
    double GetTotalFormationEnergy(
        const vector<double> &clusterVector);

<<<<<<< Updated upstream
private:
  /*! @brief Maximum Cluster Size
   */
  const size_t maxClusterSize_{};

  /*! @brief Maximum Bond Order
   */
  const size_t maxBondOrder_{};

  /*! @brief Effective cluster interactions
   */
  const VectorXd ecis_{};

  /*! @brief Basis Set
  */
  BasisSet atomicBasis_;

  /*! @brief Element set
   */
  const set<Element> elementSet_{};

  /*! @brief Equivalent Clusters Encoding
   */
  const vector<pair<vector<vector<size_t>>,
                    LatticeClusterType>>
      equivalentClustersEncoding_;
=======
    // Returns total energy
    double GetTotalEnergy(
        const vector<double> &clusterVector);

=======
    EnergyPredictor(
        const ClusterExpansionParameters &ceParams,
        const Config &supercellConfig,
        const Config &primitiveConfig);

    double ComputeLocalFormationEnergyOfSite(
        const Config &config,
        const size_t &latticeId);

    // Returns total energy of config using Symmetric site centered cluster expansion
    // E = sum_i * Ei
    double ComputeEnergyOfConfig(const Config &config);

    // Returns the dE due to atom swap
    double GetDeSwap(
        Config &config,
        const pair<size_t, size_t> &latticeIdJumpPair);


  private:
    // Returns total formation energy
    double GetTotalFormationEnergy(
        const vector<double> &clusterVector);

    // Returns total energy
    double GetTotalEnergy(
        const vector<double> &clusterVector);

>>>>>>> Stashed changes
=======
    EnergyPredictor(
        const ClusterExpansionParameters &ceParams,
        const Config &supercellConfig,
        const Config &primitiveConfig);

    double ComputeLocalFormationEnergyOfSite(
        const Config &config,
        const size_t &latticeId);

    // Returns total energy of config using Symmetric site centered cluster expansion
    // E = sum_i * Ei
    double ComputeEnergyOfConfig(const Config &config);

    // Returns the dE due to atom swap
    double GetDeSwap(
        Config &config,
        const pair<size_t, size_t> &latticeIdJumpPair);


  private:
    // Returns total formation energy
    double GetTotalFormationEnergy(
        const vector<double> &clusterVector);

    // Returns total energy
    double GetTotalEnergy(
        const vector<double> &clusterVector);

>>>>>>> Stashed changes
    static unordered_map<string, int> GetElementCountMap(
        const Config &supercellConfig);
    
    const size_t nAtoms_{};
    const vector<string> allowedElements_{};
    const vector<double> clusterCutoffs_{};
    SymmetricCE symCE_;
    const vector<vector<vector<int>>> localOrbitsEncoding_;
    const unordered_map<string, double> chemicalPotentialsMap_;
    const unordered_map<string, int> elementCountMap_{};
    const vector<double> ecis_{};
<<<<<<< Updated upstream
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
};

#endif // LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

