#ifndef LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#include <random>
#include <omp.h>
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
#include "VacancyMigrationPredictor.h"

namespace mc {
class KineticMcFirstOmp : public KineticMcFirstAbstract {
  public:
    /**
     * @brief Construct an OpenMP KMC simulation with displacement tracking.
     * @param logDumpMode Event sampling schedule: "adaptive" or "linear".
     * @param speciesDisplacements Initial species totals supplied in the parameter file.
     * @param displacementRestartFilename Optional gzip atom snapshot matching the restart configuration.
     * Other parameters are documented by KineticMcFirstAbstract.
     */
    KineticMcFirstOmp(Config config,
                      const unsigned long long int logDumpSteps,
                      const unsigned long long int configDumpSteps,
                      const unsigned long long int maximumSteps,
                      const unsigned long long int thermodynamicAveragingSteps,
                      const unsigned long long int restartSteps,
                      const double restartEnergy,
                      const double restartTime,
                      const double temperature,
                      VacancyMigrationPredictor &vacancyMigrationPredictor,
                      const string &timeTemperatureFilename,
                      const bool isRateCorrector,
                      const Eigen::RowVector3d &vacancyTrajectory,
                      const string &logDumpMode = "adaptive",
                      const string &displacementRestartFilename = "",
                      const map<Element, Eigen::RowVector3d> &speciesDisplacements = {});
    ~KineticMcFirstOmp() override;
  protected:
    void BuildEventList() override;
    double CalculateTime() override;
};
} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
