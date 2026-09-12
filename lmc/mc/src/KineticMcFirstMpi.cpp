
#include "KineticMcFirstMpi.h"

/**
 * @file KineticMcFirstMpi.cpp
 * @brief Implements First Order Residence Time Algorithm using MPI.
 *
 *            k -> i
 *            |
 *            Current Position
 */

namespace mc
{

  KineticMcFirstMpi::KineticMcFirstMpi(Config config,
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
                                       const string &logDumpMode,
                                       const string &displacementRestartFilename,
                                       const map<Element, Eigen::RowVector3d> &speciesDisplacements)
      : KineticMcFirstAbstract(move(config),
                               logDumpSteps,
                               configDumpSteps,
                               maximumSteps,
                               thermodynamicAveragingSteps,
                               restartSteps,
                               restartEnergy,
                               restartTime,
                               temperature,
                               vacancyMigrationPredictor,
                               timeTemperatureFilename,
                               isRateCorrector,
                               vacancyTrajectory,
                               logDumpMode,
                               displacementRestartFilename,
                               speciesDisplacements)
  {
    if (world_size_ != kEventListSize_)
    {
      cout << "Must use " << kEventListSize_ << " processes. Terminating...\n"
           << endl;
      MPI_Finalize();
      exit(0);
    }
    if (world_rank_ == 0)
    {
      cout << "Using " << world_size_ << " processes." << endl;
    }
  }
  KineticMcFirstMpi::~KineticMcFirstMpi() = default;

  void KineticMcFirstMpi::BuildEventList()
  {
    total_rate_k_ = 0;

    // Neighbours of Vacancy
    const auto neighbor_vacancy_id =
        config_.GetNeighborLatticeIdVectorOfLattice(vacancyLatticeId_, 1)[static_cast<size_t>(world_rank_)];

    JumpEvent local_event_k_i(
        // Jump Pair
        {vacancyLatticeId_, neighbor_vacancy_id},
        vacancyMigrationPredictor_.GetBarrierAndDeltaE(config_,
                                                       {vacancyLatticeId_,
                                                        neighbor_vacancy_id}),
        beta_);

    event_k_i_ = local_event_k_i; // Assign the local variable to the member variable

    const double rate_k_i = event_k_i_.GetForwardRate();

    MPI_Allreduce(&rate_k_i, &total_rate_k_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    event_k_i_.CalculateProbability(total_rate_k_);

    MPI_Allgather(&event_k_i_,
                  sizeof(JumpEvent),
                  MPI_BYTE,
                  event_k_i_list_.data(),
                  sizeof(JumpEvent),
                  MPI_BYTE,
                  MPI_COMM_WORLD);
    double cumulative_probability = 0.0;
    for (auto &event_it : event_k_i_list_)
    {
      cumulative_probability += event_it.GetProbability();
      event_it.SetCumulativeProbability(cumulative_probability);
    }
  }

  double KineticMcFirstMpi::CalculateTime()
  {
    static uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
    return -log(uniform_distribution(generator_)) / total_rate_k_ / constants::kPrefactor;
  }

} // mc