#include "KineticMcFirstOmp.h"

namespace mc
{

  KineticMcFirstOmp::KineticMcFirstOmp(Config config,
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
    if (world_size_ != 1)
    {
      std::cout << "Must use 1 process. Terminating...\n"
                << std::endl;
      MPI_Finalize();
      exit(0);
    }
#pragma omp parallel default(none) shared(std::cout)
    {
#pragma omp master
      {
        std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
      }
    }
  }
  KineticMcFirstOmp::~KineticMcFirstOmp() = default;

  void KineticMcFirstOmp::BuildEventList()
  {
    total_rate_k_ = 0.0;
    const auto neighbor_lattice_id_vector = config_.GetNeighborLatticeIdVectorOfLattice(vacancyLatticeId_, 1);
#pragma omp parallel default(none) shared(neighbor_lattice_id_vector) reduction(+ : total_rate_k_)
    {
#pragma omp for
      for (size_t i = 0; i < kEventListSize_; ++i)
      {
        const auto neighbor_lattice_id = neighbor_lattice_id_vector[i];
        JumpEvent lattice_jump_event(
            {vacancyLatticeId_, neighbor_lattice_id},
            vacancyMigrationPredictor_.GetBarrierAndDeltaE(config_,
                                                           {vacancyLatticeId_,
                                                            neighbor_lattice_id}),
            beta_);
        total_rate_k_ += lattice_jump_event.GetForwardRate();
        event_k_i_list_[i] = std::move(lattice_jump_event);
      }
    }
    double cumulative_probability = 0.0;
    for (auto &event_it : event_k_i_list_)
    {
      event_it.CalculateProbability(total_rate_k_);
      cumulative_probability += event_it.GetProbability();
      event_it.SetCumulativeProbability(cumulative_probability);
    }
  }
  double KineticMcFirstOmp::CalculateTime()
  {
    static std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
    return -std::log(uniform_distribution(generator_)) / total_rate_k_ / constants::kPrefactor;
  }
} // mc