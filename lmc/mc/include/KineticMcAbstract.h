/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file KineticMcAbstract.h
 * @brief File contains declaration of KineticMc abstract class.
 *
 */

#ifndef LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_

#include <random>
#include <map>
#include <omp.h>
#include <mpi.h>
#include <Eigen/Dense>
#include "McAbstract.h"
#include "JumpEvent.h"
#include "TimeTemperatureInterpolator.h"
#include "VacancyMigrationPredictor.h"

using namespace std;

namespace mc
{

  /**
   * @brief Abstract class for Kinetic Monte Carlo Simulation.
   *
   *  Based on First-Order Residence Time Algorithm.
   */
  class KineticMcFirstAbstract : public McAbstract
  {
  public:
    /**
     * Constructor for KineticMcFirstAbstract.
     *
     * Initializes the kinetic Monte Carlo simulation with the given parameters.
     *
     * @param config Simulation configuration.
     * @param supercellConfig Training configuration used for the Cluster Expansion Model.
     * @param logDumpSteps Steps between logging progress.
     * @param configDumpSteps Steps between configuration dumps.
     * @param maximumSteps Maximum simulation steps.
     * @param thermodynamicAveragingSteps Steps for thermodynamic averaging.
     * @param restartSteps Steps for restarting the simulation.
     * @param restartEnergy Restart energy.
     * @param restartTime Restart time.
     * @param temperature Simulation temperature (in Kelvin).
     * @param elementSet Set of elements involved in the simulation.
     * @param predictorFilename Path to JSON file with cluster interaction coefficients.
     * @param timeTemperatureFilename Path to time-temperature data file.
     * @param isRateCorrector Whether rate correction needs to be applied.
     * @param vacancyTrajectory Initial vacancy trajectory vector.
     * @param logDumpMode Event sampling schedule: "adaptive" or "linear".
     * @param speciesDisplacements Initial species totals supplied in the parameter file.
     * @param displacementRestartFilename Optional gzip atom snapshot matching config and restart time.
     */
    KineticMcFirstAbstract(Config config,
                           unsigned long long int logDumpSteps,
                           unsigned long long int configDumpSteps,
                           unsigned long long int maximumSteps,
                           unsigned long long int thermodynamicAveragingSteps,
                           unsigned long long int restartSteps,
                           double restartEnergy,
                           double restartTime,
                           double temperature,
                           VacancyMigrationPredictor &vacancyMigrationPredictor,
                           const string &timeTemperatureFilename,
                           bool isRateCorrector,
                           const Eigen::RowVector3d &vacancyTrajectory,
                           const string &logDumpMode = "adaptive",
                           const string &displacementRestartFilename = "",
                           const map<Element, Eigen::RowVector3d> &speciesDisplacements = {});

    /**
     * @brief Destructor for KineticMcFirstAbstract.
     */
    ~KineticMcFirstAbstract() override;

    /**
     * @brief Deleted copy constructor.
     */
    KineticMcFirstAbstract(const KineticMcFirstAbstract &) = delete;

    /**
     * @brief Deleted assignment operator.
     */
    void operator=(const mc::KineticMcFirstAbstract &) = delete;

    /**
     * @brief Starts the Kinetic Monte Carlo Simulation.
     */
    void Simulate() override;

  protected:
    /**
     * @brief Update the temperature based on the current temperature.
     */
    void UpdateTemperature();

    /**
     * @brief Dumps the current simulation state.
     */
    virtual void Dump() const;

    /**
     * @brief Accumulate unwrapped Cartesian displacements for the selected exchange.
     * Updates the vacancy, moving atom, and its species on every accepted event.
     * Must be called before LatticeJump changes the atom-to-lattice mapping.
     */
    void UpdateDisplacements();

    /**
     * @brief Return the current interval between log and atom trajectory snapshots.
     * @return logDumpSteps_ in linear mode; a bounded power of ten during the
     * early part of adaptive mode, followed by logDumpSteps_.
     */
    [[nodiscard]] unsigned long long int GetLogDumpSteps() const;

    /**
     * @brief Stream one snapshot to <step>.displacements.gz on rank zero.
     * Writes step/time and measurement-origin metadata, followed by whitespace-
     * separated atom_id, element, dx, dy, dz columns in increasing atom-ID order.
     * Vectors are cumulative and unwrapped; vacancies are omitted. Dump() supplies
     * the sampling schedule. A complete gzip file atomically replaces any snapshot
     * at the same step, including the initial state of a restarted simulation.
     * @throws std::runtime_error If compression, writing or replacement fails.
     */
    void DumpAtomDisplacements() const;

    /**
     * @brief Restore per-atom tracer vectors and their origin from a gzip snapshot.
     * Species totals and the vacancy vector remain those supplied in the parameters.
     * @param filename Snapshot matching the restart step, time and atom-ID ordering.
     * @throws std::runtime_error For invalid metadata, atom IDs/types or vectors.
     */
    void ReadAtomDisplacements(const string &filename);

    /**
     * @brief Selects an event to simulate based on rates.
     * @return The lattice Id which will jump to vacant site.
     */
    size_t SelectEvent() const;

    /** @brief Debugging utility to log one-step simulation time.
     *  @param one_step_time Time taken for a single simulation step.
     */
    void Debug(double one_step_time) const;

    /** @brief Builds the event list for possible transitions.
     *
     *  This function must be implemented by derived classes.
     */
    virtual void BuildEventList() = 0;

    /** @brief Calculates the time increment for the simulation step.
     *         This function must be implemented by derived classes.
     *
     *  @return Time increment.
     */
    virtual double CalculateTime() = 0;

    /**
     * @brief Performs one step of the simulation.
     */
    virtual void OneStepSimulation();

    /**
     * @brief  Dynamic size of the event list.
     */
    size_t kEventListSize_;

    // Helpful properties

    /**
     * @brief Vacancy Migration Energy Predictor
     */
    VacancyMigrationPredictor &vacancyMigrationPredictor_;

    /**
     * @brief Time Temperature Interpolator
     *
     */
    const TimeTemperatureInterpolator timeTemperatureInterpolator_;

    /**
     * @brief Indicates if time-temperature interpolation is used
     *
     */
    const bool isTimeTemperatureInterpolator_;

    // @brief  Rate Corrector
    // const pred::RateCorrector rate_corrector_;

    /**
     * @brief Indicates if rate correction is used.
     *
     */
    const bool isRateCorrector_;

    /**
     * @brief Vacancy Lattice Id.
     *
     */
    size_t vacancyLatticeId_;

    /**
     * @brief Vacancy Trajectory Vector
     *
     */
    Eigen::RowVector3d vacancyTrajectory_;

    /// Event-count sampling mode, validated at construction.
    const string logDumpMode_;

    /// Species totals in deterministic element order; vacancy is excluded.
    map<Element, Eigen::RowVector3d> speciesDisplacements_{};

    /// Cumulative unwrapped displacement indexed by persistent atom ID.
    vector<Eigen::RowVector3d> atomDisplacements_{};

    /// Origin of the per-atom tracer measurement; species totals are independent.
    unsigned long long int displacementOriginSteps_;
    double displacementOriginTime_;

    /// Start of this invocation, independent of a restored measurement origin.
    const unsigned long long int displacementSegmentStartSteps_;

    /// Whether this invocation has written its initial/restart snapshot.
    mutable bool firstSnapshotWritten_{false};

    /**
     * @brief Jump events vector.
     *
     */
    vector<JumpEvent> event_k_i_list_{};

    /**
     * @brief Selected jump event.
     *
     */
    JumpEvent event_k_i_{};

    /**
     * @brief Total Rate of the events.
     *
     */
    double total_rate_k_{0.0};
  };

  /**
   * @brief Abstract class for Kinetic Monte Carlo Simulation.
   *
   *  Based on Second-Order Residence Time Algorithm.
   */
  class KineticMcChainAbstract : public KineticMcFirstAbstract
  {
  public:
    /**
     * Constructor for KineticMcFirstAbstract.
     *
     * Initializes the kinetic Monte Carlo simulation with the given parameters.
     *
     * @param config Simulation configuration.
     * @param supercellConfig Training configuration used for the Cluster Expansion Model.
     * @param logDumpSteps Steps between logging progress.
     * @param configDumpSteps Steps between configuration dumps.
     * @param maximumSteps Maximum simulation steps.
     * @param thermodynamicAveragingSteps Steps for thermodynamic averaging.
     * @param restartSteps Steps for restarting the simulation.
     * @param restartEnergy Restart energy.
     * @param restartTime Restart time.
     * @param temperature Simulation temperature (in Kelvin).
     * @param elementSet Set of elements involved in the simulation.
     * @param predictorFilename Path to JSON file with cluster interaction coefficients.
     * @param timeTemperatureFilename Path to time-temperature data file.
     * @param isRateCorrector Whether rate correction needs to be applied.
     * @param vacancyTrajectory Initial vacancy trajectory vector.
     * @param logDumpMode Event sampling schedule: "adaptive" or "linear".
     * @param speciesDisplacements Initial species totals supplied in the parameter file.
     * @param displacementRestartFilename Optional gzip atom snapshot matching config and restart time.
     */
    KineticMcChainAbstract(Config config,
                           unsigned long long int logDumpSteps,
                           unsigned long long int configDumpSteps,
                           unsigned long long int maximumSteps,
                           unsigned long long int thermodynamicAveragingSteps,
                           unsigned long long int restartSteps,
                           double restartEnergy,
                           double restartTime,
                           double temperature,
                           VacancyMigrationPredictor &vacancyMigrationPredictor,
                           const string &timeTemperatureFilename,
                           bool isRateCorrector,
                           const Eigen::RowVector3d &vacancyTrajectory,
                           const string &logDumpMode = "adaptive",
                           const string &displacementRestartFilename = "",
                           const map<Element, Eigen::RowVector3d> &speciesDisplacements = {});

    /**
     * @brief Destructor for KineticMcChainAbstract.
     */
    ~KineticMcChainAbstract() override;

    /**
     * @brief Deleted copy constructor.
     */
    KineticMcChainAbstract(const KineticMcChainAbstract &) = delete;

    /**
     * @brief Deleted assignment operator.
     */
    void operator=(const mc::KineticMcChainAbstract &) = delete;

  protected:
    /**
     * @brief Performs one step of the simulation.
     */
    void OneStepSimulation() override;

    // helpful properties

    /**
     * @brief Lattice ID of the previous jump site.
     *
     */
    size_t previous_j_lattice_id_;

    /**
     * @brief Total rate for jump to nearest neighbours.
     *
     */
    double total_rate_i_{0.0};

    /**
     * @brief Neighbours of choosen lattice Id for the jump.
     *        j -> k -> i ->l
     *        k : current position
     *
     */
    vector<size_t> l_lattice_id_list_{};

    MPI_Op mpi_op_{};
    MPI_Datatype mpi_datatype_{};
  };

  /**
   * @struct MpiData
   * @brief Structure for storing data used in MPI custom reduction.
   */
  struct MpiData
  {
    double beta_bar_k{0.0};
    double beta_k{0.0};
    double gamma_bar_k_j{0.0};
    double gamma_k_j{0.0};
    double beta_k_j{0.0};
    double alpha_k_j{0.0};
    double ts_numerator{0.0};
    double ts_j_numerator{0.0};
  };

  /**
   * @brief Performs a custom reduction operation for MpiData.
   *
   * @param input_buffer Input buffer.
   * @param output_buffer Output buffer.
   * @param len Number of elements.
   * @param datatype MPI datatype.
   */
  inline void DataSum(void *input_buffer,
                      void *output_buffer,
                      int *len,
                      [[maybe_unused]] MPI_Datatype *datatype)
  {
    auto *input = static_cast<MpiData *>(input_buffer);
    auto *output = static_cast<MpiData *>(output_buffer);
    for (int i = 0; i < *len; ++i)
    {
      output[i].beta_bar_k += input[i].beta_bar_k;
      output[i].beta_k += input[i].beta_k;
      output[i].gamma_bar_k_j += input[i].gamma_bar_k_j;
      output[i].gamma_k_j += input[i].gamma_k_j;
      output[i].beta_k_j += input[i].beta_k_j;
      output[i].alpha_k_j += input[i].alpha_k_j;
      output[i].ts_numerator += input[i].ts_numerator;
      output[i].ts_j_numerator += input[i].ts_j_numerator;
    }
  }

  inline void DefineStruct(MPI_Datatype *datatype)
  {
    const int count = 8;
    int block_lens[count];
    MPI_Datatype types[count];
    MPI_Aint displacements[count];

    for (int i = 0; i < count; i++)
    {
      types[i] = MPI_DOUBLE;
      block_lens[i] = 1;
    }
    displacements[0] = offsetof(MpiData, beta_bar_k);
    displacements[1] = offsetof(MpiData, beta_k);
    displacements[2] = offsetof(MpiData, gamma_bar_k_j);
    displacements[3] = offsetof(MpiData, gamma_k_j);
    displacements[4] = offsetof(MpiData, beta_k_j);
    displacements[5] = offsetof(MpiData, alpha_k_j);
    displacements[6] = offsetof(MpiData, ts_numerator);
    displacements[7] = offsetof(MpiData, ts_j_numerator);
    MPI_Type_create_struct(count, block_lens, displacements, types, datatype);
    MPI_Type_commit(datatype);
  }
} // mc

#endif // LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
