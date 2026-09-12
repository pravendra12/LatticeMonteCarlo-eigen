/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file KineticMcAbstract.cpp
 * @brief File contains implementation of KineticMcAbstract Class.
 */

#include "KineticMcAbstract.h"
#include <filesystem>
#include <limits>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace mc
{

  KineticMcFirstAbstract::KineticMcFirstAbstract(Config config,
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
      : McAbstract(move(config),
                   logDumpSteps,
                   configDumpSteps,
                   maximumSteps,
                   thermodynamicAveragingSteps,
                   restartSteps,
                   restartEnergy,
                   restartTime,
                   temperature,
                   "kmc_log.txt"),
        kEventListSize_(
            config_.GetNeighborLatticeIdVectorOfLattice(0, 1).size()),
        vacancyMigrationPredictor_(
            vacancyMigrationPredictor),
        timeTemperatureInterpolator_(
            timeTemperatureFilename),
        isTimeTemperatureInterpolator_(!timeTemperatureFilename.empty()),
        isRateCorrector_(isRateCorrector),
        vacancyLatticeId_(config_.GetVacancyLatticeId()),
        vacancyTrajectory_(vacancyTrajectory),
        logDumpMode_(logDumpMode),
        atomDisplacements_(config_.GetNumAtoms(), Eigen::RowVector3d::Zero()),
        displacementOriginSteps_(restartSteps),
        displacementOriginTime_(restartTime),
        displacementSegmentStartSteps_(restartSteps),
        event_k_i_list_(kEventListSize_)
  {
    if (logDumpMode_ != "adaptive" && logDumpMode_ != "linear")
    {
      throw invalid_argument("log_dump_mode must be adaptive or linear");
    }
    if (logDumpSteps_ == 0 || configDumpSteps_ == 0)
    {
      throw invalid_argument("KMC log_dump_steps and config_dump_steps must be positive");
    }
    if (restartSteps > maximumSteps_)
    {
      throw invalid_argument("restart_steps must not exceed maximum_steps");
    }
    ofs_.precision(numeric_limits<double>::max_digits10);
    for (const auto &element : config_.GetAtomVector())
    {
      if (element != ElementName::X)
      {
        speciesDisplacements_.try_emplace(element, Eigen::RowVector3d::Zero());
        speciesSquaredDisplacements_.try_emplace(element, 0.0);
        // Species belongs to the persistent atom, not its current lattice site.
        ++speciesAtomCounts_[element];
      }
    }
    for (const auto &[element, displacement] : speciesDisplacements)
    {
      if (!speciesDisplacements_.count(element) || !displacement.allFinite())
      {
        throw invalid_argument("species_displacement must name a species in the configuration and contain finite values");
      }
      speciesDisplacements_.at(element) = displacement;
    }
    if (!displacementRestartFilename.empty())
    {
      ReadAtomDisplacements(displacementRestartFilename);
    }
    if (world_rank_ == 0 && is_restarted_ && displacementRestartFilename.empty())
    {
      cout << "No atom snapshot supplied; starting a new tracer measurement at step "
           << steps_ << "." << endl;
    }
  }

  KineticMcFirstAbstract::~KineticMcFirstAbstract() = default;

  void KineticMcFirstAbstract::UpdateTemperature()
  {
    if (isTimeTemperatureInterpolator_)
    {
      temperature_ = timeTemperatureInterpolator_.GetTemperature(time_);
      beta_ = 1.0 / constants::kBoltzmann / temperature_;
    }
  }

  //
  // double KineticMcFirstAbstract::GetTimeCorrectionFactor() {
  //   if (isRateCorrector_) {
  //     return rate_corrector_.GetTimeCorrectionFactor(temperature_);
  //   }
  //   return 1.0;
  // }

  void KineticMcFirstAbstract::Dump() const
  {
    if (world_rank_ != 0)
    {
      return;
    }
    const bool firstSnapshot = !firstSnapshotWritten_;
    if (firstSnapshot)
    {
      ostringstream header;
      header << "steps\ttime\taverage_time\ttemperature\tenergy\taverage_energy\tEa\tdE\tEa_Backward\tselected\tvacancy_trajectory";
      for (const auto &[element, displacement] : speciesDisplacements_)
      {
        header << "\tdR_" << element;
      }
      for (const auto &[element, displacement] : speciesDisplacements_)
      {
        header << "\tMSD_" << element;
      }
      // Never append the new scalar fields under an older or mismatched header.
      ofs_.seekp(0, ios::end);
      if (ofs_.tellp() == streampos(0))
      {
        ofs_ << header.str() << endl;
      }
      else
      {
        ifstream existingLog("kmc_log.txt");
        string existingHeader;
        if (!getline(existingLog, existingHeader) || existingHeader != header.str())
          throw runtime_error("Incompatible kmc_log.txt header; restart in a new directory or archive the old log");
      }
    }
    if (steps_ % configDumpSteps_ == 0)
    {
      config_.WriteConfig(to_string(steps_) + ".cfg.gz", config_);
    }
    if (steps_ == maximumSteps_)
    {
      config_.WriteConfig("end.cfg", config_);
    }

    // Always include the measurement's initial and final states.
    // Configuration steps also need log counters and an atom snapshot for restart.
    if (firstSnapshot || steps_ == maximumSteps_ || steps_ % configDumpSteps_ == 0 ||
        steps_ % GetLogDumpSteps() == 0)
    {
      ofs_ << steps_ << '\t'
           << time_ << '\t'
           << 1 / (total_rate_k_ * constants::kPrefactor) << "\t"
           << temperature_ << '\t'
           << energy_ << '\t'
           << thermodynamicAveraging_.GetThermodynamicAverage(beta_) << "\t"
           << event_k_i_.GetForwardBarrier() << '\t'
           << event_k_i_.GetEnergyChange() << '\t'
           << event_k_i_.GetBackwardBarrier() << "\t"
           << event_k_i_.GetIdJumpPair().second << '\t'
           << vacancyTrajectory_;
      for (const auto &[element, displacement] : speciesDisplacements_)
      {
        ofs_ << '\t' << displacement;
      }
      // The running sums already describe this state; logging needs no atom scan.
      for (const auto &[element, displacement] : speciesDisplacements_)
      {
        const auto count = speciesAtomCounts_.at(element);
        ofs_ << '\t' << (count == 0 ? 0.0 :
            speciesSquaredDisplacements_.at(element) / static_cast<double>(count));
      }
      ofs_ << endl;
      if (!ofs_)
      {
        throw runtime_error("Failed writing kmc_log.txt");
      }
      if (steps_ % configDumpSteps_ == 0 || steps_ == maximumSteps_)
      {
        DumpAtomDisplacements();
      }
      firstSnapshotWritten_ = true;
    }
  }

  unsigned long long int KineticMcFirstAbstract::GetLogDumpSteps() const
  {
    if (logDumpMode_ == "linear")
    {
      return logDumpSteps_;
    }
    // Integer arithmetic avoids negative floating-point exponents at step zero
    // and overflow in the original 10 * logDumpSteps_ threshold.
    if (steps_ > 0 && (steps_ - 1) / 10 >= logDumpSteps_)
    {
      return logDumpSteps_;
    }
    auto decade = steps_ == numeric_limits<unsigned long long int>::max()
                      ? steps_ / 10 : (steps_ + 1) / 10;
    unsigned long long int interval = 1;
    while (decade >= 10 && interval < logDumpSteps_)
    {
      decade /= 10;
      interval *= 10;
    }
    return min(interval, logDumpSteps_);
  }

  void KineticMcFirstAbstract::UpdateDisplacements()
  {
    const auto destination = event_k_i_.GetIdJumpPair().second;
    const auto atomId = config_.GetAtomIdOfLattice(destination);
    const auto element = config_.GetElementOfAtom(atomId);
    const Eigen::RowVector3d vacancyJump =
        config_.GetRelativeDistanceVectorLattice(vacancyLatticeId_, destination)
            .transpose() * config_.GetBasis();

    // The exchanging atom moves opposite to the vacancy. Accumulating jump
    // vectors preserves unwrapped displacement across periodic boundaries.
    auto &displacement = atomDisplacements_.at(atomId);
    auto &squaredDisplacements = speciesSquaredDisplacements_.at(element);
    // Replace only this atom's contribution to sum_i |dR_i|^2. Squaring each
    // jump instead would lose correlations, including cancellation on return.
    squaredDisplacements -= displacement.squaredNorm();
    displacement -= vacancyJump;
    squaredDisplacements += displacement.squaredNorm();
    speciesDisplacements_.at(element) -= vacancyJump;
    vacancyTrajectory_ += vacancyJump;
  }

  void KineticMcFirstAbstract::DumpAtomDisplacements() const
  {
    if (world_rank_ != 0)
    {
      return;
    }
    const string filename = to_string(steps_) + ".displacements.gz";
    const string temporaryFilename = filename + ".tmp";
    try
    {
      ofstream output;
      output.exceptions(ios::badbit | ios::failbit);
      output.open(temporaryFilename, ios::binary | ios::trunc);
      boost::iostreams::filtering_ostream compressed;
      compressed.push(boost::iostreams::gzip_compressor());
      compressed.push(output);
      compressed.exceptions(ios::badbit | ios::failbit);
      compressed.precision(numeric_limits<double>::max_digits10);
      compressed << "# step " << steps_
                 << "\n# time " << time_
                 << "\n# displacement_origin_step " << displacementOriginSteps_
                 << "\n# displacement_origin_time " << displacementOriginTime_
                 << "\n# segment_start_step " << displacementSegmentStartSteps_
                 << "\natom_id element dx dy dz\n";

      // Stream rows directly into gzip instead of buffering the whole snapshot.
      const Eigen::IOFormat vectorFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n");
      for (size_t atomId = 0; atomId < atomDisplacements_.size(); ++atomId)
      {
        const auto element = config_.GetElementOfAtom(atomId);
        if (element != ElementName::X)
        {
          compressed << atomId << ' ' << element << ' '
                     << atomDisplacements_[atomId].format(vectorFormat) << '\n';
        }
      }
      // Finish the gzip trailer before publishing, so a restarted run replaces
      // the old snapshot with one complete file rather than appending more rows.
      compressed.flush();
      // Detaching the filter buffer sets badbit even on successful completion.
      // The underlying file stream still throws on write/close errors.
      compressed.exceptions(ios::goodbit);
      compressed.reset();
      output.close();
      std::filesystem::rename(temporaryFilename, filename);
    }
    catch (const exception &error)
    {
      std::error_code ignored;
      std::filesystem::remove(temporaryFilename, ignored);
      throw runtime_error("Failed writing " + filename + ": " + error.what());
    }
  }

  void KineticMcFirstAbstract::ReadAtomDisplacements(const string &filename)
  {
    try
    {
      ifstream file(filename, ios::binary);
      if (!file)
        throw runtime_error("Cannot open snapshot");
      boost::iostreams::filtering_istream input;
      input.push(boost::iostreams::gzip_decompressor());
      input.push(file);
      input.exceptions(ios::badbit);

      // Read metadata written by DumpAtomDisplacements; no positions are needed.
      map<string, string> metadata;
      string line;
      while (getline(input, line) && !line.empty() && line.front() == '#')
      {
        istringstream header(line.substr(1));
        string key, value;
        if (!(header >> key >> value) || !metadata.emplace(key, value).second)
          throw runtime_error("Invalid snapshot metadata");
      }
      if (line != "atom_id element dx dy dz")
        throw runtime_error("Invalid atom displacement columns");
      const auto snapshotStep = stoull(metadata.at("step"));
      const double snapshotTime = stod(metadata.at("time"));
      const auto originStep = stoull(metadata.at("displacement_origin_step"));
      const double originTime = stod(metadata.at("displacement_origin_time"));
      if (snapshotStep != steps_ || !isfinite(snapshotTime) || !isfinite(time_) ||
          abs(snapshotTime - time_) > 1e-12 * max(abs(snapshotTime), abs(time_)) ||
          originStep > steps_ || !isfinite(originTime) || originTime > snapshotTime)
        throw runtime_error("Snapshot step, time or measurement origin does not match restart");

      // IDs follow CFG atom ordering, with the vacancy omitted from the snapshot.
      for (size_t atomId = 0; atomId < atomDisplacements_.size(); ++atomId)
      {
        const auto element = config_.GetElementOfAtom(atomId);
        if (element == ElementName::X)
          continue;
        size_t savedId;
        string savedElement;
        Eigen::RowVector3d displacement;
        if (!(input >> savedId >> savedElement >> displacement(0) >> displacement(1) >> displacement(2)) ||
            savedId != atomId || Element(savedElement) != element || !displacement.allFinite())
          throw runtime_error("Snapshot atom IDs, species or displacement vectors do not match configuration");
        atomDisplacements_[atomId] = displacement;
      }
      string extra;
      if (input >> extra)
        throw runtime_error("Unexpected extra atom snapshot data");
      // Rebuild once from the restored vectors at their original tracer origin.
      // A logged MSD scalar cannot replace these vectors for future hop updates.
      for (auto &[element, squaredDisplacements] : speciesSquaredDisplacements_)
        squaredDisplacements = 0.0;
      for (size_t atomId = 0; atomId < atomDisplacements_.size(); ++atomId)
      {
        const auto element = config_.GetElementOfAtom(atomId);
        if (element != ElementName::X)
          speciesSquaredDisplacements_.at(element) += atomDisplacements_[atomId].squaredNorm();
      }
      displacementOriginSteps_ = originStep;
      displacementOriginTime_ = originTime;
      time_ = snapshotTime;
    }
    catch (const exception &error)
    {
      throw runtime_error("Cannot restore atom displacements from " + filename + ": " + error.what());
    }
  }

  size_t KineticMcFirstAbstract::SelectEvent() const
  {
    const double random_number = unitDistribution_(generator_);
    auto it = lower_bound(
        event_k_i_list_.begin(), event_k_i_list_.end(), random_number, [](const auto &lhs, double value)
        { return lhs.GetCumulativeProbability() < value; });

    // If not find (maybe generated 1), which rarely happens, returns the last event
    if (it == event_k_i_list_.cend())
    {
      it--;
    }
    if (world_size_ > 1)
    {
      int event_id = static_cast<int>(it - event_k_i_list_.cbegin());
      MPI_Bcast(&event_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return static_cast<size_t>(event_id);
    }
    else
    {
      return static_cast<size_t>(it - event_k_i_list_.cbegin());
    }
  }

  void KineticMcFirstAbstract::Debug(double one_step_time) const
  {
    if (world_rank_ == 0)
    {
      if (isnan(one_step_time) or isinf(one_step_time) or one_step_time < 0.0)
      {

        config_.WriteConfig("debug" + to_string(steps_) + ".cfg.gz", config_);
        cerr << "Invalid time step: " << one_step_time << endl;
        cerr << "For each event: Energy Barrier, Energy Change, Probability, " << endl;
        for (auto &event : event_k_i_list_)
        {
          cerr << event.GetForwardBarrier() << '\t' << event.GetEnergyChange() << '\t'
               << event.GetCumulativeProbability() << endl;
        }
        throw runtime_error("Invalid time step");
      }
    }
  }

  void KineticMcFirstAbstract::OneStepSimulation()
  {
    UpdateTemperature();

    thermodynamicAveraging_.AddEnergy(energy_);

    BuildEventList();

    // double one_step_time = CalculateTime() * GetTimeCorrectionFactor();
    double one_step_time = CalculateTime();
    Debug(one_step_time);

    event_k_i_ = event_k_i_list_[SelectEvent()];

    Dump();

    // The final snapshot is the terminal state, not the start of another jump.
    if (steps_ == maximumSteps_)
    {
      return;
    }

    // modify
    time_ += one_step_time;
    energy_ += event_k_i_.GetEnergyChange();
    absolute_energy_ += event_k_i_.GetEnergyChange();

    UpdateDisplacements();
    config_.LatticeJump(event_k_i_.GetIdJumpPair());
    ++steps_;
    vacancyLatticeId_ = event_k_i_.GetIdJumpPair().second;
  }

  void KineticMcFirstAbstract::Simulate()
  {
    while (steps_ < maximumSteps_)
    {
      OneStepSimulation();
    }
    OneStepSimulation(); // Write the terminal state without executing another jump.
  }

  KineticMcChainAbstract::KineticMcChainAbstract(Config config,
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
                               speciesDisplacements),
        previous_j_lattice_id_(config_.GetNeighborLatticeIdVectorOfLattice(vacancyLatticeId_, 1)[0]),
        l_lattice_id_list_(kEventListSize_)
  {
    MPI_Op_create(DataSum, 1, &mpi_op_);
    DefineStruct(&mpi_datatype_);
  }

  void KineticMcChainAbstract::OneStepSimulation()
  {
    KineticMcFirstAbstract::OneStepSimulation();
    previous_j_lattice_id_ = event_k_i_.GetIdJumpPair().first;
  }

  KineticMcChainAbstract::~KineticMcChainAbstract()
  {
    MPI_Op_free(&mpi_op_);
    MPI_Type_free(&mpi_datatype_);
  }

} // namespace mc
