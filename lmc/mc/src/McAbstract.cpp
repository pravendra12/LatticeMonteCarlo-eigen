/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file McAbstract.h
 * @brief File contains implementation of McAbstract Class.
 */

#include "McAbstract.h"

McAbstract::McAbstract(Config config,
                       unsigned long long int logDumpSteps,
                       unsigned long long int configDumpSteps,
                       unsigned long long int maximumSteps,
                       unsigned long long int thermodynamicAveragingSteps,
                       unsigned long long int restartSteps,
                       double restartEnergy,
                       double restartTime,
                       double temperature,
                       const string &logFilename)
    : config_(move(config)),
      logDumpSteps_(logDumpSteps),
      configDumpSteps_(configDumpSteps),
      maximumSteps_(maximumSteps),
      steps_(restartSteps),
      energy_(restartEnergy),
      absolute_energy_(0), // Initial Energy will be Zero
      time_(restartTime),
      temperature_(temperature),
      beta_(1.0 / constants::kBoltzmann / temperature_),
      is_restarted_(steps_ > 0),
      thermodynamicAveraging_(thermodynamicAveragingSteps),
      generator_(static_cast<unsigned long long int>(
          chrono::system_clock::now().time_since_epoch().count())),
      unitDistribution_(0.0, 1.0)
{
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
  // Only the writer may open/truncate the shared log file.
  if (world_rank_ == 0)
  {
    ofs_.open(logFilename, is_restarted_ ? ofstream::app : ofstream::out);
    if (!ofs_)
    {
      throw runtime_error("Cannot open " + logFilename);
    }
    ofs_.precision(16);
  }
}
McAbstract::~McAbstract()
{
  MPI_Finalize();
}
