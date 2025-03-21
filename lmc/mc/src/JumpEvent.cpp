/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/30/23 3:09 PM                                                          *
 **************************************************************************************************/

// With Implementation of xSv Determinant

#include <cmath>
#include <utility>
#include "JumpEvent.h"

namespace mc
{
  JumpEvent::JumpEvent() = default;
  JumpEvent::JumpEvent(const std::pair<size_t, size_t> &jump_pair,
                       const double forward_barrier,
                       // const double backward_barrier,
                       // Energy change from CE
                       const double dE,
                       double beta)
      : beta_(beta),
        jump_pair_(std::move(jump_pair)),
        forward_barrier_(forward_barrier),
        // From Vacancy Migration Barrier Predictor
        // backward_barrier_(backward_barrier),
        energy_change_(dE),
        forward_rate_(std::exp(-forward_barrier_ * beta_))
  {}

  const std::pair<size_t, size_t> &JumpEvent::GetIdJumpPair() const
  {
    return jump_pair_;
  }
  double JumpEvent::GetForwardBarrier() const
  {
    return forward_barrier_;
  }
  double JumpEvent::GetForwardRate() const
  {
    return forward_rate_;
  }

  // Using Zhucong's way
  double JumpEvent::GetBackwardBarrier() const
  {
    return forward_barrier_ - energy_change_;
  }

  // Backward Barrier from predictor
  // double JumpEvent::GetTrueBackwardBarrier() const
  // {
  //   return backward_barrier_;
  // }

  double JumpEvent::GetBackwardRate() const
  {
    return std::exp(-GetBackwardBarrier() * beta_);
  }

  // From CE
  double JumpEvent::GetEnergyChange() const
  {
    return energy_change_;
  }

  // // Energy Change from barrier
  // double JumpEvent::GetdEBarrier() const
  // {
  //   return forward_barrier_ - backward_barrier_;
  // }

  double JumpEvent::GetProbability() const
  {
    return probability_;
  }

  double JumpEvent::GetCumulativeProbability() const
  {
    return cumulative_probability_;
  }

  void JumpEvent::SetProbability(double probability)
  {
    probability_ = probability;
  }

  void JumpEvent::SetCumulativeProbability(double cumulative_probability)
  {
    cumulative_probability_ = cumulative_probability;
  }

  void JumpEvent::CalculateProbability(double total_rates)
  {
    probability_ = forward_rate_ / total_rates;
  }
  JumpEvent JumpEvent::GetReverseJumpEvent() const
  {

    return JumpEvent{{jump_pair_.second, jump_pair_.first},
                     GetBackwardBarrier(),
                     // forward_barrier_,
                     -energy_change_,
                     beta_};
  }
} // mc
