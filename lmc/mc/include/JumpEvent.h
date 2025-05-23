/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 6/30/23 4:05 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_MC_INCLUDE_JUMPEVENT_H_
#define LMC_MC_INCLUDE_JUMPEVENT_H_
#include <cstddef>
#include <utility>
#include <array>
namespace mc
{
  class JumpEvent
  {
  public:
    /// Constructor
    JumpEvent();
    JumpEvent(const std::pair<size_t, size_t> &jump_pair,
              const double forward_barrier,
              const double dE,
              double beta);

    /// Getter
    [[nodiscard]] const std::pair<size_t, size_t> &GetIdJumpPair() const;
    [[nodiscard]] double GetForwardBarrier() const;
    [[nodiscard]] double GetForwardRate() const;

    // Zhucongs way using forward barrier and dE from CE
    [[nodiscard]] double GetBackwardBarrier() const;

    // From Model
    // double GetTrueBackwardBarrier() const;

    [[nodiscard]] double GetBackwardRate() const;
    // Returns Energy Change from the CE
    [[nodiscard]] double GetEnergyChange() const;

    // dE from barrier
    // double GetdEBarrier() const;

    [[nodiscard]] double GetProbability() const;
    [[nodiscard]] double GetCumulativeProbability() const;
    [[nodiscard]] JumpEvent GetReverseJumpEvent() const;
    /// Setter
    void SetProbability(double probability);
    void SetCumulativeProbability(double cumulative_probability);
    void CalculateProbability(double total_rates);

  private:
    double beta_{};
    std::pair<size_t, size_t> jump_pair_{};
    double forward_barrier_{};
    // double backward_barrier_{};
    double energy_change_{};
    double forward_rate_{};

    double probability_{};
    double cumulative_probability_{};
  };
} // mc
#endif // LMC_MC_INCLUDE_JUMPEVENT_H_
