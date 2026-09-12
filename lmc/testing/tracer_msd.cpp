/** @file tracer_msd.cpp
 * @brief Exercise production transport updates with prescribed vacancy exchanges.
 * Run through script/test_tracer_msd.py; each process owns one MPI lifetime.
 */
#include "KineticMcAbstract.h"
#include <cassert>
#include <iomanip>

class TracerTest : public mc::KineticMcFirstAbstract
{
public:
  using KineticMcFirstAbstract::KineticMcFirstAbstract;
  using KineticMcFirstAbstract::ReadAtomDisplacements;

  void BuildEventList() override
  {
    // Repeated exchanges with one atom give exactly alternating MSD values.
    event_k_i_list_.assign(1, mc::JumpEvent(
        {vacancyLatticeId_, config_.GetLatticeIdOfAtom(1)}, {0.1, 0.0}, beta_));
    event_k_i_list_[0].SetCumulativeProbability(1.0);
    total_rate_k_ = event_k_i_list_[0].GetForwardRate();
  }
  double CalculateTime() override { return 0.25; }

  void Hop(size_t atom)
  {
    const auto destination = config_.GetLatticeIdOfAtom(atom);
    event_k_i_ = mc::JumpEvent({vacancyLatticeId_, destination}, {0.1, 0.0}, beta_);
    UpdateDisplacements();
    config_.LatticeJump(event_k_i_.GetIdJumpPair());
    vacancyLatticeId_ = destination;
    ++steps_;
    time_ += 0.25;
    Check();
  }

  // Independent brute-force reference, also checking the collective invariant.
  void Check() const
  {
    assert(speciesAtomCounts_.size() == 2);
    assert(speciesAtomCounts_.at(Element("Fe")) == 2);
    assert(speciesAtomCounts_.at(Element("Ni")) == 1);
    Eigen::RowVector3d total = Eigen::RowVector3d::Zero();
    for (const auto &[element, sum] : speciesSquaredDisplacements_)
    {
      double expected = 0.0;
      Eigen::RowVector3d collective = Eigen::RowVector3d::Zero();
      for (size_t i = 0; i < atomDisplacements_.size(); ++i)
        if (config_.GetElementOfAtom(i) == element)
        {
          expected += atomDisplacements_[i].squaredNorm();
          collective += atomDisplacements_[i];
        }
      assert(abs(sum - expected) < 1e-10 * max(1.0, expected));
      assert((collective - speciesDisplacements_.at(element)).norm() < 1e-10);
      total += collective;
    }
    assert((total + vacancyTrajectory_).norm() < 1e-10);
  }

  void Expect(double fe, double ni) const
  {
    assert(abs(speciesSquaredDisplacements_.at(Element("Fe")) - fe) < 1e-10);
    assert(abs(speciesSquaredDisplacements_.at(Element("Ni")) - ni) < 1e-10);
  }

  void SaveState() const
  {
    if (world_rank_ != 0) return;
    std::ofstream out("state.txt");
    out << std::setprecision(17) << steps_ << ' ' << time_ << ' '
        << displacementOriginSteps_ << ' ' << displacementOriginTime_ << '\n';
    for (const auto &d : atomDisplacements_) out << d << '\n';
    for (const auto &[e, s] : speciesSquaredDisplacements_)
      out << s << ' ' << s / speciesAtomCounts_.at(e) << ' '
          << speciesDisplacements_.at(e) << '\n';
    out << vacancyTrajectory_ << '\n';
  }
};

Config MakeConfig()
{
  Eigen::Matrix3Xd positions(3, 4);
  positions << 0.9, 0.1, 0.3, 0.5, 0, 0, 0, 0, 0, 0, 0, 0;
  Config config(10.0 * Eigen::Matrix3d::Identity(), positions,
                {Element("X"), Element("Fe"), Element("Ni"), Element("Fe")});
  config.UpdateNeighborList({2.1});
  return config;
}

int main(int argc, char **argv)
{
  assert(argc == 2);
  const string mode = argv[1];
  const bool restart = mode == "restart" || mode == "new_origin" || mode == "bad_header";
  Config config = restart ? Config::ReadConfig("3.cfg.gz") : MakeConfig();
  if (restart) config.UpdateNeighborList({2.1});
  VacancyMigrationPredictor predictor;
  // Step three is an odd exchange: atom 1 has displacement (-2,0,0).
  map<Element, Eigen::RowVector3d> collective;
  if (restart) collective.emplace(Element("Fe"), Eigen::RowVector3d(-2, 0, 0));
  TracerTest sim(config, 1, 3, mode == "checkpoint" ? 3 : 7, 10,
                 restart ? 3 : 0, 0.0, restart ? 0.75 : 0.0, 600,
                 predictor, "", false,
                 restart ? Eigen::RowVector3d(2, 0, 0) : Eigen::RowVector3d::Zero(),
                 mode == "adaptive" ? "adaptive" : "linear",
                 restart && mode != "new_origin" ? "3.displacements.gz" : "", collective);
  if (mode == "tracking")
  {
    sim.Check(); sim.Expect(0, 0);
    sim.Hop(1); sim.Expect(4, 0);  // Single jump across a periodic boundary.
    sim.Hop(1); sim.Expect(0, 0);  // Exact return cancels the squared contribution.
    sim.Hop(1); sim.Hop(2); sim.Hop(1);
    sim.Expect(4, 4);             // Net +2, not squared jumps 4 + 16 = 20.
    sim.Hop(3);                   // More than one atom of the same species moves.
    for (int i = 0; i < 1000; ++i) sim.Hop(1 + (i * 17) % 3);
    bool rejected = false;
    try { sim.ReadAtomDisplacements("missing.displacements.gz"); }
    catch (const runtime_error &) { rejected = true; }
    assert(rejected);
  }
  else if (mode == "new_origin")
  {
    sim.Expect(0, 0);
    sim.SaveState();
  }
  else if (mode == "bad_header")
  {
    bool rejected = false;
    try { sim.Simulate(); }
    catch (const runtime_error &e)
    { rejected = string(e.what()).find("Incompatible kmc_log.txt header") != string::npos; }
    assert(rejected);
  }
  else
  {
    sim.Check();
    if (restart) sim.Expect(4, 0);
    sim.Simulate();
    sim.Check(); sim.Expect(4, 0);
    sim.SaveState();
  }
}
