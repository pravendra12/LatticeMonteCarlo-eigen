#ifndef LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#define LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_

#include <omp.h>
#include <filesystem>
#include <nlohmann/json.hpp>
#include "ShortRangeOrder.h"
#include "B2OrderParameter.h"
#include "ClusterDynamics.h"
#include "ConfigEncoding.h"
#include <memory>
#include "AnsysFlags.h"

namespace fs = std::filesystem;

namespace ansys
{

  class Traverse
  {
  public:
    Traverse(unsigned long long int initial_steps,
             unsigned long long int increment_steps,
             const std::vector<double> &cutoffs,
             const AnsysFlags &ansys_flags,
             std::string log_type,
             std::string config_type);
    virtual ~Traverse();
    void RunAnsys() const;
    //  void RunReformat() const;

    void RunAnsysOnConfig(
        const Config &config,
        const set<Element> &element_set,
        ostringstream &oss) const;

  private:
    std::string GetHeaderFrameString(const std::set<Element> &element_set) const;

    const unsigned long long initial_steps_;
    const unsigned long long increment_steps_;
    unsigned long long final_steps_;
    const std::vector<double> cutoffs_;
    const AnsysFlags ansys_flags_;
    const std::string log_type_;
    const std::string config_type_;

    using MapVariant =
        std::variant<std::unordered_map<unsigned long long, double>, std::unordered_map<unsigned long long, std::string>>;

    std::unordered_map<std::string, MapVariant> log_map_;

    mutable std::ofstream frame_ofs_;
  };

} // namespace ansys

#endif // LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_