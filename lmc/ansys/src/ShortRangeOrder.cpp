/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/22/23 10:40 PM                                                          *
 **************************************************************************************************/

#include "ShortRangeOrder.h"

#include <utility>

ShortRangeOrder::ShortRangeOrder(const Config &config, const std::set<Element> &element_set)
    : config_(config), element_set_(element_set) {}
// std::map<std::string, double> ShortRangeOrder::FindPairCorrelationCluster(
//     const size_t shell_number, const std::vector<size_t> &cluster_atom_id_list) const {
//   std::map<std::string, double> warren_cowley;
//   for (const auto &element1: element_set_) {
//     for (const auto &element2: element_set_) {
//       warren_cowley[element1.GetString() + "-" + element2.GetString()] = 0;
//     }
//   }
//   const auto element_list_map_all = config_.GetElementAtomIdVectorMap();
//   std::map<Element, double> concentration;
//   for (const auto type: element_set_) {
//     concentration[type] = static_cast<double>(element_list_map_all.at(type).size())
//         / static_cast<double>(config_.GetNumAtoms());
//   }
//
//   std::map<Element, std::vector<size_t> > element_list_map{};
//   for (const auto &element: element_set_) {
//     element_list_map[element] = {};
//   }
//   for (const auto &atom_id: cluster_atom_id_list) {
//     auto element = config_.GetAtomVector()[atom_id].GetElement();
//     if (element == ElementName::X) { continue; }
//     element_list_map.at(element).push_back(atom_id);
//   }
//   size_t num_bonds;
//   switch (shell_number) {
//     case 1:num_bonds = constants::kNumFirstNearestNeighbors;
//       break;
//     case 2:num_bonds = constants::kNumSecondNearestNeighbors;
//       break;
//     case 3:num_bonds = constants::kNumThirdNearestNeighbors;
//       break;
//     default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
//   }
//   for (const auto &element1: element_set_) {
//     size_t num_all_bonds = element_list_map.at(element1).size() * num_bonds;
//     std::map<Element, size_t> ct_this_pair_map{};
//     for (const auto &element2: element_set_) {
//       ct_this_pair_map[element2] = 0;
//     }
//     for (const auto &atom_id1: element_list_map.at(element1)) {
//       std::vector<size_t> neighbor_list;
//       switch (shell_number) {
//         case 1:neighbor_list = config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id1);
//           break;
//         case 2:neighbor_list = config_.GetSecondNeighborsAtomIdVectorOfAtom(atom_id1);
//           break;
//         case 3:neighbor_list = config_.GetThirdNeighborsAtomIdVectorOfAtom(atom_id1);
//           break;
//         default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
//       }
//       for (const auto &atom_id2: neighbor_list) {
//         const auto &this_element = config_.GetElementAtAtomId(atom_id2);
//         if (this_element == ElementName::X) {
//           continue;
//         }
//         ct_this_pair_map[this_element]++;
//       }
//     }
//     for (auto [element2, ct_this_pair]: ct_this_pair_map) {
//       std::string key = element1.GetString() + "-" + element2.GetString();
//       double pij = static_cast<double>(ct_this_pair) / static_cast<double>(num_all_bonds);
//       warren_cowley[key] = (pij - concentration[element2]) /
//           (static_cast<double>(element1 == element2) - concentration[element2]);
//     }
//   }
//   return warren_cowley;
// }
std::map<std::string, double> ShortRangeOrder::FindWarrenCowley(const size_t shell_number) const {
  std::map<std::string, double> warren_cowley;
  
  // Initializing the warren_cowley map
  for (const auto &element1 : element_set_) {
    for (const auto &element2 : element_set_) {
      if (element1 == ElementName::X || element2 ==  ElementName::X) {
        continue;
      }
      warren_cowley[element1.GetElementString() + "-" + element2.GetElementString()] = 0;
    }
  }
  
  auto element_list_map = config_.GetElementOfAtomIdVectorMap();

  auto concentration = config_.GetConcentration();
  
  auto num_bonds = config_.GetNeighborLatticeIdVectorOfLattice(0,shell_number).size();

  for (const auto &element1 : element_set_) {
    
    if (element1 == ElementName::X){
      continue;
    }
    // size_t num_all_bonds = element_list_map.at(element1).size() * num_bonds;

    size_t num_all_bonds = element_list_map.at(element1).size()*num_bonds;

    std::map<Element, size_t> ct_this_pair_map{};
    for (const auto &element2 : element_set_) {
      if (element2 == ElementName::X) {
        continue;
      }
      ct_this_pair_map[element2] = 0;
    }
    for (const auto &atom_id1 : element_list_map.at(element1)) {
      std::vector<size_t> neighbor_list = config_.GetNeighborAtomIdVectorOfAtom(atom_id1, shell_number);
      // switch (shell_number) {
      //   case 1:neighbor_list = config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id1);
      //     break;
      //   case 2:neighbor_list = config_.GetSecondNeighborsAtomIdVectorOfAtom(atom_id1);
      //     break;
      //   case 3:neighbor_list = config_.GetThirdNeighborsAtomIdVectorOfAtom(atom_id1);
      //     break;
      //   default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
      // }
      for (const auto &atom_id2 : neighbor_list) {
        const auto &this_element = config_.GetElementOfAtom(atom_id2);
        if (this_element == ElementName::X) {
          continue;
        }
        ct_this_pair_map[this_element]++;
      }
    }
    for (auto [element2, ct_this_pair] : ct_this_pair_map) {
      std::string key = element1.GetElementString() + "-" + element2.GetElementString();
      double pij = static_cast<double>(ct_this_pair) / static_cast<double>(num_all_bonds);
      warren_cowley[key] = (pij - concentration[element2]) /
          (static_cast<double>(element1 == element2) - concentration[element2]);
      // res[key] = 1 - (pij / concentration[element2]);
    }
  }
  return warren_cowley;
}

//std::map<std::string, double> ShortRangeOrder::FindProbabilityCluster(
//    size_t shell_number, const std::vector<size_t> &cluster_atom_id_list) const {
//  std::map<std::string, double> pij_map;
//  for (const auto &element1 : element_set_) {
//    for (const auto &element2 : element_set_) {
//      pij_map[element1.GetString() + "-" + element2.GetString()] = 0;
//    }
//  }
//  const auto element_list_map_all = config_.GetElementAtomIdVectorMap();
//  std::map<Element, std::vector<size_t> > element_list_map{};
//  for (const auto &element : element_set_) {
//    element_list_map[element] = {};
//  }
//  for (const auto &atom_id : cluster_atom_id_list) {
//    auto element = config_.GetAtomVector()[atom_id].GetElement();
//    if (element == ElementName::X) { continue; }
//    element_list_map.at(element).push_back(atom_id);
//  }
//  size_t num_bonds;
//  switch (shell_number) {
//    case 1:num_bonds = constants::kNumFirstNearestNeighbors;
//      break;
//    case 2:num_bonds = constants::kNumSecondNearestNeighbors;
//      break;
//    case 3:num_bonds = constants::kNumThirdNearestNeighbors;
//      break;
//    default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
//  }
//  for (const auto &element1 : element_set_) {
//    size_t num_all_bonds = element_list_map.at(element1).size() * num_bonds;
//    std::map<Element, size_t> ct_this_pair_map{};
//    for (const auto &element2 : element_set_) {
//      ct_this_pair_map[element2] = 0;
//    }
//    for (const auto &atom_id1 : element_list_map.at(element1)) {
//      std::vector<size_t> neighbor_list;
//      switch (shell_number) {
//        case 1:neighbor_list = config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id1);
//          break;
//        case 2:neighbor_list = config_.GetSecondNeighborsAtomIdVectorOfAtom(atom_id1);
//          break;
//        case 3:neighbor_list = config_.GetThirdNeighborsAtomIdVectorOfAtom(atom_id1);
//          break;
//        default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
//      }
//      for (const auto &atom_id2 : neighbor_list) {
//        const auto &element2 = config_.GetElementAtAtomId(atom_id2);
//        if (element2 == ElementName::X) {
//          continue;
//        }
//        ct_this_pair_map[element2]++;
//      }
//    }
//    for (auto [element2, ct_this_pair] : ct_this_pair_map) {
//      std::string key = element1.GetString() + "-" + element2.GetString();
//      double pij = static_cast<double>(ct_this_pair) / static_cast<double>(num_all_bonds);
//      pij_map[key] = pij;
//    }
//  }
//  return pij_map;
//}
//std::map<std::string, double> ShortRangeOrder::FindProbability(size_t shell_number) const {
//  std::map<std::string, double> pij_map;
//  for (const auto &element1 : element_set_) {
//    for (const auto &element2 : element_set_) {
//      pij_map[element1.GetString() + "-" + element2.GetString()] = 0;
//    }
//  }
//  const auto &element_list_map = config_.GetElementAtomIdVectorMap();
//
//  size_t num_bonds;
//  switch (shell_number) {
//    case 1:num_bonds = constants::kNumFirstNearestNeighbors;
//      break;
//    case 2:num_bonds = constants::kNumSecondNearestNeighbors;
//      break;
//    case 3:num_bonds = constants::kNumThirdNearestNeighbors;
//      break;
//    default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
//  }
//
//  for (const auto &element1 : element_set_) {
//    size_t num_all_bonds = element_list_map.at(element1).size() * num_bonds;
//    std::map<Element, size_t> ct_this_pair_map{};
//    for (const auto &element2 : element_set_) {
//      ct_this_pair_map[element2] = 0;
//    }
//    for (const auto &atom_id1 : element_list_map.at(element1)) {
//      std::vector<size_t> neighbor_list;
//      switch (shell_number) {
//        case 1:neighbor_list = config_.GetFirstNeighborsAtomIdVectorOfAtom(atom_id1);
//          break;
//        case 2:neighbor_list = config_.GetSecondNeighborsAtomIdVectorOfAtom(atom_id1);
//          break;
//        case 3:neighbor_list = config_.GetThirdNeighborsAtomIdVectorOfAtom(atom_id1);
//          break;
//        default:throw std::invalid_argument("Unknown shell number: " + std::to_string(shell_number));
//      }
//      for (const auto &atom_id2 : neighbor_list) {
//        const auto &element2 = config_.GetElementAtAtomId(atom_id2);
//        if (element2 == ElementName::X) {
//          continue;
//        }
//        ct_this_pair_map[element2]++;
//      }
//    }
//    for (auto [element2, ct_this_pair] : ct_this_pair_map) {
//      std::string key = element1.GetString() + "-" + element2.GetString();
//      double pij = static_cast<double>(ct_this_pair) / static_cast<double>(num_all_bonds);
//      pij_map[key] = pij;
//    }
//  }
//  return pij_map;
//}
//// std::map<std::string, double> ShortRangeOrder::FindWarrenCowley(size_t shell_number) const {
////   std::map<std::string, double> res;
////   size_t num_bonds;
////   switch (shell_number) {
////     case 1:num_bonds = constants::kNumFirstNearestNeighbors;
////       break;
////     case 2:num_bonds = constants::kNumSecondNearestNeighbors;
////       break;
////     case 3:num_bonds = constants::kNumThirdNearestNeighbors;
////       break;
////     default:num_bonds = 0;
////       break;
//   }
//   for (const auto &element1: element_set_) {
//     for (const auto &element2: element_set_) {
//       res[element1.GetString() + "-" + element2.GetString()] = 0;
//     }
//   }
//   const auto &element_list_map = config_.GetElementAtomIdVectorMap();
//   std::map<Element, double> concentration;
//   std::map<Element, double> count;
//   for (const auto type: element_set_) {
//     count[type] = static_cast<double>(element_list_map.at(type).size());
//     concentration[type] = static_cast<double>(element_list_map.at(type).size())
//         / static_cast<double>(config_.GetNumAtoms());
//   }
//   for (const auto &atom1: config_.GetAtomVector()) {
//     const auto element1 = atom1.GetElement();
//     if (element1 == ElementName::X) {
//       continue;
//     }
//     std::vector<size_t> neighbor_list;
//     switch (shell_number) {
//       case 1:neighbor_list = config_.GetFirstNeighborsAtomIdVectorOfAtom(atom1.GetId());
//         break;
//       case 2:neighbor_list = config_.GetSecondNeighborsAtomIdVectorOfAtom(atom1.GetId());
//         break;
//       case 3:neighbor_list = config_.GetThirdNeighborsAtomIdVectorOfAtom(atom1.GetId());
//         break;
//       default:break;
//     }
//     for (const auto &element2: element_set_) {
//       double bond_count = 0;
//       for (auto neighbor_index: neighbor_list) {
//         if (config_.GetAtomVector()[neighbor_index].GetElement() == element2) {
//           bond_count += 1;
//         }
//       }
//       double pij = bond_count / static_cast<double>(num_bonds);
//       double aij = (pij - concentration[element2])
//           / (static_cast<double>(atom1.GetElement() == element2) - concentration[element2]);
//       res[element1.GetString() + "-" + element2.GetString()] += aij / count[element1];
//     }
//   }
//   return res;
// }
 // ansys