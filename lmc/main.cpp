/*! \file  main.cpp
 *  \brief File for the main function.
 */

#include "Home.h"
/*
int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);

  
}
  */

int main()
{

  // auto cfg = Config::ReadXYZ("//media/sf_Phd/Structures/TiMo_Supercell.xyz");
  auto cfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");
  Config::WriteXyzExtended("//media/sf_Phd/Structures/Mo_test.xyz.gz", cfg, {}, {});
  Config::WriteConfig("//media/sf_Phd/Structures/Mo_test.cfg.gz", cfg);
  auto cfg2 = Config::ReadXYZ("//media/sf_Phd/Structures/Mo_test.xyz.gz");

}