#include "WidomInsertion.h"

vector<WidomInsertionRecord> WidomInsertion(
    Config &config,
    SymmetricCEPredictor &symCEPredictor,
    Element ghostElement)
{
  const auto numSites = config.GetNumLattices();
  vector<WidomInsertionRecord> records;
  records.reserve(numSites);

  for (size_t latticeId = 0; latticeId < numSites; latticeId++)
  {
    auto elementAtLatticeId = config.GetElementOfLattice(latticeId);

    double localFormationEnergyBefore =
        symCEPredictor.ComputeLocalFormationEnergyOfSite(config,
                                                         latticeId);

    // assign the ghostElement
    config.SetElementOfLattice(latticeId, ghostElement);

    // compute the local formation energy
    double localFormationEnergyAfter =
        symCEPredictor.ComputeLocalFormationEnergyOfSite(config,
                                                         latticeId);

    // reassign the previous element
    config.SetElementOfLattice(latticeId, elementAtLatticeId);

    double localFormationEnergyChange = localFormationEnergyAfter - localFormationEnergyBefore;

    records.push_back({latticeId, elementAtLatticeId, localFormationEnergyChange});
  }

  return records;
}

void WidomInsertionHelper(
    unsigned long long initialSteps,
    unsigned long long incrementSteps,
    SymmetricCEPredictor &symCEPredictor,
    Element ghostElement)
{
  fs::create_directories("widomInsertion");

  unsigned long long finalStep = initialSteps;

  while (true)
  {
    string f = to_string(finalStep) + ".cfg.gz";
    if (!fs::exists(f))
      break;
    finalStep += incrementSteps;
  }

  // Handle "no files exist at initialSteps"
  if (finalStep == initialSteps)
  {
    cerr << "[Widom] No config files found starting at step "
         << initialSteps << "\n";
    return;
  }

  // finalStep is now the first missing step
  finalStep -= incrementSteps;

  for (unsigned long long idx = initialSteps; idx <= finalStep; idx += incrementSteps)
  {
    string configFile = to_string(idx) + ".cfg.gz";

    if (!fs::exists(configFile))
    {
      cerr << "[Widom] Skipping step " << idx
           << " because config file does not exist: " << configFile << "\n";
      continue;
    }

    try
    {
      Config cfg = Config::ReadConfig(configFile);

      // If needed (recommended if CE uses neighbors):
      // cfg.UpdateNeighborList(cutoffs);

      auto widomInsertionRecords = WidomInsertion(cfg, symCEPredictor, ghostElement);

      std::string outPath = "widomInsertion/" + std::to_string(idx) + ".txt.gz";

      boost::iostreams::filtering_ostream outfile;
      outfile.push(boost::iostreams::gzip_compressor());
      outfile.push(boost::iostreams::file_sink(outPath, std::ios_base::binary));

      if (!outfile.good())
      {
        throw std::runtime_error("Cannot open output file: " + outPath);
      }

      outfile << "latticeId\telement\\tdeltaE_ghost_insert\n";
      outfile << setprecision(16);

      for (const auto &entry : widomInsertionRecords)
      {
        outfile << entry.latticeId << "\t"
                << entry.originalElement.GetElementString() << "\t"
                << entry.localFormationEnergyChange << "\n";
      }

      outfile.flush();

      cout << "Done for " << configFile << endl;
    }
    catch (const exception &e)
    {
      cerr << "[Widom] step " << idx << " failed: " << e.what() << "\n";
    }
  }
}
