#include "Home.h"
#include "Config.h"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <nlohmann/json.hpp> // https://github.com/nlohmann/json

namespace fs = std::filesystem;
using namespace std;
using json = nlohmann::json;

int main()
{
  string pathOutput = "/home/pravendra3/Documents/ceData";
  string jsonFilePath = pathOutput + "/cutoffs.json";

  // ---------------- Read JSON ----------------
  ifstream inFile(jsonFilePath);
  if (!inFile.is_open())
  {
    cerr << "Cannot open JSON file: " << jsonFilePath << endl;
    return 1;
  }

  json cutoffDictJson;
  try
  {
    inFile >> cutoffDictJson;
  }
  catch (const json::parse_error &e)
  {
    cerr << "Failed to parse JSON: " << e.what() << endl;
    return 1;
  }
  inFile.close();

  map<string, vector<double>> cutoffDict;
  for (auto &[key, value] : cutoffDictJson.items())
  {
    try
    {
      vector<double> cutoffs = value.get<vector<double>>();
      cutoffDict[key] = cutoffs;
    }
    catch (const json::type_error &e)
    {
      cerr << "Invalid data for key " << key << ": " << e.what() << endl;
    }
  }

  // ---------------- Print map ----------------
  for (const auto &[key, cutoffs] : cutoffDict)
  {
    cout << key << ": ";
    for (double c : cutoffs)
      cout << c << " ";
    cout << endl;
  }

  // ---------------- Iterate over subdirectories ----------------
  for (const auto &entry : fs::directory_iterator(pathOutput))
  {
    if (entry.is_directory())
    {
      string configFolder = entry.path().string();
      string folderName = entry.path().filename().string(); // e.g., "ceConfig_0"
      string cfgFilePath = configFolder + "/" + folderName + ".POSCAR";

      if (fs::exists(cfgFilePath))
      {
        // Optional: check file size before reading
        if (fs::file_size(cfgFilePath) == 0)
        {
          cerr << "Warning: Empty cfg file: " << cfgFilePath << endl;
          continue;
        }

        try
        {
          Config cfg = Config::ReadPoscar(cfgFilePath);
          cout << "Read configuration from: " << cfgFilePath << endl;

          auto cutoffs = cutoffDict.at(folderName);

          cfg.UpdateNeighborList(cutoffs);

          cout << cfg.GetNeighborLatticeIdVectorOfLattice(0, 1).size() << endl;
          cout << cfg.GetNeighborLatticeIdVectorOfLattice(0, 2).size() << endl;
          cout << cfg.GetNeighborLatticeIdVectorOfLattice(0, 3).size() << endl;

          // Access cfg data here
        }
        catch (const std::exception &e)
        {
          cerr << "Failed to read cfg: " << cfgFilePath << " Error: " << e.what() << endl;
        }
      }
      else
      {
        cout << "Config file not found: " << cfgFilePath << endl;
      }
    }
  }

  return 0;
}
