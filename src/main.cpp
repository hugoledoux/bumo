
#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include "json.hpp"
#include "definitions.h"
#include "geomtools.h"
#include "Shell.h"

#include <boost/program_options.hpp>

using json = nlohmann::json;

void    list_all_vertices(json& j);
std::vector<Point3> get_coordinates(const json& j, bool translate = true);
void    calculate_metrics(std::vector<Point3>& lspts, const json &j);

std::set<std::string> metrics = {
  "area",
  "circumference",
  "cohesion",
  "convexity",
  "cubeness",
  "cuboidindex",
  // "depth",
  "dispersion",
  "fractality",
  "girth",
  "hemisphericality",
  "proximity",
  "range",
  "rectangularity",
  "roughness",
  "spin",
  "volume"
};

int main(int argc, const char * argv[]) {
  std::string ifile; 
  bool bTranslate = false;
  bool bAutorepair = false;
  bool bVerbose = false;

  try {
    namespace po = boost::program_options;
    po::options_description pomain("Allowed options");
    pomain.add_options()
      ("help", "View all options")
      ("metrics", po::bool_switch(), "List the metrics calculated")
      ("translate", po::bool_switch(), "Use transform/translate (default=false)")
      ("autorepair", po::bool_switch(), "Try to autorepair (default=false)")
      ("verbose", po::bool_switch(), "Verbose output")
      ;
    po::options_description pohidden("Hidden options");
    pohidden.add_options()
      ("inputfile", po::value<std::string>(&ifile), "Input CityJSON file")
      ;        
    po::positional_options_description popos;
    popos.add("inputfile", -1);
    
    po::options_description poall;
    poall.add(pomain).add(pohidden);
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(poall).positional(popos).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << "Usage: bumo myfile.city.json" << std::endl;
      std::cout << pomain << std::endl;
      return 1;
    }
    if (vm.count("license")) {
      std::cout << "GPLv3" << std::endl;
      return 1;
    }
    if (vm.count("version")) {
      std::cout << "bumo v0.1" << std::endl;
      return 1;
    }
    if (vm["metrics"].as<bool>() == true) {
      for (auto& m : metrics) {
        std::cout << m << std::endl;
      }
      return 1;
    }
    if (vm.count("inputfile") == 0) {
      std::cerr << "Error: one input CityJSON file must be specified." << std::endl;
      std::cout << std::endl << pomain << std::endl;
      return 0;  
    }
    //-- store params
    if (vm["translate"].as<bool>() == true) {
      bTranslate = true;
    }
    if (vm["autorepair"].as<bool>() == true) {
      bAutorepair = true;
    }
    if (vm["verbose"].as<bool>() == true) {
      bVerbose = true;
    }
  } 
  catch(std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  } 

  std::ifstream input(ifile);
  json j;
  input >> j;
  input.close();

  std::vector<Point3> lspts = get_coordinates(j, false);

  calculate_metrics(lspts, j);

  return 0;
}


void calculate_metrics(std::vector<Point3>& lspts, const json &j) {
  
  //-- header CSV output
  std::cout << "id[lod],";
  for (auto& metric : metrics) {
    std::cout << metric << ",";
  }
  std::cout << std::endl;

  //-- process each CityObjects (and each of its geoms)
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      if (g["type"] != "Solid") {
        continue;
      }
      std::vector<std::vector<int>> trs;
      for (int i = 0; i < g["boundaries"].size(); i++) {
        for (int j = 0; j < g["boundaries"][i].size(); j++) {
          std::vector<std::vector<int>> gb = g["boundaries"][i][j];
          //-- save the triangles
          std::vector<std::vector<int>> tris = construct_ct_one_face(gb, lspts);
          for (auto& each : tris) {
            trs.push_back(each);
          }
        }
      }
      if (trs.empty() == false) {
        std::cout << std::setprecision(3) << std::fixed;
        Shell s = Shell(trs, lspts);
        std::cout << co.key() << "[" << g["lod"].get<std::string>() << "]" << "," << std::endl;   
        std::cout << "area: " << s.area() << "," << std::endl;
        std::cout << "circumference: " << s.circumference() << "," << std::endl;
        std::cout << "cohesion: " << s.cohesion() << "," << std::endl;
        std::cout << "convexity: " << s.convexity() << "," << std::endl;
        std::cout << "cubeness: " << s.cubeness() << "," << std::endl;
        std::cout << "cuboidindex: " << s.cuboidindex() << "," << std::endl;
        std::cout << "depth: " << s.depth() << "," << std::endl;  
        std::cout << "dispersion: " << s.dispersion() << "," << std::endl;
        std::cout << "fractality: " << s.fractality() << "," << std::endl;
        std::cout << "girth: " << s.girth() << "," << std::endl;
        std::cout << "hemisphericality: " << s.hemisphericality() << "," << std::endl;
        std::cout << "proximity: " << s.proximity() << "," << std::endl;
        std::cout << "range: " << s.range() << "," << std::endl;
        std::cout << "rectangularity: " << s.rectangularity() << "," << std::endl;
        std::cout << "roughness: " << s.roughness() << "," << std::endl;
        std::cout << "spin: " << s.spin() << "," << std::endl;
        std::cout << "volume: " << s.volume() << "," << std::endl;
        std::cout << std::endl;     
        // std::string output_name = "/Users/hugo/temp/" + co.key() + ".off";
        // // std::string output_name = "/Users/hugo/temp/" + co.key() + ".wrap.off";
        // s.write_off(output_name);
      }
    }
  }
}


//--
std::vector<Point3> get_coordinates(const json& j, bool translate) {
  std::vector<Point3> lspts;
  std::vector<std::vector<int>> lvertices = j["vertices"];
  if (translate) {
    for (auto& vi : lvertices) {
      double x = (vi[0] * j["transform"]["scale"][0].get<double>()) + j["transform"]["translate"][0].get<double>();
      double y = (vi[1] * j["transform"]["scale"][1].get<double>()) + j["transform"]["translate"][1].get<double>();
      double z = (vi[2] * j["transform"]["scale"][2].get<double>()) + j["transform"]["translate"][2].get<double>();
      lspts.push_back(Point3(x, y, z));
    } 
  } else {
    for (auto& vi : lvertices) {
      double x = (vi[0] * j["transform"]["scale"][0].get<double>());
      double y = (vi[1] * j["transform"]["scale"][1].get<double>());
      double z = (vi[2] * j["transform"]["scale"][2].get<double>());
      lspts.push_back(Point3(x, y, z));
    }
  }
  return lspts;
}


