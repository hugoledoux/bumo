
#include <iostream>
#include <fstream>
#include <string>
#include "json.hpp"
#include "definitions.h"
#include "geomtools.h"

#include <boost/program_options.hpp>

using json = nlohmann::json;

void    list_all_vertices(json& j);
void    save2obj(std::string filename, std::vector<Point3>& lspts, std::vector<std::vector<int>>& trs);
std::vector<Point3> get_coordinates(const json& j, bool translate = true);
std::vector<std::vector<int>> triangulate_all(const std::vector<Point3>& lspts, const json &j);
void    calculate_metrics(std::vector<Point3>& lspts, const json &j);



int main(int argc, const char * argv[]) {
  std::string ifile; 
  try {
    namespace po = boost::program_options;
    po::options_description pomain("Allowed options");
    pomain.add_options()
      ("help", "View all options")
      ("no_duplicates", "Do *not* remove the duplicate vertices")
      ("no_orphans", "Do *not* remove the orphans")
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
      std::cout << "=== cjcompress help ===" << std::endl;
      std::cout << pomain << std::endl;
      return 0;
    }
    if (vm.count("inputfile") == 0) {
      std::cerr << "Error: one input CityJSON file must be specified." << std::endl;
      std::cout << std::endl << pomain << std::endl;
      return 0;  
    }
  } 
  catch(std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  } 


  std::cout << "Processing: " << ifile << std::endl;

  std::ifstream input(ifile);
  json j;
  input >> j;
  input.close();

  std::vector<Point3> lspts = get_coordinates(j, false);

  calculate_metrics(lspts, j);

  std::vector<std::vector<int>> trs = triangulate_all(lspts, j);
  save2obj("/Users/hugo/temp/z.obj", lspts, trs);

  return 0;
}


void calculate_metrics(std::vector<Point3>& lspts, const json &j) {
  for (auto& co : j["CityObjects"].items()) {
    std::vector<std::vector<int>> trs;
    std::vector<Point3> shellpts;
    bool hasgeom = false;
    for (auto& g : co.value()["geometry"]) {
      hasgeom = true;
      for (int i = 0; i < g["boundaries"].size(); i++) {
        for (int j = 0; j < g["boundaries"][i].size(); j++) {
          std::vector<std::vector<int>> gb = g["boundaries"][i][j];
          //-- save the triangles
          std::vector<std::vector<int>> tris = construct_ct_one_face(gb, lspts);
          for (auto& each : tris) {
            trs.push_back(each);
          }
          //-- save the unique points
          std::set<size_t> uids;
          for (auto& ring : gb) {
            for (auto& pi: ring) {
              uids.insert(pi);
            }
          }
          for (auto& each : uids) {
            shellpts.push_back(lspts[each]);
          }
        }
      }
    }
    if (hasgeom) {
      std::cout << co.key() << " ";
      // std::cout << "==> " << co.key() << std::endl;
      
      double vol = volume_shell(trs, lspts);
      std::cout << std::setprecision(2) << std::fixed << vol << " ";
      // std::cout << "vol: " << std::setprecision(2) << std::fixed << vol << std::endl;

      double area = area_shell(trs, lspts);
      // std::cout << "area: " << area << std::endl;

      double voloobb = volume_oobb(shellpts);
      // std::cout << "vol oobb: " << voloobb << std::endl;

      std::cout << (vol / voloobb) << " ";
      // std::cout << "rectangularity: "   << std::setprecision(3) << std::fixed << (vol / voloobb) << std::endl;

      std::cout << (3*sqrt(2*3.14159)*vol) / pow(area, 1.5) << " ";
      // std::cout << "hemisphericality: " << (3*sqrt(2*3.14159)*vol) / pow(area, 1.5) << std::endl;

      double mu1 = mu(shellpts, trs, lspts);
      // std::cout << "mu: " << mu1 << std::endl;
      double roughness = pow(mu1, 3.0) * 48.735 / (vol + pow(area, 1.5));
      std::cout << roughness << std::endl;
      // std::cout << "roughness: " << roughness << std::endl;

      // std::cout << "orientations:" << std::endl;
      // roof_orientation(co.key(), lspts, j);
          
      // std::cout << "---" << std::endl;
    }
  }
}


void roof_orientation(std::string coid, std::vector<Point3>& lspts, const json &j) {
  auto co = j["CityObjects"][coid];
  for (auto& g : co["geometry"]) {
    if (g["type"] == "Solid") { //-- if MultiSurface then one less for-loop
      for (int i = 0; i < g["boundaries"].size(); i++) {
        for (int j = 0; j < g["boundaries"][i].size(); j++) {
          if (g.contains("semantics")) {
            int sem_index = g["semantics"]["values"][i][j];
            if (g["semantics"]["surfaces"][sem_index]["type"].get<std::string>().compare("RoofSurface") == 0) {
//                std::cout << "RoofSurface: " << g["boundaries"][i][j] << std::endl;
              std::vector<std::vector<int>> tris = construct_ct_one_face(g["boundaries"][i][j], lspts);
              if (tris.empty()) {
                std::cout << "crash" << std::endl;
                continue;
              }
              auto lastr = tris.back();
              auto n = CGAL::unit_normal(lspts[lastr[0]], lspts[lastr[1]], lspts[lastr[2]]);
              if ( (n[2] > 0.99) || (n[2] < -0.99) ) {
                std::cout << "horizontal" << std::endl;
              } else {
                double a = atan2(n[1], n[0]);
                // std::cout << "atan2 " << a << std::endl;
                if ( (a >= 0.0) && (a < (M_PI/4)) )
                  std::cout << "EN" << std::endl;
                else if ( (a >= (M_PI/4)) && (a < (M_PI/2)) )
                  std::cout << "NE" << std::endl;
                else if ( (a >= (M_PI/2)) && (a < (M_PI*3/4)) )
                  std::cout << "NW" << std::endl;
                else if ( (a >= (M_PI*3/4)) && (a <= (M_PI)) )
                  std::cout << "WN" << std::endl;
                else if ( (a < 0.0) && (a >= -(M_PI/4)) )
                  std::cout << "ES" << std::endl;
                else if ( (a < -(M_PI/4)) && (a >= -(M_PI/2)) )
                  std::cout << "SE" << std::endl;
                else if ( (a < -(M_PI/2)) && (a >= -(M_PI*3/4)) )
                  std::cout << "SW" << std::endl;
                else if ( (a < -(M_PI*3/4)) && (a >= -(M_PI)) )
                  std::cout << "WS" << std::endl;  
              }
            }
          }
        }
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



std::vector<std::vector<int>> triangulate_all(const std::vector<Point3>& lspts, const json &j) {
  std::vector<std::vector<int>> trs;
  // auto g = j["CityObjects"]["NL.IMBAG.Pand.0503100000018416-0"]["geometry"][0];
  // auto g = j["CityObjects"]["id-1"]["geometry"][0];
  for (auto& co : j["CityObjects"].items()) {
    for (auto& g : co.value()["geometry"]) {
      // if (g["type"] == "Solid") {
      if ( (g["type"] == "Solid") && (g["lod"] == "2.2") ) {   //-- LoD2.2 only!!!!!
        for (int i = 0; i < g["boundaries"].size(); i++) {
          for (int j = 0; j < g["boundaries"][i].size(); j++) {
            std::vector<std::vector<int>> gb = g["boundaries"][i][j];
            // std::cout << "surface: " << g["boundaries"][i][j][0] << std::endl;
            std::vector<std::vector<int>> tris = construct_ct_one_face(gb, lspts);
            for (auto& each : tris) {
              trs.push_back(each);
            }
          }
        }
      }
    }
  }
  return trs;
}


//-- write the OBJ file
void save2obj(std::string filename, std::vector<Point3>& lspts, std::vector<std::vector<int>>& trs) {
  std::ofstream ofile(filename);
  for (auto& p : lspts)
    ofile << std::setprecision(5) << std::fixed << "v " << p.x() << " " << p.y() << " " << p.z() << std::endl;
  for (auto& f : trs) {
    ofile << "f " << (f[0] + 1) << " " << (f[1] + 1) << " " << (f[2] + 1) << std::endl;
  }
  ofile.close();
}



