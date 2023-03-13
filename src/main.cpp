
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
void test1();

std::set<std::string> metrics = {
  "area", 
  "cubeness", 
  "hemisphericality",
  "rectangularity",
  "roughness",
  "total_vertices",
  "total_surfaces",
  "volume", 
  "volume_aabb", 
  "volume_oobb" 
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



  // std::cout << "Processing: " << ifile << std::endl;

  std::ifstream input(ifile);
  json j;
  input >> j;
  input.close();

  std::vector<Point3> lspts = get_coordinates(j, false);

  calculate_metrics(lspts, j);
  // test1();

  return 0;
}


void test1() {
  Mesh mesh;
  // CGAL::IO::read_polygon_mesh("/Users/hugo/temp/NL.IMBAG.Pand.0503100000031316-0.off", mesh);
  // CGAL::IO::read_polygon_mesh("/Users/hugo/temp/NL.IMBAG.Pand.0503100000032914-0_goodvol.off", mesh);
  CGAL::IO::read_polygon_mesh("/Users/hugo/temp/NL.IMBAG.Pand.0503100000032914-0_badvol.off", mesh);
  std::vector<Point3> lspts;
  std::vector<std::vector<int>> trs;
  CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(mesh, lspts, trs); 
  std::cout << lspts.size() << "  " << trs.size() << std::endl;
  double v = volume_shell(trs, lspts);
  std::cout << "vol: " << v << std::endl;

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
      std::vector<std::vector<int>> trs;
      // std::vector<Point3> shellpts;
      for (int i = 0; i < g["boundaries"].size(); i++) {
        for (int j = 0; j < g["boundaries"][i].size(); j++) {
          std::vector<std::vector<int>> gb = g["boundaries"][i][j];
          //-- save the triangles
          std::vector<std::vector<int>> tris = construct_ct_one_face(gb, lspts);
          for (auto& each : tris) {
            trs.push_back(each);
          }
          // //-- save the unique points
          // std::set<size_t> uids;
          // for (auto& ring : gb) {
          //   for (auto& pi: ring) {
          //     uids.insert(pi);
          //   }
          // }
          // for (auto& each : uids) {
          //   shellpts.push_back(lspts[each]);
          // }
        }
      }
      if (trs.empty() == false) {
        std::cout << co.key() << "[" << g["lod"].get<std::string>() << "]" << "," << std::endl;   
        Shell s = Shell(trs, lspts);
        std::cout << std::setprecision(3) << std::fixed;

        double area = area_shell(trs, lspts);
        double volume = volume_shell(trs, lspts);

        double area1 = s.area();
        double volume1 = s.volume();

        // double d = s.largest_sphere_inside_mesh();
        // std::cout << "largest_sphere_inside_mesh: " << d << std::endl;

        double cohesion = s.cohesion();
        std::cout << "cohesion: " << cohesion << std::endl;
        
        // if ( CGAL::is_closed(s.get_mesh()) == false ) {
          std::cout << "area: " << area << " | " << area1 << std::endl;
          std::cout << "volume: " << volume << " | " << volume1 << std::endl;
          std::cout << "rectangularity: " << s.rectangularity() << std::endl;
        // }

        std::cout << "range: " << s.range() << std::endl;

        // s.compute_wrap_mesh();
        // s.use_wrap_mesh(true);
        // std::cout << "area wrap: " << s.area()<< std::endl;
        // std::cout << "volume wrap: " << s.volume() << std::endl;

        // s.fill_holes();
        // std::cout << "is_closed: " << s.is_closed() << std::endl;

        // std::cout << "area: " << s.area()<< std::endl;
        // std::cout << "volume: " << s.volume() << std::endl;
        // std::cout << "rectangularity: " << s.rectangularity() << std::endl;

        std::string output_name = "/Users/hugo/temp/" + co.key() + ".off";
        // std::string output_name = "/Users/hugo/temp/" + co.key() + ".wrap.off";
        s.write_off(output_name);

        // double vol_oobb = volume_oobb(shellpts);

        // std::cout << std::setprecision(2) << std::fixed;
        // for (auto& metric : metrics) {   
        //   if (metric == "volume") {
        //     std::cout << volume << ",";
        //   } else if (metric == "area") {
        //     std::cout << area << ",";
        //   } else if (metric == "roughness") {
        //     double mu1 = mu(shellpts, trs, lspts);
        //     double roughness = pow(mu1, 3.0) * 48.735 / (volume + pow(area, 1.5));
        //     std::cout << roughness << ",";
        //   }

        // // std::cout << (vol / voloobb) << " ";
              
        // }
        // std::cout << std::endl;
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


