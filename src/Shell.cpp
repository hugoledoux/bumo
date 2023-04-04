
#include "Shell.h"
#include "geomtools.h"


Shell::Shell(std::vector<std::vector<int>> trs, std::vector<Point3> lspts) {
  CGAL::Polygon_mesh_processing::repair_polygon_soup(lspts, trs);
  CGAL::Polygon_mesh_processing::orient_polygon_soup(lspts, trs);
  _lspts = lspts;
  _trs = trs;
  Mesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(_lspts, _trs, mesh);   
  _mesh_original = mesh;
  if (CGAL::is_closed(_mesh_original) == false) {
    //-- attempt to fill holes
    std::vector<halfedge_descriptor> border_cycles;
    //-- collect one halfedge per boundary cycle
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(_mesh_original, std::back_inserter(border_cycles));
    for(halfedge_descriptor h : border_cycles)
    {
      std::vector<face_descriptor>  patch_facets;
      std::vector<vertex_descriptor> patch_vertices;
      bool success = std::get<0>(
        CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(_mesh_original,
                                                                        h,
                                                                        std::back_inserter(patch_facets),
                                                                        std::back_inserter(patch_vertices)));
    }
  }
  //-- if still not closed: create the alph-wrap mesh to compute in/out
  _mesh_wrap = Mesh();
  if (CGAL::is_closed(_mesh_original) == false) {
    CGAL::alpha_wrap_3(_mesh_original, 1.3, 0.3, _mesh_wrap); //-- values of Ivan
    _mesh = &_mesh_wrap;
    std::cout << "use_wrap_mesh!" << std::endl; // TODO: should we use wrap-alpha if invalid?

  } else {
    _mesh = &_mesh_original;
  }
  _mesh = &_mesh_original;
  //-- area+volume
  _area = CGAL::Polygon_mesh_processing::area(_mesh_original);
  _volume = CGAL::Polygon_mesh_processing::volume(_mesh_original);
  //-- samples_surfaces
  CGAL::Polygon_mesh_processing::sample_triangle_soup(_lspts, 
                          _trs, 
                          std::back_inserter(_samples_surface), 
                          CGAL::parameters::number_of_points_per_area_unit(2).
                          use_random_uniform_sampling(true)
  );
  // std::cout << "_samples_surface: " << _samples_surface.size() << std::endl;

  //-- samples_volume
  auto bbox = this->get_aabb();
  auto rand = CGAL::Random();
  Mesh* m;
  if (CGAL::is_closed(_mesh_original) == true) {
    m = &_mesh_original;
  } else {
    m = &_mesh_wrap;
  }
  CGAL::Side_of_triangle_mesh<Mesh, K> inside(*m);
  int n = 0;
  int total = int(_volume * 4.0); //-- 4pts/m^3
  while (n < total) {
    double x = rand.uniform_real(bbox.xmin(), bbox.xmax());
    double y = rand.uniform_real(bbox.ymin(), bbox.ymax());
    double z = rand.uniform_real(bbox.zmin(), bbox.zmax());
    Point3 p(x, y, z);
    if (inside(p) == CGAL::ON_BOUNDED_SIDE) { 
      _samples_volume.push_back(p);
      n++;
    }
  }
  // std::cout << "_samples_volume.size() " << _samples_volume.size() << std::endl;
  // for (auto& p : _samples_volume)
    // std::cout << p << std::endl;
}




//-- returns the radius
double Shell::largest_sphere_inside_mesh() {
  // https://github.com/CGAL/cgal/blob/master/Polygon_mesh_processing/test/Polygon_mesh_processing/test_pmp_distance.cpp
  double re = CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Parallel_if_available_tag>(
    _samples_volume, 
    _mesh_original);
  return re;

}


void
Shell::compute_wrap_mesh() {
  // const double relative_alpha = 300.0;
  // const double relative_offset = 1200.0;
  // auto bbox = this->get_aabb();
  // const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
  //                                      CGAL::square(bbox.ymax() - bbox.ymin()) +
  //                                      CGAL::square(bbox.zmax() - bbox.zmin()));
  // const double alpha = diag_length / relative_alpha;
  // const double offset = diag_length / relative_offset;
  // CGAL::alpha_wrap_3(_mesh_original, alpha, offset, _mesh_wrap);
  CGAL::alpha_wrap_3(_mesh_original, 1.3, 0.3, _mesh_wrap); //-- values of Ivan
}


void 
Shell::use_wrap_mesh(bool b) {
  if (b == true) { 
    _mesh = &_mesh_wrap;
  } else {
    _mesh = &_mesh_original;
  }
}

Polyhedron            
Shell::get_convex_hull() {
  Polyhedron poly;
  CGAL::convex_hull_3(_lspts.begin(), _lspts.end(), poly);
  return poly;
}

std::array<Point3,8> 
Shell::get_oobb() {
  std::array<Point3, 8> obb_points;
  CGAL::oriented_bounding_box(_lspts, obb_points);
  return obb_points;
}

K::Iso_cuboid_3 
Shell::get_aabb() {
  return CGAL::bounding_box(_lspts.begin(), _lspts.end());
}

Mesh* 
Shell::get_mesh() {
  return _mesh;
}

bool 
Shell::is_closed() {
  return CGAL::is_closed(*_mesh);
}

double 
Shell::volume() {
  return _volume;
}

double 
Shell::area() {
  return _area;
}

Point3
Shell::centroid() {
  return CGAL::centroid(_samples_surface.begin(), _samples_surface.end());
}

double
Shell::circumference() {
  double re = 4 * 3.14159 * pow(3 * _volume / (4 * 3.14159), 2.0/3.0) / _area;
  return re;
}

double
Shell::cohesion() {
  double totaldistance = 0.0;
  for (auto i = 0; i < _samples_volume.size(); i+=10) {
    for (auto j = 0; j < _samples_volume.size(); j+=10) {
      totaldistance += sqrt(CGAL::squared_distance(_samples_volume[i], _samples_volume[j]));
    }
  }
  double re = 36 / 35 * pow(3 * _volume / (4 * 3.14159), 1.0/3.0) / (1 / pow(_samples_volume.size()/10.0, 2.0) * totaldistance);
  return re;
}


double
Shell::convexity() {
  Polyhedron p = this->get_convex_hull();
  double vol = CGAL::Polygon_mesh_processing::volume(p);
  return (_volume / vol);
}

double
Shell::cubeness() {
  return ( 6 * pow(_volume, 2.0/3.0) / _area );
}

double
Shell::cuboidindex() {
  auto o = this->get_oobb();
  double voloobb = oobb_volume(o);
  double areaoobb = oobb_area(o);
  return pow(_volume / voloobb, 2.0/3.0) * areaoobb / _area;
}

double
Shell::depth() {
  return (4 * this->avg_dist_samples_volume_surface() / pow(3 * _volume / 4 / 3.14159, 1.0/3.0));
}


double
Shell::dispersion() {
  double re = 1 - (this->avg_dist_samples_surface_radius_sphere() / this->avg_dist_samples_surface_centroid());
  return re;
}


double
Shell::fractality() {
  double re = 1 - (std::log(_volume) / 1.5 / std::log(_area));
  return re;
}

double
Shell::girth() {
  double re = this->largest_sphere_inside_mesh() / pow(3 * _volume / 4 / 3.14159, 1.0/3.0);
  return re;
}


double
Shell::hemisphericality() {
  return ( (3 * sqrt(2 * 3.14159) * _volume) / pow(_area, 1.5) );
}

double
Shell::proximity() {
  double re = 0.75 * pow(3 * _volume / 4 / 3.14159, 1.0/3.0) / this->avg_dist_samples_volume_centroid();
  return re;
}

double
Shell::range() {
  //-- get smallest enclosing spheres
  Min_sphere ms(_lspts.begin(), _lspts.end());
  double re = pow(3 * _volume / (4 * 3.14159), 1.0/3.0) / ms.radius();
  return re;
}

double
Shell::rectangularity() {
  auto o = this->get_oobb();
  double voloobb = oobb_volume(o);
  return (_volume / voloobb);
}

double
Shell::roughness() {
  return (pow(this->avg_dist_samples_surface_centroid(), 3.0) * 48.735 / (_volume + pow(_area, 1.5)));
}

double
Shell::spin() {
  double re = 0.6 * pow(3 * _volume / 4 / 3.14159, 2.0/3.0) / this->avg_sq_dist_samples_volume_centroid();
  return re;
}

double
Shell::avg_dist_samples_surface_radius_sphere() {
  Point3 c = CGAL::centroid(_samples_surface.begin(), _samples_surface.end());
  double r = get_sphere_radius_from_volume(_volume);
  double distance = 0.0;
  int total = 0;   
  for (auto& s : _samples_surface) {
    total += 1;
    distance += abs(sqrt(CGAL::squared_distance(c, s)) - r);
  }
  return (distance / total);
}


double
Shell::avg_dist_samples_surface_centroid() {
  Point3 c = CGAL::centroid(_samples_surface.begin(), _samples_surface.end());
  double distance = 0.0;
  int total = 0;     
  for (auto& s : _samples_surface) {
    total += 1;
    distance += sqrt(CGAL::squared_distance(c, s));
  }
  return (distance / total);
}

double
Shell::avg_dist_samples_volume_centroid() {
  Point3 c = CGAL::centroid(_samples_surface.begin(), _samples_surface.end());
  double distance = 0.0;
  int total = 0;     
  for (auto& s : _samples_volume) {
    total += 1;
    distance += sqrt(CGAL::squared_distance(c, s));
  }
  return (distance / total);
}

double
Shell::avg_dist_samples_volume_surface() {
  double distance = 0.0;
  int count = 0;
  KDTree kdtree(_samples_surface.begin(), _samples_surface.end());
  for (auto i = 0; i < _samples_volume.size(); i+=10) {
    Neighbor_search search(kdtree, _samples_volume[i], 1);
    distance += std::sqrt(search.begin()->second);
    count++;
  }
  return (distance / count);
}

// Shell::avg_dist_samples_volume_surface() {
//   double distance = 0.0;
//   int count = 0;
//   // for (auto& s : _samples_volume) {
//   for (auto i = 0; i < _samples_volume.size(); i+=10) {
//     std::vector<Point3> v;
//     v.push_back(_samples_volume[i]);
//     double d = CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Parallel_if_available_tag>(
//                 v, 
//                 _mesh_original);
//     distance += d;
//     count += 1;
//   }
//   return (distance / count);
// }

double
Shell::avg_sq_dist_samples_volume_centroid() {
  Point3 c = CGAL::centroid(_samples_surface.begin(), _samples_surface.end());
  double distance = 0.0;
  int total = 0;     
  for (auto& s : _samples_volume) {
    total += 1;
    distance += CGAL::squared_distance(c, s);
  }
  return (distance / total);
}

void
Shell::write_off(std::string s) {
   CGAL::IO::write_polygon_mesh(s, *_mesh, CGAL::parameters::stream_precision(17));
}