
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
  _mesh_wrap = Mesh();
  _mesh = &_mesh_original;
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
Shell::fill_holes() {
  std::vector<halfedge_descriptor> border_cycles;
  //-- collect one halfedge per boundary cycle
  CGAL::Polygon_mesh_processing::extract_boundary_cycles(*_mesh, std::back_inserter(border_cycles));
  unsigned int nb_holes = 0;
  for(halfedge_descriptor h : border_cycles)
  {
    std::vector<face_descriptor>  patch_facets;
    std::vector<vertex_descriptor> patch_vertices;
    bool success = std::get<0>(
      CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(*_mesh,
                                                                      h,
                                                                      std::back_inserter(patch_facets),
                                                                      std::back_inserter(patch_vertices)));
    std::cout << "* Number of facets in constructed patch: " << patch_facets.size() << std::endl;
    std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
    std::cout << "  Is fairing successful: " << success << std::endl;
    ++nb_holes;
  }
  std::cout << std::endl;
  std::cout << nb_holes << " holes have been filled" << std::endl;
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
  return CGAL::Polygon_mesh_processing::volume(*_mesh);
}

double 
Shell::area() {
  return CGAL::Polygon_mesh_processing::area(*_mesh);
}

Point3
Shell::centroid() {
  std::vector<Point3> samples;
  CGAL::Polygon_mesh_processing::sample_triangle_soup(_lspts, 
                          _trs, 
                          std::back_inserter(samples), 
                          CGAL::parameters::do_sample_vertices(true).
                          use_random_uniform_sampling(true)
  );
  return CGAL::centroid(samples.begin(), samples.end());
}

double
Shell::convexity() {
  Polyhedron p = this->get_convex_hull();
  double vol = CGAL::Polygon_mesh_processing::volume(p);
  return (this->volume() / vol);
}

double
Shell::cubeness() {
  return ( 6 * pow(this->volume(), 2.0/3.0) / this->area() );
}

double
Shell::dispersion() {
  double re = 1 - (this->mu_surface_sphere() / this->mu_surface_centroid());
  return re;
}


double
Shell::fractality() {
  double re = 1 - (std::log(this->volume()) / 1.5 / std::log(this->area()));
  return re;
}

double
Shell::hemisphericality() {
  return ( (3 * sqrt(2 * 3.14159) * this->volume()) / pow(this->area(), 1.5) );
}

double
Shell::range() {
  //-- get smallest enclosing spheres
  Min_sphere ms(_lspts.begin(), _lspts.end());
  double re = pow(3 * this->volume() / (4 * 3.14159), 1.0/3.0) / ms.radius();
  return re;
}

double
Shell::rectangularity() {
  auto o = this->get_oobb();
  double voloobb = oobb_volume(o);
  return (this->volume() / voloobb);
}

double
Shell::roughness() {
  return (pow(this->mu_surface_centroid(), 3.0) * 48.735 / (this->volume() + pow(this->area(), 1.5)));
}

double
Shell::mu_surface_sphere() {
  std::vector<Point3> samples;
  CGAL::Polygon_mesh_processing::sample_triangle_soup(_lspts, 
                          _trs, 
                          std::back_inserter(samples), 
                          CGAL::parameters::do_sample_vertices(true).
                          use_random_uniform_sampling(true)
                          );
  Point3 c = CGAL::centroid(samples.begin(), samples.end());
  double r = get_sphere_radius_from_volume(this->volume());
  double distance = 0.0;
  int total = 0;   
  for (auto& s : samples) {
    total += 1;
    distance += abs(sqrt(CGAL::squared_distance(c, s)) - r);
  }
  return (distance / total);
}


double
Shell::mu_surface_centroid() {
  std::vector<Point3> samples;
  CGAL::Polygon_mesh_processing::sample_triangle_soup(_lspts, 
                          _trs, 
                          std::back_inserter(samples), 
                          CGAL::parameters::do_sample_vertices(true).
                          use_random_uniform_sampling(true)
                          );
  Point3 c = CGAL::centroid(samples.begin(), samples.end());
  double distance = 0.0;
  int total = 0;     
  for (auto& s : samples) {
    total += 1;
    distance += sqrt(CGAL::squared_distance(c, s));
  }
  return (distance / total);
}

void
Shell::write_off(std::string s) {
   CGAL::IO::write_polygon_mesh(s, *_mesh, CGAL::parameters::stream_precision(17));
}