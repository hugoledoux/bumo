
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

  CGAL::alpha_wrap_3(_mesh_original, 1.0, 0.3, _mesh_wrap); //-- values of Ivan
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


double
Shell::rectangularity() {
  auto o = this->get_oobb();
  double voloobb = oobb_volume(o);
  return (this->volume() / voloobb);
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


void
Shell::write_off(std::string s) {
   CGAL::IO::write_polygon_mesh(s, *_mesh, CGAL::parameters::stream_precision(17));
}