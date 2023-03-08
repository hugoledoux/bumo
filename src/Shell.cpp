
#include "Shell.h"


Shell::Shell(std::vector<std::vector<int>> trs, std::vector<Point3> lspts) {
  CGAL::Polygon_mesh_processing::repair_polygon_soup(lspts, trs);
  CGAL::Polygon_mesh_processing::orient_polygon_soup(lspts, trs);
  _lspts = lspts;
  _trs = trs;
  Mesh mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(_lspts, _trs, mesh);   
  _mesh = mesh;
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


Mesh 
Shell::get_mesh() {
  return _mesh;
}

bool 
Shell::is_closed() {
  return CGAL::is_closed(_mesh);
}

double 
Shell::volume() {
  return CGAL::Polygon_mesh_processing::volume(_mesh);
}

double 
Shell::area() {
  return CGAL::Polygon_mesh_processing::area(_mesh);
}


