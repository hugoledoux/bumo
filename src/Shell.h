

#include "definitions.h"

class Shell {
public:
  Shell(std::vector<std::vector<int>> trs, std::vector<Point3> lspts);

  void                  compute_wrap_mesh();
  void                  use_wrap_mesh(bool b);

  void                  fill_holes();

  Mesh*                 get_mesh();
  Polyhedron            get_convex_hull();
  std::array<Point3,8>  get_oobb();
  K::Iso_cuboid_3       get_aabb();
  
  bool                  is_closed();

  double                area();
  Point3                centroid();
  double                cohesion();
  double                convexity();
  double                cubeness();
  double                dispersion();
  double                fractality();
  double                hemisphericality();
  double                range();
  double                rectangularity();
  double                roughness();
  double                volume();


  void                  write_off(std::string s);

private:
  std::vector<Point3>           _lspts;
  std::vector<std::vector<int>> _trs;
  
  Mesh                          _mesh_original;
  Mesh                          _mesh_wrap;
  Mesh*                         _mesh;

  double                        _area;
  double                        _volume;
  std::vector<Point3>           _samples_surface;
  std::vector<Point3>           _samples_volume;

  double                largest_sphere_inside_mesh();
  double                mu_surface_centroid();
  double                mu_surface_sphere();
};