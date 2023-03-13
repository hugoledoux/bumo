

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
  double                circumference();
  double                cohesion();
  double                convexity();
  double                cubeness();
  double                cuboidindex();
  double                depth();
  double                dispersion();
  // double                exchange();   TODO? NEEDS Nef for intersection of 3D mesh
  double                fractality();
  double                girth();
  double                hemisphericality();
  double                proximity();
  double                range();
  double                rectangularity();
  double                roughness();
  double                spin();
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

  double                avg_dist_samples_surface_centroid();
  double                avg_dist_samples_volume_centroid();
  double                avg_dist_samples_volume_surface();
  double                avg_sq_dist_samples_volume_centroid();
  double                largest_sphere_inside_mesh();
  double                avg_dist_samples_surface_radius_sphere();
};