

#include "definitions.h"

class Shell {
public:
  Shell(std::vector<std::vector<int>> trs, std::vector<Point3> lspts);

  Mesh                  get_mesh();
  Polyhedron            get_convex_hull();
  std::array<Point3, 8> get_oobb();
  K::Iso_cuboid_3       get_aabb();
  
  bool                  is_closed();
  double                volume();
  double                area();

private:
  std::vector<Point3>           _lspts;
  std::vector<std::vector<int>> _trs;
  Mesh                          _mesh;

};