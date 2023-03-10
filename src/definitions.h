#ifndef DEFS_H
#define DEFS_H

// CGAL kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/squared_distance_3.h> 

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

//-- for mark_domain()
struct FaceInfo2
{
  FaceInfo2() {}
  int nesting_level;
  bool in_domain() {
    return nesting_level % 2 == 1;
  }
};  

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2                  Point2;
typedef K::Point_3                  Point3;
typedef K::Tetrahedron_3            Tetrahedron;
typedef CGAL::Polyhedron_3<K>       Polyhedron;
typedef CGAL::Polygon_2<K>          Polygon2;
typedef K::Plane_3                  Plane;
typedef CGAL::Surface_mesh<Point3>  Mesh;

// CGAL typedefs
typedef CGAL::Triangulation_vertex_base_with_id_2 <K>             Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K>   Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb>       Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_intersections_tag                             Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CT;

typedef boost::graph_traits<Mesh>::vertex_descriptor        vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor      halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;

typedef CGAL::Min_sphere_of_points_d_traits_3<K,double>     MSPT;
typedef CGAL::Min_sphere_of_spheres_d<MSPT>                 Min_sphere;

#endif
