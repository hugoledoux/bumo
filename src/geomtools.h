#ifndef __geomtools__
#define __geomtools__

#include "definitions.h"


Polyhedron            convex_hull(const std::vector<Point3>& lspts);
Plane                 get_best_fitted_plane(const std::vector<Point3> &lspts);
double                mu(std::vector<Point3>& shellpts, std::vector<std::vector<int>>& trs, const std::vector<Point3>& lspts);
double                area_shell(std::vector<std::vector<int>>& trs, const std::vector<Point3>& lspts);
double                volume_shell(std::vector<std::vector<int>>& trs, const std::vector<Point3>& lspts);
double                oobb_volume(std::array<Point3, 8> oobbpts);
K::Iso_cuboid_3       aabb(const std::vector<Point3>& lspts);
std::array<Point3, 8> oobb(const std::vector<Point3>& lspts);

void                  mark_domains(CT& ct);
void                  mark_domains(CT& ct, CT::Face_handle start, int index, std::list<CT::Edge>& border);
std::vector<std::vector<int>>
construct_ct_one_face(const std::vector<std::vector<int>>& lsRings, 
                      const std::vector<Point3>& lspts);

#endif 