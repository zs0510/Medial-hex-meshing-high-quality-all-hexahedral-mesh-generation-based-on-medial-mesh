//
// Created by su on 2022/3/13.
//

#ifndef CROSSFIELD__CROSS_FIELD_H_
#define CROSSFIELD__CROSS_FIELD_H_

#include <vector>
#include <array>
#include <unordered_map>
#include <utility>

#include "basic_types.h"

//output one dir of cross on edge
//singularity:1:vertex,2 faceID*3+edgeId in face,3 faceId

int BuildBackgroundMeshAndGuidingField(const std::vector<Vec3> &points,
                                       const std::vector<ID3> &triangles,
                                       const std::vector<ID2> &lines,
                                       std::vector<Vec3> &global_edge_directions,
                                       std::vector<std::array<double,
                                                              9>> &global_triangle_directions,
                                       std::vector<ID3> &global_singularity_list);

#endif //CROSSFIELD__CROSS_FIELD_H_
