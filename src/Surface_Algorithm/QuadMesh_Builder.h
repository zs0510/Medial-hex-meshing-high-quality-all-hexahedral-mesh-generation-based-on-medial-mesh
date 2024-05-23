#pragma once
#include "Kernel/Mesh.h"

namespace QuadMesh_Builder {

	void build_quad_mesh_from_hex_mesh(MeshKernel::VolumeMesh& hexmesh, MeshKernel::SurfaceMesh& quadmesh) {

		quadmesh.destory();

		std::unordered_map<VH, VH> hex2quad;
		for (auto& vp : hexmesh.allvertices()) {
			if (hexmesh.isOnBoundary(vp.first)) {
				VH qvh = quadmesh.AddVertex(vp.second);
				hex2quad[vp.first] = qvh;
			}
		}

		for (auto& fp : hexmesh.allfaces()) {
			if (hexmesh.isOnBoundary(fp.first)) {
				FH hex_fh = fp.first;
				Vex face_center = hexmesh.getFaceCenter(hex_fh);
				CH ch = *(hexmesh.NeighborCh(fp.first).begin());
				Vec face_normal = hexmesh.getFaceNormal(hex_fh);
				Vex cell_center = hexmesh.getCellCenter(ch);
				Vec dir_out = (face_center - cell_center).normalized();
				auto vhs = fp.second.getVertexHandle();
				if (dir_out.dot(face_normal) < 0) std::reverse(vhs.begin(), vhs.end());
				for (auto& vh : vhs) {
					vh = hex2quad[vh];
				}
				quadmesh.AddFace(vhs);
			}
		}

	}

};