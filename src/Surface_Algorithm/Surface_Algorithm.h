#pragma once
#ifndef IGAME_SURFACE_ALGORITHM
#define IGAME_SURFACE_ALGORITHM

#include "Surface_Algorithm/Parameterization.h"
#include "Surface_Algorithm/Remeshing.h"
#include "Surface_Algorithm/Subdivision_Loop.h"
#include "Surface_Algorithm/Denoising_BNF.h"
#include "Surface_Algorithm/Simplification_QEM.h"
#include "Surface_Algorithm/Simplification_Powered_QEM.h"
#include "Surface_Algorithm/Deformation.h"
#include "Surface_Algorithm/Analysis_TriMesh.h"
#include "Surface_Algorithm/TriMesh_Generator.h"
#include "Surface_Algorithm/BoundaryVertexMerger.h"
#include "Surface_Algorithm/LSCM_Parameterization.h"
#include "Surface_Algorithm/Subdivision_Butterfly.h"
#include "Surface_Algorithm/QuadDomain_Builder.h"
#include "Surface_Algorithm/QuadMeshSubdivision_CC.h"
#include "Surface_Algorithm/Parameterization_ARAP.h"
#include "Surface_Algorithm/DomainFinder.h"

// OpenMesh
#include "Surface_Algorithm/OpenMesh_Algorithm/OpenMesh_Smoothing.h"
#include "Surface_Algorithm/OpenMesh_Algorithm/OpenMesh_Remeshing.h"
#include "Surface_Algorithm/OpenMesh_Algorithm/OpenMesh_Simplification.h"
#include "Surface_Algorithm/OpenMesh_Algorithm/OpenMesh_Denoising.h"
#include "Surface_Algorithm/OpenMesh_Algorithm/OpenMesh_Subdivision.h"

#endif // !IGAME_SURFACE_ALGORITHM
