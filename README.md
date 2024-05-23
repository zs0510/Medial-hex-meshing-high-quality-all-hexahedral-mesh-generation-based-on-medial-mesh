# Medial hex‑meshing: high‑quality all‑hexahedral mesh generation based on medial mesh

This repository contains the data and source code of the paper [Medial hex-meshing: high-quality all-hexahedral mesh generation based on medial mesh](https://link.springer.com/article/10.1007/s00366-023-01925-5).

![pipeline](./figures/pipeline.png)

In this paper, based on the curve/surface-mixed medial axis representation, we propose a new high-quality all-hex mesh generation method. Given an input watertight model, we first compute the corresponding medial mesh via the medial axis transform simplification method. Then, we build the all-hexahedral layer on the surface skeleton of the medial mesh via cross-field-guided quad-meshing and extrusion, and construct the hex-mesh elements for the curve skeleton of the medial mesh via sweeping approach. Based on the topological and geometrical information of the medial mesh, the initial hexahedral mesh can be obtained. Furthermore, with the iterative volumetric subdivision fitting approach, the hexahedral mesh is fitted to the input model. Finally, padding refinement and mesh optimization method are used to improve the element quality. Compared with existing hexahedral mesh generation methods, our method can generate all-hex meshes with simpler singular structure, better element quality, and smaller element numbers.

# Erratum

1. In Fig.14, the figures for *(a) Non-manifold medial mesh* and *(b) Manifold medial mesh* need to be interchanged.

# Citation
```angular2html
@article{zhang2024medial,
  title={Medial hex-meshing: high-quality all-hexahedral mesh generation based on medial mesh},
  author={Zhang, Sheng and Xu, Gang and Wu, Haiyan and Gu, Renshu and Qi, Long and Pang, Yufei},
  journal={Engineering with Computers},
  pages={1--21},
  year={2024},
  publisher={Springer}
}
```
