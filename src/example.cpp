MeshKernel::VolumeMesh* hexMesh = new MeshKernel::VolumeMesh();

// 1. Generate hex-mesh for face skeleton
QMessageBox::StandardButton ret1 = QMessageBox::question(this, " Question ", " Is non-manifold ? ", QMessageBox::Ok | QMessageBox::Cancel);
auto start_app = std::chrono::high_resolution_clock::now();
if (ret1 == QMessageBox::Cancel) {
	HexMeshing_FaceSkeleton app0(*skeletal_mesh, *hexMesh);
	app0.hexmeshing();
	updateDrawData();
} else {
	HexMeshing_FaceSkeleton_Nonmanifold app0(*skeletal_mesh, *hexMesh);
	app0.hexmeshing();
	updateDrawData();
}

VH algorithm_insert(-1);
// 2. Generate branching elements
HexMeshing_Branch_Processor app1(*skeletal_mesh, *hexMesh);
app1.optimize(algorithm_insert, true);
updateDrawData();
auto point1 = std::chrono::high_resolution_clock::now();
QMessageBox::StandardButton ret = QMessageBox::question(this, " Question ", " Continue sweep hex ? ", QMessageBox::Ok | QMessageBox::Cancel);
if (ret == QMessageBox::Cancel) return;
auto point2 = std::chrono::high_resolution_clock::now();

// 3. Generate hex-mesh for curve skeleton
HexMeshing_ConformingTessellation app2(*skeletal_mesh, *hexMesh);
app2.conforming_split(algorithm_insert);

auto end_app = std::chrono::high_resolution_clock::now();

std::chrono::duration<double> duration = end_app - start_app;
std::chrono::duration<double> duration_req = point2 - point1;

std::cout << "[HexMeshing]: cost " << (duration.count() - duration_req.count()) << "s\n";