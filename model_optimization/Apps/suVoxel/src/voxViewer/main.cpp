
#include "MeshViewer.h"
#include <iostream>

#ifdef _DEBUG

#pragma comment(lib,"OpenMeshCored.lib")
#pragma comment(lib,"OpenMeshToolsd.lib")

#else

#pragma comment(lib,"OpenMeshCore.lib")
#pragma comment(lib,"OpenMeshTools.lib")

#endif // DEBUG

void help()
{
	std::cout << "\nThis is a debug program to test volumetric assembling";

	std::cout << "\n\nHot keys: \n"
		"\tESC - quit the program\n"
		"\tS - save depth of hand in .PNG format\n"
		"\ts - save original depth and color data into a .PNG format\n"
		"To initialize tracking, select the object with mouse\n";
}


int main(int argc, char **argv)
{
  glutInit(&argc, argv);

  MeshViewer window("Volume Viewer", 512, 512);

  if (argc>1)
    window.open_mesh(argv[1]);

  glutMainLoop();
}
