#include "config.h"
#include <iostream>
#include <string>
#include <map>
#include "../voxlize/suMesh.h"
#include "../voxlize/suHybridMesh.h"

using namespace std;


UTFUNC(readmodel)
{
	//1.volumize the mesh
	suHybridMesh& mesh;
	mesh.LoadMeshFile("");
	mesh.PartitionSpace(4);

}
int _tmain(int argc, _TCHAR* argv[])
{
	
	bool state = suUnitTest::gOnly().run ();
	suUnitTest::gOnly().dumpResults (std::cout);
	

	return state;
}
