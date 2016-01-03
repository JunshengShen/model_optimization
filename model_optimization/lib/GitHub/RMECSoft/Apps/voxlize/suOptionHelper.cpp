#include "config.h"
#include "suOptionHelper.h"


std::string getHelp()
{
	std::ostringstream msg;
	msg << "���ܣ��������뵥Ԫ�ߴ����ػ����루stl���ļ�" << std::endl;
	msg << "������ʽ��" << std::endl;
	msg << "voxel.exe  [-command] [-D][deep value] [-L] [Unit size] [-I][Input filename] [-O][Output filename]" << std::endl;

	msg << "==========================================================" << std::endl;
	msg << "1.ʵ���� ����octree��Ȼ��ֿռ�" << std::endl;
	msg << "         voxel -ByTreeDeep  -D 6 -I r:/input.stl -O output.vtk" << std::endl;
	msg << "2.ʵ���� ���ݻ��ֵ�Ԫ�ߴ绮�ֿռ�" << std::endl;
	msg << "         RF -ByUnitSize -l 0.5 -I r:/input.stl -O output.vtk" << std::endl;
	msg << "3.ʵ���� �������(mm^2)" << std::endl;
	msg << "         RF -V" << std::endl;
	msg << "4.ʵ���� ��ӡ����" << std::endl;
	msg << "         RF -h" << std::endl;
	msg << "==========================================================" << std::endl;
	return msg.str();
}

bool suOptionValue::isReady(){
	if (sInput.empty())
	{
		errMsg << "error: need input file! \n";
		return false;
	}

	if (sOutput.empty())
	{
		int lastindex = sInput.find_last_of(".");
		sOutput = sInput.substr(0, lastindex) + ".vtk";
	}

	if (eDevideType == eByTreeDeep){
		if (nTreeDeep < 1){
			errMsg << "error: need tree deep value! \n";
			return false;
		}
	}
	if (eDevideType == eByUnitSize){
		if (nUnitSize < 1){
			errMsg << "error: need tree unit size value! \n";
			return false;
		}
	}

	return true;
}

bool suOptionValue::getParams()
{	
	if (pCMDParser_->findOption("-h") || pCMDParser_->noParams())
	{
		std::cout << getHelp();
		return true;
	}

	if (pCMDParser_->findOption("-v"))
	{
		eFunction = eCOMPUTE_VOLUME;
	}

	//ByTreeDeep
	if (pCMDParser_->findOption("-d"))
	{
		eDevideType = eByTreeDeep;
		nTreeDeep = pCMDParser_->findOptionValue("-d").ConvertToInteger();

		eFunction = eVOXELIZE;
	}

	//ByUnitSize
	if (pCMDParser_->findOption("-l"))
	{
		eDevideType = eByUnitSize;
		nUnitSize = pCMDParser_->findOptionValue("-l").ConvertToInteger();

		eFunction = eVOXELIZE;
	}

	if (pCMDParser_->findOption("-i")){
		sInput = pCMDParser_->findOptionValue("-i").CString();
	}

	if (pCMDParser_->findOption("-o")){
		sOutput = pCMDParser_->findOptionValue("-o").CString();
	}

	if (!isReady())
	{
		std::cout << errMsg.str();
		return false;
	}
	return true;
}
