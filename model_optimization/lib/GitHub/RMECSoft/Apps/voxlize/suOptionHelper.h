#pragma once

#include "config.h"

/*!
 * \file suOptionHelper.h
 * \date 2015/01/20 9:48
 *
 * \author Yuan Yao
 * Contact: fly2mars@gmail.com
 *
 * \brief �����в���������Ҫ�����ض����������±�дHelper class
 *        ���������в����洢���������
 * \example
 *        suCMDParser option(argv, argc);
 *        suOptionValue optionValue(&option);
 *        optionValue.getParams();
 *
 * TODO: ���ݡ�-�� �� ��/��
 *
 * \note
**/
enum gChoice{ eByTreeDeep = 0, eByUnitSize = 1 };
enum gFunction{eCOMPUTE_VOLUME, eVOXELIZE};
class suOptionValue{
public:
	suOptionValue(suCMDParser *pCMD) :nTreeDeep(0), nUnitSize(0), eDevideType(eByTreeDeep), pCMDParser_(NULL)
	{
		pCMDParser_ = pCMD;
	}

	/*Check all parameter is ready*/
	bool isReady();
	bool getParams();
	
public:
	std::string sInput;
	std::string sOutput;

	int nTreeDeep;
	int nUnitSize;

	
	gChoice eDevideType;
	gFunction eFunction;

	std::ostringstream errMsg;

private:
	suCMDParser *pCMDParser_;

};