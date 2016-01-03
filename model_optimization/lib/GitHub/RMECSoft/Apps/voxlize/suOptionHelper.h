#pragma once

#include "config.h"

/*!
 * \file suOptionHelper.h
 * \date 2015/01/20 9:48
 *
 * \author Yuan Yao
 * Contact: fly2mars@gmail.com
 *
 * \brief 命令行参数处理，需要根据特定的命令重新编写Helper class
 *        所有命令行参数存储在类变量中
 * \example
 *        suCMDParser option(argv, argc);
 *        suOptionValue optionValue(&option);
 *        optionValue.getParams();
 *
 * TODO: 兼容“-” 和 “/”
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