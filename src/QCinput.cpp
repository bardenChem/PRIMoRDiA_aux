//QCinput.cpp

/*********************************************************************/
/* This source code file is part of LQQCMMtools software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
/*********************************************************************/

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <experimental/filesystem>

#include "../include/global.h"
#include "../include/system.h" 
#include "../include/geometry.h"
#include "../include/gamessInput.h"
#include "../include/orcaInput.h"
#include "../include/mopacInput.h"
#include "../include/QCinput.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;
namespace fs = std::experimental::filesystem;

/************************************************************/
QCPinput::QCPinput()	:
	program(0)			,
	base_multi(1)		,
	base_charge(0)		,
	QMmethod("HF")		,
	base_basis("3-21G")	,
	geo_ext("xyz")		{
		
}
/************************************************************/
QCPinput::QCPinput(string _folder	,
				   string _geo_ext	, 
				   string _basis	, 
				   string _method)	{
					   

}
/************************************************************/
QCPinput::~QCPinput(){}
/************************************************************/
void QCPinput::make_input_from_folder(	int QCP				,
										unsigned int bs_mlt	,
										int bs_chg)			{
											
}
/************************************************************/								
void QCPinput::make_input_from_folder_FD(int QCP			,
										 unsigned int bs_mlt,
										 int bs_chg			,
										 int chg_diff		){
											 
}
//////////////////////////////////////////////////////////////