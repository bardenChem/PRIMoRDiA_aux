//atom.h

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
#include <cmath>
#include <iostream>

#include "../include/global.h"
#include "../include/atom.h"

using std::cout;
using std::endl;
using std::move;

/**********************************************************/
atom::atom()	:
	xc(0.00)	,
	yc(0.00)	,
	zc(0.00)	,
	pCharge(0.0),
	aMass(1.0)	,
	aNmb(1)		,
	element("H"){
}
/**********************************************************/
atom::~atom(){}
/**********************************************************/
atom::atom(double x			,
		   double y			, 
		   double z			, 
		   std::string type):
		   xc(x)			,
		   yc(y)			,
		   zc(z)			,
		   pCharge(0.0)		,
		   element(type)	{
	
	aMass = get_atom_mass(element);
	aNmb  = get_atomic_number(element);
}
/**********************************************************/
atom::atom(const atom& rhs)	:
	xc(rhs.xc)				,
	yc(rhs.yc)				,
	zc(rhs.zc)				,
	element(rhs.element)	,
	pCharge(rhs.pCharge)	,
	aMass(rhs.aMass)		,
	aNmb(rhs.aNmb)			{
	
}
/**********************************************************/
atom& atom::operator=(const atom& rhs){
	if ( this != &rhs ){
		xc		= rhs.xc;
		yc		= rhs.yc;
		zc		= rhs.zc;
		element	= rhs.element;
		pCharge	= rhs.pCharge;
		aMass	= rhs.aMass;
		aNmb	= rhs.aNmb;
	}
	return *this;
}
/**********************************************************/
atom::atom(atom&& rhs) noexcept :
	xc(rhs.xc)					,
	yc(rhs.yc)					,
	zc(rhs.zc)					,
	element( move(rhs.element) ),
	pCharge(rhs.pCharge)		,
	aMass(rhs.aMass)			,
	aNmb(rhs.aNmb)				{
}
/**********************************************************/
atom& atom::operator=(atom&& rhs) noexcept{
	if ( this != &rhs ){
		xc		= rhs.xc;
		yc		= rhs.yc;
		zc		= rhs.zc;
		element	= move(rhs.element);
		pCharge	= rhs.pCharge;
		aMass	= rhs.aMass;
		aNmb	= rhs.aNmb;
	}
	return *this;
}
/**********************************************************/
void atom::set_element( std::string Type ){
	element = Type;
	aMass 	= get_atom_mass(element);
	aNmb  	= get_atomic_number(element);
}
/**********************************************************/
void atom::set_coord(double x, double y, double z){
	xc = x;
	yc = y;
	zc = z;	
}
/**********************************************************/
double atom::get_distance(const atom& a2){
	double dist = 0.0;
	dist = (xc - a2.xc)*(xc - a2.xc);
	dist+= (yc - a2.yc)*(yc - a2.yc);
	dist+= (zc - a2.zc)*(zc - a2.zc);
	return sqrt(dist);
}
/**********************************************************/
void atom::set_pCharge(double chg){
	pCharge = chg;
}
/**********************************************************/
std::string atom::print(){
	
	std::string p_info = " Printing atom information\n";
	p_info +=  "\tElement name: ";
	p_info += element;
	p_info += "\n\tx axis coordinate: ";
	p_info += std::to_string(xc);
	p_info += "\n\yx axis coordinate: ";
	p_info += std::to_string(yc);
	p_info += "\n\zx axis coordinate: ";
	p_info += std::to_string(zc);
	p_info += "\n\tPartial charge: ";
	p_info += std::to_string(pCharge);
	p_info += "\n\tAtomic mass: ";
	p_info += std::to_string(aMass);
	p_info += "\n\tAtomic number: ";
	p_info += std::to_string(aNmb);
	
	cout << p_info << endl;
	return p_info;
}

///////////////////////////////////////////////////////////
void UnitTest_atom(){
	atom _Atom_A; // default constructor 
	_Atom_A.print();
	atom _Atom_B(0.000,0.000,0.000, "H" ); // info constructor
	_Atom_B.print();
	atom _Atom_C(_Atom_B); // copy constructor
	
	ut_log.input_line("Printing Copied atom: ");
	ut_log.input_line( _Atom_C.print() );
	
	atom _Atom_D = _Atom_C; // operator overloading 
	cout << "Printing atomic copied information with operator overloading:" << endl;
	_Atom_D.print();
	
	atom _Atom_E( move(_Atom_D) ); // move constructor
	cout << "Printing atom after being moved: " << endl;
	_Atom_D.print();
	cout << "Printing atom after receive infor from move : " << endl;
	_Atom_E.print();
	
	atom _Atom_F = move(_Atom_E); // move assign operator oveloading
	cout << "Printing atom after being moved: " << endl;
	_Atom_E.print();
	cout << "Printing atom after receive infor from move : " << endl;
	_Atom_F.print();
	
	cout << "Testing the charge setting " << endl;
	_Atom_F.set_pCharge( -1.00 );
	cout << "Testing coordinates setting " << endl;
	_Atom_F.set_coord(1.000,0.000,0.000);
	cout << "Testing elment setting " << endl;
	_Atom_F.set_element("C");
	cout << "Printing new information" << endl;
	_Atom_F.print();
	
	cout << "Calculate distance with other atom : " << endl;
	double distance = _Atom_F.get_distance(_Atom_B);
	cout << "Calculated distance: " << distance << endl;
	cout << "Finishing Unit test for class 'atom'" << endl;
}
/**********************************************************/