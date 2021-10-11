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
std::ostream& operator<<(std::ostream& out, const atom& obj){
	
	out << "Outputting atomic information \n" 
		<< "\tElement name: " << obj.element
		<< "\n\tx axis coordinate: " << obj.xc
		<< "\n\ty axis coordinate: " << obj.yc
		<< "\n\tz axis coordinate: " << obj.zc
		<< "\n\tPartial charge: " << obj.pCharge
		<< "\n\tAtomic mass: " << obj.aMass
		<< "\n\tAtomic number: " << obj.aNmb;
	
	return out;
}

/**********************************************************/
void atom::print(){
	cout << *this << endl;
}
///////////////////////////////////////////////////////////
void UnitTest_atom(){
	ut_log.input_line("Starting 'atom' class unit test!");
	ut_log.input_line("================================");
	atom _Atom_A; // default constructor 
	ut_log.data << _Atom_A << endl;
	atom _Atom_B(0.000,0.000,0.000, "H" ); // info constructor
	ut_log.data << _Atom_B << endl;
	atom _Atom_C(_Atom_B); // copy constructor
	
	ut_log.input_line("Printing Copied atom: ");
	ut_log.data << _Atom_C << endl;
	
	atom _Atom_D = _Atom_C; // operator overloading 
	
	ut_log.input_line( "Printing atomic copied information with operator overloading:");
	ut_log.data << _Atom_D << endl;
	
	atom _Atom_E( move(_Atom_D) ); // move constructor
	
	ut_log.input_line( "Printing atom after being moved: ");
	ut_log.data <<_Atom_D << endl;
	ut_log.input_line( "Printing atom after receive infor from move : ");
	ut_log.data <<_Atom_E << endl;
	
	atom _Atom_F = move(_Atom_E); // move assign operator oveloading
	ut_log.input_line("Printing atom after being moved: ");
	ut_log.data <<_Atom_E  << endl;
	ut_log.input_line("Printing atom after receive infor from move : ");
	ut_log.data <<_Atom_F  << endl;
	
	ut_log.input_line("Testing the charge setting ");
	_Atom_F.set_pCharge( -1.00 );
	ut_log.input_line("Testing coordinates setting ");
	_Atom_F.set_coord(1.000,0.000,0.000);
	ut_log.input_line("Testing elment setting ");
	_Atom_F.set_element("C");
	ut_log.input_line("Printing new information");
	ut_log.data <<_Atom_F  << endl;
	
	ut_log.input_line("Calculate distance with other atom : ");
	double distance = _Atom_F.get_distance(_Atom_B);
	ut_log.input_line("Calculated distance: "); 
	ut_log.data << distance << endl;
	ut_log.input_line("Finishing Unit test for class 'atom'");
}
/**********************************************************/