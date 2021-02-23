//system.cpp

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

#include <"../include/global.h"
#include <"../include/system.h"

using std::vector;
using std::string;
using std::move;

/**********************************************************/
atom::atom()	:
	xc(0.00)	,
	yc(0.00)	,
	zc(0.00)	,
	pCharge(0.0),
	aMass(1.0)	,
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
		   element(type)	{
	
	aMass = get_atom_mass(element);
}
/**********************************************************/
atom::atom(const atom& rhs)	:
	xc(rhs.xc)				,
	yc(rhs.yc)				,
	zc(rhs.zc)				,
	element(rhs.element)	,
	pCharge(rhs.pCharge)	,
	aMass(rhs.aMass)		{
	
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
	}
	return *this;
}
/**********************************************************/
atom::atom(atom&& rhs) noexcept:
	xc(rhs.xc)					,
	yc(rhs.yc)					,
	zc(rhs.zc)					,
	element( move(rhs.element) ),
	pCharge(rhs.pCharge)		,
	aMass(rhs.aMass)			{	
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
	}
	return *this;
}
/**********************************************************/
void atom::set_element( std::string type ){
	element = type;
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
////////////////////////////////////////////////////////////
system::system()	:
	nAtoms(0)		,
	nElectrons(0)	,
	name("nonamed")	,
	type("none")	,
	fCharge(0.0)	,
	{
	
}
/**********************************************************/
system::~system(){
	
}
/**********************************************************/
system::system(const system& rhs){
	
}
/*********************************************************/
system::system(const system& rhs){
	
}
/*********************************************************/
system& system::operator=(const system& rhs){
	
}
/*********************************************************/
system::system(system&& rhs) noexcept:

{
	
}
/*********************************************************/
system& system::operator=(system&& rhs) noexcept:

{
	
}
/*********************************************************/
void system::add_atom(atom a){
	
}
/*********************************************************/
void system::add_atom(double x,double y, double z, std::string el){
	
}
/*********************************************************/
void system::remove_atom(unsigned int i){
	
}
////////////////////////////////////////////////////////////
pdbAtom::pdbAtom(){
	
}
/*********************************************************/
pdbAtom::~pdbAtom(){
	
}
/*********************************************************/
pdbAtom::pdbAtom(const pdbAtom& rhs){
	
}
/*********************************************************/
pdbAtom& pdbAtom::operator=(const pdbAtom& rhs){
	
}
/*********************************************************/
pdbAtom::pdbAtom(pdbAtom&& rhs) noexcept:
{
	
}
/*********************************************************/
pdbAtom& pdbAtom::operator=(pdbAtom&& rhs) noexcept:
{
	
}
/*********************************************************/
bool  pdbAtom::operator!=(const pdbAtom& rhs){
	
}
////////////////////////////////////////////////////////////
residue::residue(){
	
}
/*********************************************************/
residue::~residue(){
	
}
/*********************************************************/
residue::residue(const residue& rhs){
	
}
/*********************************************************/
residue& residue::operator=(const residue& rhs){
	
}
/*********************************************************/
residue::residue( residue&& rhs) noexcept{
	
}
/*********************************************************/
residue& residue::operator=( residue&& rhs) noexcept{
	
}
/*********************************************************/
void residue::set_charge(){
	
}
////////////////////////////////////////////////////////////
pdbModel::pdbModel(){
	
}
/*********************************************************/
pdbModel::pdbModel(std::vector<residue> residues){
	
}
/*********************************************************/
pdbModel::~pdbModel(){}
/*********************************************************/
pdbModel::pdbModel(const pdbModel& rhs){
	
}
/*********************************************************/
pdbModel& pdbModel::operator=(const pdbModel& rhs){
	
}
/*********************************************************/
pdbModel::pdbModel(pdbModel&& rhs) noexcept:
{
	
}
/*********************************************************/
pdbModel& pdbModel::operator=(pdbModel&& rhs) noexcept:
{
	
}
/*********************************************************/
void pdbModel::write_model(std::string out_name){
	
}
/*********************************************************/
void pdbModel::prune_atoms(){
	
}
/*********************************************************/
void pdbModel::remove_waters(){
	
}
/*********************************************************/
void pdbModel::remove_waters(double radius){
	
}
/*********************************************************/
void pdbModel::remove_ions(){
	
}
/*********************************************************/
void pdbModel::remove_residue(unsigned int i){
	
}
/*********************************************************/
void pdbModel::update_residues(){
	
}
/*********************************************************/
void pdbModel::split_complex(std::string mol){
	
}
/*********************************************************/
void pdbModel::built_complex(const char* pdb_mol){
	
}
//////////////////////////////////////////////////////////