/*
 *   Author:   $Author: hexc $
 *   Date:     $Date: 2014/08/06 15:59:09 $
 *   Revision: $Revision: 1.2 $
 *   Tag Name: $Name:  $
 *   
 *   $Log: SpenvisCSV.cc,v $
 *   Revision 1.2  2014/08/06 15:59:09  hexc
 *   Compiled for Geant4.9.6
 *
 *   Revision 1.1.1.1  2014/06/02 18:39:07  hexc
 *   magnetocosmics project
 *
 *   Revision 1.1.1.1  2011/07/14 13:39:40  hexc
 *   Cosmic ray air shower simulation software
 *
 * 
 *   Modification  2004/12/04   	ldesorgher
 *		the map of meta variables is replaced by a vector of meta variables 
 *		some cleaning to avoid warnings in G4	   
 *   Revision 1.2  2004/09/16 12:45:29  hevans
 *   added reading and parsing of SpenvisCVS files
 *
 *   Revision 1.1.1.1  2004/09/08 09:23:45  hevans
 *   First Release
 * 
 *
 */ 
#include "SpenvisCSV.hh"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

SpenvisCSV::SpenvisCSV ()
{

  nvar_cols = 0;
}

SpenvisCSV::~SpenvisCSV()
{
}

void SpenvisCSV::AddComment( STRING val) {
  Comments c;
  c.val = val;

  comments.push_back( c);

}

void SpenvisCSV::AddAnnotation( STRING val) {

  annotations.push_back( val);
}

void SpenvisCSV::AddMetaVariable   ( STRING name, int dims, float *vals, STRING format, STRING unit){
  
  MetaVar mv;
  int i;

  mv.name = name;
  mv.dims = dims;
  mv.format = format;
  mv.unit = unit;
  for (i=0; i<dims; i++) 
    mv.fvals.push_back( vals[i]);

//Modified by ldesorgher
/*  metaVar[name] = mv;
*/
 unsigned int index = metaVar.size();
 for (unsigned int i=0;i<metaVar.size();i++){
 	if (metaVar[i].name == name ){
		index = i;
		metaVar[i]=mv;
		i=metaVar.size();
	}
 
 }
 if (index == metaVar.size()) metaVar.push_back(mv);
  
}

void SpenvisCSV::AddMetaVariable   ( STRING name, std::vector <double> vals, STRING format, STRING unit){
  MetaVar mv;
 
  mv.name = name;
  mv.dims = vals.size();
  mv.format = format;
  mv.fvals  = vals;
  mv.unit = unit;

//Modified by ldesorgher
//  metaVar[name] = mv;
 
 unsigned int index = metaVar.size();
 for (unsigned int i=0;i<metaVar.size();i++){
 	if (metaVar[i].name == name) {
		index = i;
		metaVar[i]=mv;
		i=metaVar.size();
	}
 
 }
 if (index == metaVar.size()) metaVar.push_back(mv);
  
  
}

void SpenvisCSV::AddMetaVariableStr( STRING name, STRING val){
  MetaVar mv;
  mv.name = name;
  mv.dims = -1;
  mv.svals = val;

//Modified by ldesorgher
//  metaVar[name] = mv;

 unsigned int index = metaVar.size();
 for (unsigned int i=0;i<metaVar.size();i++){
 	if (metaVar[i].name == name) {
		index = i;
		metaVar[i]=mv;
		i=metaVar.size();
	}
 
 }
 if (index == metaVar.size()) metaVar.push_back(mv);
   

}

void SpenvisCSV::AddVariable( STRING name, STRING unit, int dims, STRING description, STRING fmt)
{
  Vars  var;
  int   i;

  var.name = name;
  var.unit = unit;
  var.dims = dims;
  var.descr = description;
  var.format= fmt;

  vars.push_back( var);
  nvar_cols+= dims;

  for (i=0; i<dims; i++)
    col_formats.push_back( fmt);
}

void SpenvisCSV::AddDataRow( double *rec){
  int i;
  Vals v;

  for (i=0; i< nvar_cols; i++) 
    v.values.push_back( rec[i]);
  vals.push_back(v);
  
}

void SpenvisCSV::AddDataRow( std::vector <double> rec){
  Vals v;

  v.values = rec;
  if ( int (rec.size()) < nvar_cols) 
    std::cout << "*** AddDataRow: Number of values in input vector less than expected\n";
  vals.push_back( v);
}


STRING SpenvisCSV::GetComment( int i){
  
  return comments[i].val;
}

STRING SpenvisCSV::GetMetaVariableName( int i) {

//Modified by ldesorgher
/*  typedef map<string, MetaVar>::const_iterator CI;
  CI p;

  for ( p=metaVar.begin(); (p != metaVar.end() ) && ( i > 0); p++ ){
    i--;
  }
  return p->first;
*/
  if (i>0 && i<int(metaVar.size())) return  metaVar[i].name;
  else return ""; 
}

STRING SpenvisCSV::GetMetaVariableLine( int i) {

//Modified by ldesorgher
/*  MetaVar mv;
  typedef map<string, MetaVar>::const_iterator CI;
  CI p;
  char sdim[3], sdbl[15];
  STRING outstr;

  int j = 0;

  for ( p = metaVar.begin(); (p != metaVar.end()) && (i > 0); ++p ) {
    i--;
  }
  mv = p->second;
*/
  std::cout<<"MV1"<<std::endl;	
  MetaVar mv =metaVar[i];
  char sdim[3], sdbl[15];
  STRING outstr,outstr1;
  sprintf( sdim, "%3d", mv.dims);
  outstr = "'" + mv.name + "'," + sdim;
  std::cout<<"MV2"<<std::endl;	
  if ( mv.dims < 0) {
    outstr = outstr +",'"+ mv.svals+"'";
  } else {
    for (i=0; i<mv.dims; i++) {
      sprintf( sdbl, mv.format.c_str(), mv.fvals[i]);
      outstr = outstr + ", " + sdbl;
    }
    if (mv.unit != "" && mv.dims >0)
    		outstr = outstr + ",'" + mv.unit + "'";
  }
  std::cout<<"MV3"<<std::endl;
  std::cout<<outstr<<std::endl;
  outstr1=outstr;
  return ""+outstr;
}

STRING SpenvisCSV::GetVariableLine( int i) {
  Vars var;
  char sdim[6];
  STRING outs;

  var = vars[i];
  sprintf( sdim, "%5d", var.dims);

  outs = "'" + var.name + "','" + var.unit + "'," + sdim + ",'" + var.descr +"'";

  return outs;
}

STRING SpenvisCSV::GetVariableName( int i) {
  if ( i < int(vars.size()) ) {
    return vars[i].name;
  }
  else return "";
}

STRING SpenvisCSV::GetDataLine( int i) {
  Vals v;
  STRING outs ="";
  char str[80];
  int j;

  v = vals[i];
  sprintf( str, col_formats[0].c_str(), v.values[0]);
  outs = str;

  for (j=1; j< nvar_cols; j++) {
    sprintf(str, col_formats[j].c_str(), v.values[j]);
    outs = outs + ", " + str;
  }
  return outs;
}

std::vector <double> SpenvisCSV::GetDataRecord( int recordNum) {
  if ( recordNum < int(vals.size()) ){
    return vals[recordNum].values;
  }
  std::vector<double> a;
  a.clear();
  return a; 
}

int SpenvisCSV::GetVariableNumber( STRING name) {
  int i=0;
  //int retval = -1;

  while ( (i < int(vars.size())) && (vars[i].name != name)) i++;
  if ( i >= int(vars.size()) ) return -1;
  if (name == vars[i].name) return i;
  return -1;
}

std::vector <double> SpenvisCSV::GetDataRecord( int recordNum, int varNum) {
  int i=0, startcol=0;
  std::vector <double> data;

  if ((varNum < 0 ) || (varNum >= int(vars.size()))) return data;

  // Get the column indices for the required variable
  startcol = 0;
  for (i=0; i< varNum; i++) {
    startcol += vars[i].dims;
  }

  // Found the variable, now get the data
  for (int j=startcol; j < (startcol+vars[varNum].dims); j++) {
    data.push_back( vals[ recordNum].values[ j]);
  }
  return data;

}

std::vector <double> SpenvisCSV::GetDataRecord( STRING variable, int recordNum){
  return GetDataRecord( recordNum, GetVariableNumber( variable));
}

void SpenvisCSV::ResetBlock() {

  comments.erase( comments.begin(), comments.end());
  annotations.erase( annotations.begin(), annotations.end());
  metaVar.erase( metaVar.begin(), metaVar.end());
  vars.erase( vars.begin(), vars.end());
  col_formats.erase( col_formats.begin(), col_formats.end());
  EraseData();

  nvar_cols = 0;

  headerLine.nHeader =0;
  headerLine.nText = 0;
  headerLine.nMeta = 0;
  headerLine.nFut  = 0;
  headerLine.nVar  = 0;
  headerLine.nCol  = 0;
  headerLine.nBody = 0;
  headerLine.npart = 0;
}

void SpenvisCSV::ParseBlock( std::vector <STRING> ins){
  int i=0;
  int nrecs=0;

  ResetBlock();

  // At this point, the headerLine record structure has the various line
  // counts for the blocks, and the various GetNum...() methods will return 0.

  ParseBlockHeaderLine( ins[i++]);

  for (int j=0; j<headerLine.nText; j++)
    ParseCommentLine( ins[i++]);

  for (int j=0; j<headerLine.nMeta; j++)
    ParseMetaDataLine( ins[i++]);

  for (int j=0; j<headerLine.nFut; j++)
    ParseAnnotationLine( ins[i++]);

  for (int j=0; j<headerLine.nVar; j++)
    ParseVariableLine( ins[i++]);

  // If the number of points is specified, then read that many points.
  if (headerLine.nBody > 0 ) {
    nrecs = headerLine.nBody;
    for (int j=0; j<nrecs; j++) 
      ParseDataRow( ins[i++]);
  }

  // In this case, then the number of lines is unknown. Keep reading until the
  // end of the string vector, or until an apostrophy is encountered.
  if (headerLine.nBody < 0) {
    while ( (i < int(ins.size())) && (ins[i].find('\'') == ins[i].npos) ) 
      ParseDataRow( ins[i++] );
    headerLine.nBody = vals.size();
  }
}

void SpenvisCSV::ParseBlockHeaderLine( STRING ins) {
  char delim = ',';
  std::vector <STRING> args;
  int i;

  args = strsplit( ins, delim);

  if ( args[0] != "'*'" ) {
    std::cout << "ParseDataLine: ERROR, header doesn't start with '*'\n";
  }
  i=1;
  headerLine.nHeader = atoi( args[i++].c_str());
  headerLine.nText   = atoi( args[i++].c_str());
  headerLine.nMeta   = atoi( args[i++].c_str());
  headerLine.nFut    = atoi( args[i++].c_str());
  headerLine.nVar    = atoi( args[i++].c_str());
  headerLine.nCol    = atoi( args[i++].c_str());
  headerLine.nBody   = atoi( args[i++].c_str());
  headerLine.npart   = atoi( args[i++].c_str());
}

void SpenvisCSV::ParseBlockHeaderLine( STRING ins, HeaderLine *hl) {

  ParseBlockHeaderLine( ins);
  hl->nHeader = headerLine.nHeader ;
  hl->nText   = headerLine.nText   ;
  hl->nMeta   = headerLine.nMeta   ;
  hl->nFut    = headerLine.nFut    ;
  hl->nVar    = headerLine.nVar    ;
  hl->nCol    = headerLine.nCol    ;
  hl->nBody   = headerLine.nBody   ;
  hl->npart   = headerLine.npart   ;
}

void SpenvisCSV::ParseMetaDataLine( STRING ins){
  std::vector <STRING> args;
  char   apostrophe = '\'';
  std::vector <double> fvals;
  int dims;

  args = strsplit( ins, ',');

  dims = atoi( args[1].c_str());

  // If it is a string, then just store the string, otherwise loop over the 
  // remaining arguments, converting them to floating point numbers.
  if ( dims < 0) {
    AddMetaVariableStr( strtrim( args[0], apostrophe),
			strtrim( args[2], apostrophe) 
			); 
  }  
  else {

//Modified by ldesorgher 
    if ( int(args.size())-2 < dims) dims = args.size() - 2;
    for (int i=2; i < 2+dims; i++) {
      fvals.push_back( atof( args[i].c_str()) );
    }
    STRING unit="";
    if ( dims== int(args.size())-3) unit = strtrim( args.back(), apostrophe);
    AddMetaVariable( strtrim( args[0], apostrophe),
		     fvals, "%e",unit);
  }
}


void SpenvisCSV::ParseVariableLine( STRING ins) {
  std::vector <STRING> args;
  char   apostrophe = '\'';

  args = strsplit( ins, ',');

  AddVariable( strtrim(args[0], apostrophe), 
	       strtrim(args[1], apostrophe), 
	       atoi( args[2].c_str()), 
	       strtrim(args[3], apostrophe),
	       "%e"
	       );
}

void SpenvisCSV::ParseCommentLine ( STRING ins) {
  char   apostrophe = '\'';
  AddComment( strtrim(ins, apostrophe));
}
void SpenvisCSV::ParseAnnotationLine( STRING ins) {
  char   apostrophe = '\'';
  AddAnnotation( strtrim(ins, apostrophe));
}

void SpenvisCSV::ParseDataRow( STRING ins) {
  std::vector <double> f_vals;
  std::vector <STRING> s_vals;

  s_vals = strsplit( ins, ',');
  for (int i=0; i<int(s_vals.size()); i++){
    f_vals.push_back( atof( s_vals[i].c_str()) );
  }
  AddDataRow( f_vals);
}

STRING SpenvisCSV::GetBlockHeaderLine( int nblocks) {
  int nh, nc, nm, na, nv, nd, n1;
  char c[80];
  STRING outs;

  // Numb comments
  nc = comments.size();
  nm = metaVar.size();
  na = annotations.size();
  nv = vars.size();
  nd = nvar_cols;
  n1 = vals.size();

  nh = 1 + nc + nm + na + nv;
  
  sprintf( c, "'*', %d, %d, %d, %d, %d, %d, %d, %d",
	   nh, nc, nm, na, nv, nd, n1, nblocks);

  //std::cout << "GetBlockHeaderLine(): " << c << "\n";

  outs = c;
  return outs;
}


void SpenvisCSV::OutputCSV( ofstream &outs, int nFollowingBlocks ){
  STRING output;
  int i;
  std::cout <<"OCSV1"<<std::endl;
  outs << GetBlockHeaderLine( nFollowingBlocks) << "\n";
  std::cout <<"OCSV2"<<std::endl;
  for ( i=0; i< int(comments.size()); i++) 
    outs <<"'"<< comments[i].val<< "'" << "\n";
  std::cout <<"OCSV3"<<std::endl;
  for ( i=0; i< int(metaVar.size()); i++){
    	std::cout <<i<<std::endl;
	outs << GetMetaVariableLine(i)<<"\n";
  	std::cout <<"OCSV"<<std::endl;
  }	
  std::cout <<"OCSV4"<<std::endl;
  for ( i=0; i< int(annotations.size()); i++) 
    outs <<"'"<< annotations[i] << "'" << "\n";
  std::cout <<"OCSV5"<<std::endl;  
  for ( i=0; i< int(vars.size()); i++) 
    outs << GetVariableLine(i) <<"\n";
  std::cout <<"OCSV6"<<std::endl;
  for ( i=0; i< int(vals.size()); i++)
    outs << GetDataLine( i) << "\n";
  std::cout <<"OCSV7"<<std::endl;
  outs << "'End of Block'\n" ;
}

void SpenvisCSV::OutputCSV         (int nFollowingBlocks){

  STRING output;
  int i;

  std::cout << GetBlockHeaderLine( nFollowingBlocks) << "\n";
  
  for ( i=0; i< int(comments.size()); i++) 
    std::cout <<"'"<< comments[i].val<< "'" << "\n";

  for ( i=0; i< int(metaVar.size()); i++)
    std::cout << GetMetaVariableLine(i)<<"\n";

  for ( i=0; i< int(annotations.size()); i++) 
    std::cout <<"'"<< annotations[i] << "'" << "\n";
    
  for ( i=0; i< int(vars.size()); i++) 
    std::cout << GetVariableLine(i) <<"\n";

  for ( i=0; i< int(vals.size()); i++)
    std::cout << GetDataLine( i) << "\n";

  std::cout << "'End of Block'\n" ;
}


std::vector <STRING> SpenvisCSV::strsplit( STRING ins, char delim) {
  int i, currStart, currEnd;
  std::vector <STRING> outs;

  i = 0;
  while ( i < int(ins.size()) ) {
    // Loop until a delimiter or the end of string is found.
    currStart = i;
    while ( i < int(ins.size()) && ins[i] != delim) i++;
    currEnd = i - 1;
    outs.push_back( ins.substr( currStart, currEnd-currStart+1) );
    i++;
  }
  return outs;
}

STRING SpenvisCSV::strremove( STRING ins, char target) {
  STRING outs;
  int i;

  outs = ins;

  for (i=0; i< int(outs.size()); i++) {
    if ( outs[i] == target) outs.erase(i--,1);
  }

  return outs;
}

STRING SpenvisCSV::strtrim( STRING ins, char target) {
  STRING outs;
  //int endchar;

  outs = ins;
  if ( outs.find( target) == 0)           
    outs.replace(  0, 1, "");
  if ( outs.rfind( target) == (outs.size()-1) ) 
    outs.replace( outs.rfind( target), 1, "");

  return outs;
}
