/*
 *   Author:   $Author: hexc $
 *   Date:     $Date: 2014/08/06 15:59:09 $
 *   Revision: $Revision: 1.2 $
 *   Tag Name: $Name:  $
 *   
 *   $Log: SpenvisCSVCollection.cc,v $
 *   Revision 1.2  2014/08/06 15:59:09  hexc
 *   Compiled for Geant4.9.6
 *
 *   Revision 1.1.1.1  2014/06/02 18:39:07  hexc
 *   magnetocosmics project
 *
 *   Revision 1.1.1.1  2011/07/14 13:39:40  hexc
 *   Cosmic ray air shower simulation software
 *
 *   Revision 1.2  2004/09/16 12:45:29  hevans
 *   added reading and parsing of SpenvisCVS files
 *
 *   Revision 1.1.1.1  2004/09/08 09:23:45  hevans
 *   First Release
 * 
 *
 */

#include "SpenvisCSVCollection.hh"
#include <stdio.h>

// ========================================================================
// ========================================================================


void SpenvisCSVCollection::OutputCollection( ofstream &outs) {
  typedef std::map<std::string, SpenvisCSV>::iterator MI;
  int nfol = collection.size()-1;
  std::cout<<"COL1"<<std::endl;
  for (MI p= collection.begin(); p!=collection.end(); ++p) {
    std::cout<<"COL2"<<std::endl;	
    p->second.OutputCSV( outs, nfol);
    std::cout<<"COL3"<<std::endl;
    nfol--;
  }
  std::cout<<"COL2"<<std::endl;
  outs << "'End of File'\n";
}

void SpenvisCSVCollection::OutputCollection( ) {
  typedef std::map<std::string, SpenvisCSV>::iterator MI;
  int nfol = collection.size()-1;

  for (MI p= collection.begin(); p!=collection.end(); ++p) {
    p->second.OutputCSV( nfol);
    nfol--;
  }
  std::cout << "'End of File'\n";
}

// ========================================================================
// ========================================================================

void SpenvisCSVCollection::ReadCollection( ifstream &ins){
  std::vector <STRING> inStrings;
  STRING          line, name;
  SpenvisCSV inCSV;
  int i;
  char   blockName[5];

  std::cout << "SpenvisCSVCollection::ReadCollection\n";
  while ( getline( ins, line) ) {
    inStrings.push_back( line);
  }

 /* for (int i=0; i< inStrings.size(); i++) 
    std::cout << i << "  '" << inStrings[i] << "'\n";*/

  i = 1;
  do {
    inCSV = GetBlock( inStrings, i);
    if (inCSV.GetNumHeaderLines() > 1) {
      sprintf(blockName,"%5.5d", i);
      name = blockName;
      AddCSVBlock( name, inCSV);
    }
    i++;
  }
  while( inCSV.GetNumHeaderLines() > 1);

}

void SpenvisCSVCollection::ReadCollection( ){

}

void SpenvisCSVCollection::ReadCollection( STRING filename){
  ifstream ins;
  
  ins.open( filename.c_str());
  if (!ins) 
    std::cerr << "SpenvisCSVCollection::ReadCollection('" 
	      << filename 
	      << "', unable to open file for input.\n";
  else {
    ReadCollection( ins);
    ins.close();
  }
}


SpenvisCSV SpenvisCSVCollection::GetBlock( std::vector <STRING> inStrings, int blockNo) {

  int blockStart, blockEnd, currBlockNo;
  STRING startOfBlock = "'*'";
  std::vector <STRING> block;
  SpenvisCSV      outBlock;
  SpenvisCSV::HeaderLine      hl;

  // Search for the start of blocks by the '*' sequence in the header line.

  blockStart = 0;
  currBlockNo = 0;
  outBlock.ResetBlock();

  // Find the start of the required block in the input strings
  while( (blockStart < int(inStrings.size()) ) && (currBlockNo < blockNo) ){
    if ( inStrings[ blockStart].find( startOfBlock) == 0)
      currBlockNo++;
    if (currBlockNo < blockNo) blockStart++;
  }

  if (currBlockNo < blockNo) return outBlock;

  outBlock.ParseBlockHeaderLine( inStrings[ blockStart],  &hl );

  blockEnd = blockStart + hl.nHeader;

  // Handle the case where the number of body rows not know.
  // In this case, we assume that the block ends either at the end of the vector, or
  // with the first string that contains an apostrophy.
  if (hl.nBody < 0) {
    while ( ( blockEnd < int(inStrings.size())) && 
	    (inStrings[blockEnd].find('\'') ==  inStrings[blockEnd].npos))
      blockEnd++;
  } else 
    blockEnd += hl.nBody;

  for (int i=blockStart; i<blockEnd; i++) 
    block.push_back( inStrings[i] );

  outBlock.ParseBlock( block);
  
  return outBlock;
}
