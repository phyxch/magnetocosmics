// Updated on 8/6/2014: Hexc, Olesya
// Clean up old header file inclusions (i.e., removing ".h"
/*
 *   Author:   $Author: hexc $
 *   Date:     $Date: 2014/08/06 15:59:44 $
 *   Revision: $Revision: 1.2 $
 *   Tag Name: $Name:  $
 *   
 *   $Log: SpenvisCSVCollection.hh,v $
 *   Revision 1.2  2014/08/06 15:59:44  hexc
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


/*


SpenvisCSVCollection class - 

provides an encapsulated container class for the creation and writing of SPENVIS Comma
Separated Value blocks (also known as the MK format). Details of the file format
can be found on the SPENVIS website help pages (http://www.spenvis.oma.be/spenvis/)

The methods of this class are:
------------------------------
AddCSVBlock                   - Adds a preconstructed SpenvisCSV block to the container
DelCSVBlock                   - Deletes a previously added SpenvisCSV block from the container.
GetCSVBlock                   - returns a SpenvisCSV block, previously added to the container.
OutputCollection              - Outputs the entire collection in SPENVIS CSV format to either a
                                previously opened file (ofstream) or std::cout, depending on the
				method called.
ReadCollection                - Reads an entire collection of SPENVIS CSV blocks from an opened 
                                file, std::cin or a file name, depending on the method called. 
				The blocks are given a name representative of the order in which
				they were read, e.g. '00001', '00002'.

The blocks are added with a name and stored in a <map> construct, so they will be outputted in
alphabetical order of their names!!

PLEASE Note, that very little, if any effort was expended to validate inputs to the various methods.
It is highly likely that segmentation faults, etc. could occur due to improperly formulated calls
to the methods.

H. Evans, 
September 7, 2004


Change History:
---------------
v1.0 - released 2004/9/7.


*/

#ifndef SpenvisCSVCollection_h
#define SpenvisCSVCollection_h 1

//#include <map.h>
#include <map>
//#include <string.h>
#include <iostream>
#include <string>
#include <fstream>
//#include <fstream.h>
#include <cstdlib>
#include "SpenvisCSV.hh"

using namespace std;

class SpenvisCSVCollection
{
public:
  SpenvisCSVCollection (){};
  ~SpenvisCSVCollection(){};

  void AddCSVBlock( STRING name, SpenvisCSV &newblock){
    if (collection.find(name) != collection.end()) {
      collection[name] = newblock;
    } else {
      DelCSVBlock( name);
      collection[name] = newblock;
    }
  };
  void DelCSVBlock( STRING name) {
    int n_erased;
    n_erased = collection.erase(name);
  };

  SpenvisCSV GetCSVBlock( STRING name) {
  SpenvisCSV sp;

    if ( collection.find(name) != collection.end())
      sp = collection.find(name)->second;
    
    return sp;
  }

  void OutputCollection( ofstream &outs);           // Output to an ofstream (previously opened file
  void OutputCollection( );                         // Output to std::cout
  void OutputCollection( STRING filename) {         // Output to a file and close the file.
    ofstream outstream;
    outstream.open( filename.c_str());
    OutputCollection( outstream);
    outstream.close();
  }

  void ReadCollection( ifstream &ins);            // Read from an ifstream, previously opened file
  void ReadCollection( );                         // Read from std:cin
  void ReadCollection( STRING filename);          // Read from a file.

private:
  std::map<STRING, SpenvisCSV> collection;

  // Given a previously read in vector of strings, which should contain a series of 
  // SpenvisCSV blocks in their textual representation, this method parses the vector
  // and returns a SpenvisCSV block with the contents.
  SpenvisCSV GetBlock( std::vector <STRING> inStrings, int blockNo);

};

#endif
