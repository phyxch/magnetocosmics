// Updated on 8/6/2014: Hexc, Olesya
// Clean up old header file inclusions (i.e., removing ".h"
//
/*
 *   Author:   $Author: hexc $
 *   Date:     $Date: 2014/08/06 15:59:44 $
 *   Revision: $Revision: 1.2 $
 *   Tag Name: $Name:  $
 *   
 *   $Log: SpenvisCSV.hh,v $
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


SpenvisCSV class - 

provides an encapsulated class for the creation and writing of SPENVIS Comma
Separated Value blocks (also known as the MK format). Details of the file format
can be found on the SPENVIS website help pages (http://www.spenvis.oma.be/spenvis/)

The methods of this class are:

AddComment            - adds a comment line to the CSV file
AddMetaVariable       - adds a metavariable to the CSV file. Metavariables
                        are specified with a name, dimension of the meta-variable,
                        and the data. 
AddMetaVariableStr    - As AddMetaVariable, but adds a string type metavariable - the
                        dimension is not required in this case.
AddAnnotation         - The annotation feature is currently not in use, but is included
                        for future use by the SPENVIS project. To ensure compatability
			please refrain from its use.
AddVariable           - Specify a variable that is to be stored in the file. The variables
                        are specified with: a name; unit, e.g. 'cm-2 s-1'; dimension; a
			textual description; and the C-style format descriptor, which is 
			used when the data is output.
AddDataRow            - This method adds a record to the CSV data structure. At this time,
                        only a double array is permitted. The size of this array must match
			the number of variables and associated dimensions. So if, for example,
			two variables, epoch and data, with respective record dimensions of 1
			and 5 are specified, then the data array provided to this method MUST be
			at least 6 elements. NO RANGE CHECKING is preformed.
EraseData             - This method deletes all data previously entered to the class using the
                        AddDataRow method. The comments, metavariables, annotations and variable
			definitions are not affected by this call.

GetNumHeaderLines     - Returns the total number of lines in the block.
GetNumTextLines       - Returns the number of comment/text lines in the block.
GetNumMetaVariables   - Returns the number of metaVariables in the block.
GetNumAnnotationLines - Returns the number of annotation lines in the block.
GetNumVariables       - Returns the number of variables in the block.
GetNumDataColumns     - Returns the number of data columns in the data section of the block.
GetNumDataLines       - Returns the number of data lines in the data section.
GetNumFollowingBlocks - Returns the number of blocks in the data file after the current one.

GetComment            - Returns the i-th comment.
GetMetaVariableLine   - Returns a string with the MetaVariable information in the format specified
                        by the SPENVIS CSV format.
GetMetaVariableName   - Returns the name of the i-th metavariable.
GetVariableLine       - Returns a string with the Variable information in the format specified
                        by the SPENVIS CSV format.
GetVariableName       - Returns the name of the i-th variable.
GetVariableValue      - Returns the value of the variable in the given row of the data section.
GetDataLine           - Returns a string with the i-th Variable data information in the format specified
                        by the SPENVIS CSV format.
GetDataRecord         - Returns a vector of <double>s containing the data from the i-th row of the
                        data section.
GetBlockHeaderLine    - Returns a string of the start of block line in the format specified
                        by the SPENVIS CSV format. The nblocks argument is to specify the number
			of blocks to follow in the file; for a file containing a single block, this
			should be 0, for a file with two blocks, it should be 1, etc..
OutputCSV             - Writes the output to either the ofstream specified or std::cout, depending
                        on which method is called.
ResetBlock            - Erases the contents of the block and resets the various internal data structures.
ParseBlock            - Given a vector <STRING> input, converts the contents, which should be in the
                        SpenvisCSV MK format, and sets the internal data structures to the corresponding
			contents.
ParseBlockHeaderLine  - given a string containing a header line, the contents is parsed and used to set
                        the headerLine variable is initialised. A corresponding method also returns
			a HeaderLine structure with the parsed contents.
ParseCommentLine      - Given a string in the comment line format of the SpenvisCSV files, this method
                        parses it and stores the comment in the block.
ParseMetaDataLine     - Given a string in the MetaData line format of the SpenvisCSV files, this method
                        parses it and stores it in the block.
ParseAnnotationLine   - Given a string in the annotation line format of the SpenvisCSV files, this method
                        parses it and stores it in the block.
ParseVariableLine     - Given a string in the variable line format of the SpenvisCSV files, this method
                        parses it and stores it in the block.
ParseDataRow          - Given a string in the data row line format of the SpenvisCSV files, this method
                        parses it and stores it in the data section of the block.

PLEASE Note, that very little, if any effort was expended to validate inputs to the various methods.
It is highly likely that segmentation faults, etc. could occur due to improperly formulated calls
to the methods.


H. Evans, 
September 7, 2004


Change History:
---------------
v1.0 - released 2004/9/7.


*/


#ifndef SpenvisCSV_h
#define SpenvisCSV_h 1

//#include "fstream.h"
//#include <vector.h>
//#include <map.h>
//#include <string.h>

#include "fstream"
#include <vector>
#include <map>
#include <string>

#include"G4String.hh"

using namespace std;


#define STRING std::basic_string<char>
//#define STRING G4String

class SpenvisCSV
{
public:
  SpenvisCSV        ();
  ~SpenvisCSV       ();

  void AddComment        ( STRING comment);
  void AddMetaVariable   ( STRING name, int dims, float vals[], STRING format, STRING unit="");
  void AddMetaVariable   ( STRING name, std::vector <double> vals, STRING format, STRING unit="");
  void AddMetaVariableStr( STRING name, STRING val);
  void AddAnnotation     ( STRING annotation );
  void AddVariable       ( STRING name, STRING unit, int dims, STRING description, STRING fmt );
  void AddDataRow        ( double *rec);
  void AddDataRow        ( std::vector <double> rec);
  void EraseData         (){
    vals.erase( vals.begin(), vals.end());
  };

  int GetNumHeaderLines() {
    return 1 + GetNumTextLines() + 
               GetNumMetaVariables() + 
               GetNumAnnotationLines() + 
               GetNumVariables();
  };
  int GetNumTextLines() {
    return comments.size();
  };
  int GetNumMetaVariables() {
    return metaVar.size();
  };
  int GetNumAnnotationLines() {
    return annotations.size();
  };
  int GetNumVariables(){
    return vars.size();
  };
  int GetNumDataColumns() {
    return nvar_cols;
  };
  int GetNumDataLines() {
    return vals.size();
  };
  int GetNumFollowingBlocks() {
    return headerLine.npart;
  };

  STRING GetComment( int i);
  STRING GetMetaVariableLine( int i);
  STRING GetMetaVariableName( int i);
  STRING GetVariableLine( int i);
  STRING GetVariableName( int i);
  STRING GetVariableValue( STRING name, int row);
  STRING GetDataLine(int i);
  std::vector<double> GetDataRecord( int recordNum);
  std::vector<double> GetDataRecord( int recordNum, int varNum);
  std::vector<double> GetDataRecord( STRING variable, int recordNum);
  
  STRING GetBlockHeaderLine( int nblocks);
  
  void OutputCSV         ( ofstream &outs, int nFollowingBlocks);
  void OutputCSV         (int nFollowingBlocks);

  struct HeaderLine {
    int nHeader;        // total number of records in the header section (including header line)
    int nText;          // number of text lines included in the first part of the header section
    int nMeta;          // number of metavariables
    int nFut;           // num of records used for application specific annotation
    int nVar;           // number of variables in the body section
    int nCol;           // number of columns in the body section
    int nBody;          // number of records in the body section (-1 indicates this is unknown)
    int npart;          // number of remaining blocks associated to the current namelist section
  };
  void ResetBlock();
  void ParseBlock( std::vector <STRING> ins );
  void ParseBlockHeaderLine( STRING ins);
  void ParseBlockHeaderLine( STRING ins, HeaderLine *hl);
  void ParseCommentLine ( STRING ins);
  void ParseMetaDataLine( STRING ins);
  void ParseAnnotationLine( STRING ins);
  void ParseVariableLine( STRING ins);
  void ParseDataRow( STRING ins);


private:

  std::vector <STRING> strsplit( STRING ins, char delim);
  STRING strremove( STRING ins, char target);   // removes all target characters from the string
  STRING strtrim( STRING ins, char target);     // Trims leading/trailing target characters

  int GetVariableNumber( STRING name);

  struct Comments {
    STRING val;
  };

  struct MetaVar {
    STRING name;
    int    dims;
    std::vector<double> fvals;

    STRING svals;   
    STRING  unit;
    STRING format;
  };

  struct Vars {
    STRING  name;
    STRING  unit;
    int     dims;
    STRING  descr;
    STRING  format;
  };

  // A single record in the data section
  struct Vals {
    std::vector<double> values;
  };


  int nvar_cols;
  HeaderLine headerLine;

  std::vector<Comments>      comments;
  std::vector<STRING>        annotations;
//	changes by ldesorgher  
//  std::map<string, MetaVar>  metaVar;
  std::vector<MetaVar>	     metaVar;
  std::vector<Vars>          vars;
  std::vector<STRING>        col_formats;  // list of formats for each column.
  std::vector<Vals>          vals;


};

#endif
