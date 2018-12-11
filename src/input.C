
#include "input.h"
#include "cluster.h"

// Shouldn't need to define this here again!?
void Rewind(ifstream& infile) {
  infile.clear();
  infile.seekg(0,ios::beg);
}

string ReadG09Section(ifstream& infile) {
  string rem = "\n";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$G09") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }
  //rem += "$end\n\n";

  infile.clear();
  return rem;
}

string ReadPSI4Section(ifstream& infile) {
  string rem = "\n";  
  string line;        
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$PSI4") {
      for (int i=0;;i++) {
        getline(infile,line);
        if (line.substr(0,4) != "$end") {
          rem += line;
          rem += "\n";
 	} else {
	  break;
        }
      }
    }
  }

  infile.clear();
  return rem;    


}

string ReadDaltonSection(ifstream& infile) {
  string rem = "\n";
  string line;
  Rewind(infile);
  
  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$DALTON") {
      for (int i=0;;i++) {
	getline(infile,line);
	if (line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();
  return rem;


}


// Custom Basis Sections:
string ReadHBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$H" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "")  rem.erase(rem.length()-1);

  return rem;
}


string ReadCBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$C" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string ReadNBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$N" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string ReadOBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$O" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string ReadSBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$S" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string ReadClBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Cl" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string ReadIBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$I" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string ReadSnBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Sn" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string ReadPBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$P" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string ReadKBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$K" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


string ReadNaBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Na" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string ReadBrBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$Br" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}

string ReadFBasis( ifstream& infile) {
  string rem = "";
  string line;
  Rewind(infile);

  while ( !infile.eof() ) {
    getline(infile,line);
    if ( line.substr(0,9)=="$F" ) {
      for( int i=0;;i++) {
	getline(infile,line);
	if ( line.substr(0,4) != "$end") {
	  rem += line;
	  rem += "\n";
	} else {
	  break;
	}
      }
    }
  }

  infile.clear();

  // remove blank line at the end of the string
  if ( rem != "") rem.erase(rem.length()-1);

  return rem;
}


