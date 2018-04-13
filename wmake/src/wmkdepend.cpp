
#line 1 "wmkdepend.rl"
/*---------------------------------*- C -*-----------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2018 OpenCFD Ltd.
    \\/      M anipulation   |
------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    wmkdepend

Description
    A fast dependency list generator that emulates the behaviour and
    output of cpp -M.

    The lexer for include statements uses a Ragel-generated lexer and
    then searches the files found. Each file is visited only once,
    which helps make this faster than cpp. It also handles missing
    files somewhat more gracefully.

    The goto-based finite state machine (FSM) is generated with

        ragel -G2 -o wmkdepend.cpp wmkdepend.rl

Usage
    wmkdepend [-Idir..] [-iheader...] [-eENV...] [-oFile] [-q] filename

Note
    May not capture all possible corner cases or line continuations such as

        #include \
            "file.H"

\*---------------------------------------------------------------------------*/
/*
 * With cpp:
 *
 * cpp -x c++ -std=c++11
 *     -M -DWM_$(WM_PRECISION_OPTION) -DWM_LABEL_SIZE=$(WM_LABEL_SIZE) |
 * sed -e 's,^$(WM_PROJECT_DIR)/,$$(WM_PROJECT_DIR)/,' \
 *     -e 's,^$(WM_THIRD_PARTY_DIR)/,$$(WM_THIRD_PARTY_DIR)/,'
*/

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <forward_list>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Length of the input read buffer
#define READ_BUFLEN 16384

// The executable name (for messages), without requiring access to argv[]
#define EXENAME  "wmkdepend"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//- Program usage
void usage()
{
    std::cerr
        <<
        "\nUsage: " EXENAME
        " [-Idir...] [-iheader...] [-eENV...] [-oFile] [-q]"
        " filename\n\n"
        "  -Idir     Directories to be searched for headers.\n"
        "  -iheader  Headers to be ignored.\n"
        "  -eENV     Environment variable path substitutions.\n"
        "  -oFile    Write output to File.\n"
        "  -q        Suppress 'No such file' warnings.\n"
        "  -v        Report each include file to stderr.\n"
        "\nDependency list generator, similar to 'cpp -M'\n\n";
}

//- Suppress some error messages
bool optQuiet = false;

//- Verbose progress
bool optVerbose = false;

//- The top-level source file being processed
std::string sourceFile;


//- All file opening and writing
namespace Files
{
    //- Include directories
    static std::vector<std::string> dirs;

    //- Set of files already visited
    static std::unordered_set<std::string> visited;

    //- Environment substitution
    struct envEntry
    {
        std::string name;
        std::string value;
        size_t len;

        envEntry(std::string&& k, std::string&& v)
        :
            name(std::move(k)),
            value(std::move(v)),
            len(value.size())
        {}
    };

    //- List of environ variables to substitute
    std::forward_list<envEntry> envlist;

    //- Clear, reset all variables
    void reset(unsigned ndirs = 0u)
    {
        dirs.clear();
        visited.clear();
        envlist.clear();
        if (ndirs)
        {
            dirs.reserve(ndirs);
        }
    }


    //- Add environ replacement
    //
    // Eg for 'WM_PROJECT_DIR' stores the pair
    //     '$(WM_PROJECT_DIR)/'
    //     '/path/openfoam/project/'
    //
    // No-op if the environment does not exists or is empty
    void addEnv(const char* key)
    {
        const size_t keylen = key ? strlen(key) : 0;
        if (!keylen) return;

        const char *val = ::getenv(key);

        if (!val || !*val) return;

        // The "$(ENV)/" portion
        std::string oldText;
        oldText.reserve(keylen+4);
        oldText.append("$(").append(key,keylen).append(")/");

        // The "/env/value/" portion
        std::string newText(val);
        if (newText.back() != '/')
        {
            newText.append(1, '/');
        }

        envlist.emplace_front(std::move(oldText), std::move(newText));
    }


    //- Open a file for reading and emit its qualified name to stdout.
    //
    //  Uses env substitutions at the beginning of the path
    //
    //  Eg,
    //      /path/openfoam/project/directory/name
    //   -> $(WM_PROJECT_DIR)/directory/name
    //
    //  \return nullptr on failure
    FILE* openAndEmit(const std::string& fileName)
    {
        const auto len = fileName.size();
        const char *fname = fileName.c_str();

        FILE *infile = ::fopen(fname, "r");

        if (infile)
        {
            // Mark as having been visited
            visited.insert(fileName);

            for (const auto& entry : envlist)
            {
                if
                (
                    len > entry.len
                 && !fileName.compare(0, entry.len, entry.value)
                )
                {
                    fname += entry.len;  // Now positioned after the '/'
                    fputs(entry.name.c_str(), stdout);
                    break;
                }
            }

            fputs(fname, stdout);
            fputs(" \\\n", stdout);
        }
        else if (errno == EMFILE)
        {
            std::cerr
                << EXENAME ": too many open files while opening '"
                << fileName << "'\n"
                << "Please change your open descriptor limit\n";
        }

        return infile;
    }


    //- Open a not previously visited file for reading,
    //  using the include dirs as required.
    //
    //  On success, emits the resolved name on stdout
    //
    //  \return nullptr on failure
    FILE* open(const std::string& fileName)
    {
        // Bad file name, or already visited
        if (fileName.empty() || visited.find(fileName) != visited.end())
        {
            return nullptr;
        }

        FILE* infile = openAndEmit(fileName);
        if (!infile)
        {
            std::string fullName;

            for (const auto& d : dirs)
            {
                if (d.empty())
                {
                    continue;
                }

                fullName.clear();
                fullName.reserve(d.size() + fileName.size() + 1);

                fullName.append(d);
                if (fullName.back() != '/')
                {
                    fullName.append("/");
                }
                fullName.append(fileName);

                infile = openAndEmit(fullName);

                if (infile)
                {
                    break;
                }
            }
        }

        // Mark as having been visited.
        // This also avoids re-testing and multiple messages.
        visited.insert(fileName);

        // Report failues
        if (!infile && !optQuiet)
        {
            std::cerr
                << EXENAME ": could not open '"
                << fileName << "' for source file '"
                << sourceFile << "'";

            if (dirs.size())
            {
                std::cerr << ": " << strerror(errno);
            }

            std::cerr << '\n' << std::flush;
        }

        return infile;
    }

} // end of namespace Files


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names


#line 324 "wmkdepend.rl"



//
// FSM globals
//


#line 316 "wmkdepend.cpp"
static const int wmkdep_start = 167;
static const int wmkdep_error = 0;

static const int wmkdep_en_comment = 165;
static const int wmkdep_en_main = 167;


#line 332 "wmkdepend.rl"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//
// Open a file and process
//
void processFile(const std::string& fileName)
{
    FILE* infile = Files::open(fileName);
    if (optVerbose)
    {
        std::cerr << fileName << '\n';
    }
    if (!infile) return;

    // ts, te  = Ragel token start/end points (default naming)
    // act, cs = Ragel action, code state, respectively (default naming)
    char *ts, *te;
    int act, cs;

    
#line 347 "wmkdepend.cpp"
	{
	cs = wmkdep_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 353 "wmkdepend.rl"
   /* ^^^ FSM initialization here ^^^ */;

    // Local token start
    char *tok = nullptr;

    // Buffering
    char inbuf[READ_BUFLEN];
    size_t pending = 0;

    // Processing loop (as per Ragel pdf example)
    for (bool good = true; good; /*nil*/)
    {
        char *data = inbuf + pending;                // current data buffer
        const size_t buflen = READ_BUFLEN - pending; // space left in buffer

        if (!buflen)
        {
            // We overfilled the buffer while trying to scan a token...
            std::cerr
                << EXENAME ": buffer full while scanning '"
                << fileName << "'\n";
            break;
        }

        // p,pe = Ragel parsing point and parsing end (default naming)
        // eof  = Ragel EOF point (default naming)

        const size_t gcount = ::fread(data, 1, buflen, infile);

        char *p = data;
        char *pe = p + gcount;
        char *eof = nullptr;

        if (::feof(infile))
        {
            // Tag 'pe' as being the EOF for the FSM as well
            eof = pe;
            good = false;
        }
        else if (!gcount)
        {
            break;
        }

        
#line 401 "wmkdepend.cpp"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr0:
#line 1 "NONE"
	{	switch( act ) {
	case 0:
	{{goto st0;}}
	break;
	case 4:
	{{p = ((te))-1;} {goto st165;} }
	break;
	default:
	{{p = ((te))-1;}}
	break;
	}
	}
	goto st167;
tr3:
#line 321 "wmkdepend.rl"
	{te = p+1;}
	goto st167;
tr20:
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{te = p+1;}
	goto st167;
tr82:
#line 321 "wmkdepend.rl"
	{{p = ((te))-1;}}
	goto st167;
tr84:
#line 317 "wmkdepend.rl"
	{te = p+1;}
	goto st167;
tr112:
#line 316 "wmkdepend.rl"
	{te = p+1;}
	goto st167;
tr172:
#line 320 "wmkdepend.rl"
	{te = p+1;}
	goto st167;
tr227:
#line 313 "wmkdepend.rl"
	{te = p;p--;}
	goto st167;
tr228:
#line 321 "wmkdepend.rl"
	{te = p;p--;}
	goto st167;
tr229:
#line 316 "wmkdepend.rl"
	{te = p;p--;}
	goto st167;
st167:
#line 1 "NONE"
	{ts = 0;}
#line 1 "NONE"
	{act = 0;}
	if ( ++p == pe )
		goto _test_eof167;
case 167:
#line 1 "NONE"
	{ts = p;}
#line 470 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 34: goto st57;
		case 35: goto st3;
		case 39: goto st79;
		case 47: goto st125;
		case 76: goto st164;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( (*p) == 10 )
		goto tr3;
	goto st2;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	switch( (*p) ) {
		case 10: goto tr6;
		case 32: goto st4;
		case 35: goto st3;
		case 105: goto st24;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st4;
	goto st1;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 10: goto tr6;
		case 32: goto st4;
		case 105: goto st15;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st4;
	goto st2;
tr6:
#line 1 "NONE"
	{te = p+1;}
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st168;
tr175:
#line 1 "NONE"
	{te = p+1;}
#line 320 "wmkdepend.rl"
	{act = 5;}
	goto st168;
st168:
	if ( ++p == pe )
		goto _test_eof168;
case 168:
#line 543 "wmkdepend.cpp"
	switch( (*p) ) {
		case 32: goto st5;
		case 105: goto st6;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st5;
	goto tr0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 32: goto st5;
		case 105: goto st6;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st5;
	goto tr0;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 110 )
		goto st7;
	goto tr0;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 99 )
		goto st8;
	goto tr0;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 108 )
		goto st9;
	goto tr0;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 117 )
		goto st10;
	goto tr0;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 100 )
		goto st11;
	goto tr0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( (*p) == 101 )
		goto st12;
	goto tr0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 32: goto st12;
		case 34: goto st13;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st12;
	goto tr0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 34 )
		goto tr0;
	goto tr18;
tr18:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 630 "wmkdepend.cpp"
	if ( (*p) == 34 )
		goto tr20;
	goto st14;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	switch( (*p) ) {
		case 10: goto tr3;
		case 110: goto st16;
	}
	goto st2;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	switch( (*p) ) {
		case 10: goto tr3;
		case 99: goto st17;
	}
	goto st2;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 10: goto tr3;
		case 108: goto st18;
	}
	goto st2;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	switch( (*p) ) {
		case 10: goto tr3;
		case 117: goto st19;
	}
	goto st2;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	switch( (*p) ) {
		case 10: goto tr3;
		case 100: goto st20;
	}
	goto st2;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	switch( (*p) ) {
		case 10: goto tr3;
		case 101: goto st21;
	}
	goto st2;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	switch( (*p) ) {
		case 10: goto tr27;
		case 32: goto st21;
		case 34: goto st22;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st21;
	goto st2;
tr27:
#line 1 "NONE"
	{te = p+1;}
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st169;
tr184:
#line 1 "NONE"
	{te = p+1;}
#line 320 "wmkdepend.rl"
	{act = 5;}
	goto st169;
st169:
	if ( ++p == pe )
		goto _test_eof169;
case 169:
#line 716 "wmkdepend.cpp"
	switch( (*p) ) {
		case 32: goto st12;
		case 34: goto st13;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st12;
	goto tr0;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	switch( (*p) ) {
		case 10: goto tr30;
		case 34: goto st2;
	}
	goto tr29;
tr29:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 741 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
	}
	goto st23;
tr32:
#line 1 "NONE"
	{te = p+1;}
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st170;
tr30:
#line 1 "NONE"
	{te = p+1;}
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st170;
tr134:
#line 1 "NONE"
	{te = p+1;}
#line 316 "wmkdepend.rl"
	{act = 2;}
	goto st170;
tr189:
#line 1 "NONE"
	{te = p+1;}
#line 320 "wmkdepend.rl"
	{act = 5;}
	goto st170;
tr187:
#line 1 "NONE"
	{te = p+1;}
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
#line 320 "wmkdepend.rl"
	{act = 5;}
	goto st170;
st170:
	if ( ++p == pe )
		goto _test_eof170;
case 170:
#line 785 "wmkdepend.cpp"
	if ( (*p) == 34 )
		goto tr20;
	goto st14;
tr33:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st171;
tr79:
#line 1 "NONE"
	{te = p+1;}
#line 317 "wmkdepend.rl"
	{act = 3;}
	goto st171;
tr108:
#line 1 "NONE"
	{te = p+1;}
#line 316 "wmkdepend.rl"
	{act = 2;}
	goto st171;
st171:
	if ( ++p == pe )
		goto _test_eof171;
case 171:
#line 813 "wmkdepend.cpp"
	if ( (*p) == 10 )
		goto tr3;
	goto st2;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
		case 110: goto st25;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
		case 99: goto st26;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
		case 108: goto st27;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
		case 117: goto st28;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
		case 100: goto st29;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
		case 101: goto st30;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 10: goto tr27;
		case 32: goto st21;
		case 34: goto st31;
		case 35: goto st3;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st21;
	goto st1;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 10: goto tr30;
		case 32: goto tr29;
		case 34: goto st1;
		case 35: goto tr42;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr29;
	goto tr41;
tr41:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st32;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
#line 929 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
tr44:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st172;
tr76:
#line 1 "NONE"
	{te = p+1;}
#line 317 "wmkdepend.rl"
	{act = 3;}
	goto st172;
tr106:
#line 1 "NONE"
	{te = p+1;}
#line 316 "wmkdepend.rl"
	{act = 2;}
	goto st172;
tr169:
#line 1 "NONE"
	{te = p+1;}
#line 319 "wmkdepend.rl"
	{act = 4;}
	goto st172;
st172:
	if ( ++p == pe )
		goto _test_eof172;
case 172:
#line 969 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
tr42:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st33;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
#line 986 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr47;
		case 32: goto st34;
		case 34: goto tr44;
		case 35: goto st33;
		case 105: goto st50;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st34;
	goto st32;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 10: goto tr47;
		case 32: goto st34;
		case 34: goto tr33;
		case 105: goto st43;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st34;
	goto st23;
tr47:
#line 1 "NONE"
	{te = p+1;}
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st173;
tr204:
#line 1 "NONE"
	{te = p+1;}
#line 320 "wmkdepend.rl"
	{act = 5;}
	goto st173;
st173:
	if ( ++p == pe )
		goto _test_eof173;
case 173:
#line 1026 "wmkdepend.cpp"
	switch( (*p) ) {
		case 32: goto st35;
		case 34: goto tr20;
		case 105: goto st36;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st35;
	goto st14;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 32: goto st35;
		case 34: goto tr20;
		case 105: goto st36;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st35;
	goto st14;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 34: goto tr20;
		case 110: goto st37;
	}
	goto st14;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 34: goto tr20;
		case 99: goto st38;
	}
	goto st14;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 34: goto tr20;
		case 108: goto st39;
	}
	goto st14;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 34: goto tr20;
		case 117: goto st40;
	}
	goto st14;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 34: goto tr20;
		case 100: goto st41;
	}
	goto st14;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 34: goto tr20;
		case 101: goto st42;
	}
	goto st14;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 32: goto st42;
		case 34: goto tr58;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st42;
	goto st14;
tr58:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st174;
st174:
	if ( ++p == pe )
		goto _test_eof174;
case 174:
#line 1124 "wmkdepend.cpp"
	if ( (*p) == 34 )
		goto tr227;
	goto tr18;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
		case 110: goto st44;
	}
	goto st23;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
		case 99: goto st45;
	}
	goto st23;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
		case 108: goto st46;
	}
	goto st23;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
		case 117: goto st47;
	}
	goto st23;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
		case 100: goto st48;
	}
	goto st23;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
		case 101: goto st49;
	}
	goto st23;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 10: goto tr65;
		case 32: goto st49;
		case 34: goto tr66;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st49;
	goto st23;
tr65:
#line 1 "NONE"
	{te = p+1;}
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st175;
tr213:
#line 1 "NONE"
	{te = p+1;}
#line 320 "wmkdepend.rl"
	{act = 5;}
	goto st175;
st175:
	if ( ++p == pe )
		goto _test_eof175;
case 175:
#line 1216 "wmkdepend.cpp"
	switch( (*p) ) {
		case 32: goto st42;
		case 34: goto tr58;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st42;
	goto st14;
tr66:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st176;
tr95:
#line 1 "NONE"
	{te = p+1;}
#line 317 "wmkdepend.rl"
	{act = 3;}
	goto st176;
st176:
	if ( ++p == pe )
		goto _test_eof176;
case 176:
#line 1242 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr30;
		case 34: goto st2;
	}
	goto tr29;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
		case 110: goto st51;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
		case 99: goto st52;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
		case 108: goto st53;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
		case 117: goto st54;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
		case 100: goto st55;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
		case 101: goto st56;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 10: goto tr65;
		case 32: goto st49;
		case 34: goto tr73;
		case 35: goto st33;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st49;
	goto st32;
tr73:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st177;
tr102:
#line 1 "NONE"
	{te = p+1;}
#line 317 "wmkdepend.rl"
	{act = 3;}
	goto st177;
st177:
	if ( ++p == pe )
		goto _test_eof177;
case 177:
#line 1363 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr30;
		case 32: goto tr29;
		case 34: goto st1;
		case 35: goto tr42;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr29;
	goto tr41;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st58;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 10: goto tr3;
		case 34: goto tr79;
		case 92: goto st59;
	}
	goto st58;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	if ( (*p) == 10 )
		goto tr81;
	goto st58;
tr81:
#line 1 "NONE"
	{te = p+1;}
	goto st178;
st178:
	if ( ++p == pe )
		goto _test_eof178;
case 178:
#line 1412 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr228;
		case 34: goto tr84;
		case 92: goto st61;
	}
	goto st60;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 10: goto tr82;
		case 34: goto tr84;
		case 92: goto st61;
	}
	goto st60;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	goto st60;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	switch( (*p) ) {
		case 10: goto tr6;
		case 32: goto st63;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
		case 105: goto st72;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st63;
	goto st57;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 10: goto tr6;
		case 32: goto st63;
		case 34: goto tr79;
		case 92: goto st59;
		case 105: goto st64;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st63;
	goto st58;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	switch( (*p) ) {
		case 10: goto tr3;
		case 34: goto tr79;
		case 92: goto st59;
		case 110: goto st65;
	}
	goto st58;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 10: goto tr3;
		case 34: goto tr79;
		case 92: goto st59;
		case 99: goto st66;
	}
	goto st58;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 10: goto tr3;
		case 34: goto tr79;
		case 92: goto st59;
		case 108: goto st67;
	}
	goto st58;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	switch( (*p) ) {
		case 10: goto tr3;
		case 34: goto tr79;
		case 92: goto st59;
		case 117: goto st68;
	}
	goto st58;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 10: goto tr3;
		case 34: goto tr79;
		case 92: goto st59;
		case 100: goto st69;
	}
	goto st58;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 10: goto tr3;
		case 34: goto tr79;
		case 92: goto st59;
		case 101: goto st70;
	}
	goto st58;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 10: goto tr27;
		case 32: goto st70;
		case 34: goto tr95;
		case 92: goto st59;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st70;
	goto st58;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 10: goto tr81;
		case 32: goto st58;
		case 35: goto st62;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st58;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
		case 110: goto st73;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st58;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
		case 99: goto st74;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st58;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
		case 108: goto st75;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st58;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
		case 117: goto st76;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st58;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
		case 100: goto st77;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st58;
		case 34: goto tr76;
		case 35: goto st62;
		case 92: goto st71;
		case 101: goto st78;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st58;
	goto st57;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 10: goto tr27;
		case 32: goto st70;
		case 34: goto tr102;
		case 35: goto st62;
		case 92: goto st71;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st70;
	goto st57;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
	}
	goto st80;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	if ( (*p) == 10 )
		goto tr110;
	goto st80;
tr110:
#line 1 "NONE"
	{te = p+1;}
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st179;
tr133:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st179;
st179:
	if ( ++p == pe )
		goto _test_eof179;
case 179:
#line 1707 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr0;
		case 39: goto tr112;
		case 92: goto st83;
	}
	goto st82;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	switch( (*p) ) {
		case 10: goto tr0;
		case 39: goto tr112;
		case 92: goto st83;
	}
	goto st82;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	goto st82;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 10: goto tr6;
		case 32: goto st85;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
		case 105: goto st99;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st85;
	goto st79;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 10: goto tr6;
		case 32: goto st85;
		case 39: goto tr108;
		case 92: goto st81;
		case 105: goto st86;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st85;
	goto st80;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
		case 110: goto st87;
	}
	goto st80;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
		case 99: goto st88;
	}
	goto st80;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
		case 108: goto st89;
	}
	goto st80;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
		case 117: goto st90;
	}
	goto st80;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
		case 100: goto st91;
	}
	goto st80;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
		case 101: goto st92;
	}
	goto st80;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	switch( (*p) ) {
		case 10: goto tr27;
		case 32: goto st92;
		case 34: goto st93;
		case 39: goto tr108;
		case 92: goto st81;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st92;
	goto st80;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	switch( (*p) ) {
		case 10: goto tr30;
		case 34: goto st80;
		case 39: goto tr125;
		case 92: goto tr126;
	}
	goto tr124;
tr124:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st94;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
#line 1857 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
	}
	goto st94;
tr128:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st180;
st180:
	if ( ++p == pe )
		goto _test_eof180;
case 180:
#line 1877 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr3;
		case 39: goto tr108;
		case 92: goto st81;
	}
	goto st80;
tr129:
#line 1 "NONE"
	{te = p+1;}
#line 316 "wmkdepend.rl"
	{act = 2;}
	goto st181;
tr125:
#line 1 "NONE"
	{te = p+1;}
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
#line 316 "wmkdepend.rl"
	{act = 2;}
	goto st181;
st181:
	if ( ++p == pe )
		goto _test_eof181;
case 181:
#line 1902 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr33;
	}
	goto st23;
tr126:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st95;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
#line 1916 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr131;
		case 34: goto tr128;
	}
	goto st94;
tr131:
#line 1 "NONE"
	{te = p+1;}
#line 321 "wmkdepend.rl"
	{act = 6;}
	goto st182;
st182:
	if ( ++p == pe )
		goto _test_eof182;
case 182:
#line 1932 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto st14;
		case 34: goto tr133;
		case 39: goto tr134;
		case 92: goto st97;
	}
	goto st96;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	switch( (*p) ) {
		case 10: goto st14;
		case 34: goto tr133;
		case 39: goto tr134;
		case 92: goto st97;
	}
	goto st96;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	if ( (*p) == 34 )
		goto tr133;
	goto st96;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
	switch( (*p) ) {
		case 10: goto tr110;
		case 32: goto st80;
		case 35: goto st84;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
		case 110: goto st100;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
		case 99: goto st101;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
		case 108: goto st102;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
		case 117: goto st103;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
		case 100: goto st104;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
		case 101: goto st105;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	switch( (*p) ) {
		case 10: goto tr27;
		case 32: goto st92;
		case 34: goto st106;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st92;
	goto st79;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
	switch( (*p) ) {
		case 10: goto tr30;
		case 32: goto tr124;
		case 34: goto st79;
		case 35: goto tr144;
		case 39: goto tr145;
		case 92: goto tr146;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr124;
	goto tr143;
tr143:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st107;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
#line 2098 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
tr148:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st183;
st183:
	if ( ++p == pe )
		goto _test_eof183;
case 183:
#line 2122 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st80;
		case 35: goto st84;
		case 39: goto tr106;
		case 92: goto st98;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st79;
tr144:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st108;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
#line 2141 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr47;
		case 32: goto st109;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
		case 105: goto st118;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st109;
	goto st107;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
	switch( (*p) ) {
		case 10: goto tr47;
		case 32: goto st109;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
		case 105: goto st110;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st109;
	goto st94;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
		case 110: goto st111;
	}
	goto st94;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
		case 99: goto st112;
	}
	goto st94;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
		case 108: goto st113;
	}
	goto st94;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
		case 117: goto st114;
	}
	goto st94;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
		case 100: goto st115;
	}
	goto st94;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
	switch( (*p) ) {
		case 10: goto tr32;
		case 34: goto tr128;
		case 39: goto tr129;
		case 92: goto st95;
		case 101: goto st116;
	}
	goto st94;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
	switch( (*p) ) {
		case 10: goto tr65;
		case 32: goto st116;
		case 34: goto tr161;
		case 39: goto tr129;
		case 92: goto st95;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st116;
	goto st94;
tr161:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st184;
st184:
	if ( ++p == pe )
		goto _test_eof184;
case 184:
#line 2267 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr30;
		case 34: goto st80;
		case 39: goto tr125;
		case 92: goto tr126;
	}
	goto tr124;
tr150:
#line 1 "NONE"
	{te = p+1;}
#line 316 "wmkdepend.rl"
	{act = 2;}
	goto st185;
tr145:
#line 1 "NONE"
	{te = p+1;}
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
#line 316 "wmkdepend.rl"
	{act = 2;}
	goto st185;
st185:
	if ( ++p == pe )
		goto _test_eof185;
case 185:
#line 2293 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st23;
		case 34: goto tr44;
		case 35: goto st33;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st23;
	goto st32;
tr146:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st117;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
#line 2311 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr131;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
		case 110: goto st119;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
		case 99: goto st120;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
		case 108: goto st121;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
		case 117: goto st122;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
		case 100: goto st123;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	switch( (*p) ) {
		case 10: goto tr32;
		case 32: goto st94;
		case 34: goto tr148;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
		case 101: goto st124;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st94;
	goto st107;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
	switch( (*p) ) {
		case 10: goto tr65;
		case 32: goto st116;
		case 34: goto tr168;
		case 35: goto st108;
		case 39: goto tr150;
		case 92: goto st117;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st116;
	goto st107;
tr168:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st186;
st186:
	if ( ++p == pe )
		goto _test_eof186;
case 186:
#line 2444 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr30;
		case 32: goto tr124;
		case 34: goto st79;
		case 35: goto tr144;
		case 39: goto tr145;
		case 92: goto tr146;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr124;
	goto tr143;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 35: goto st3;
		case 42: goto tr169;
		case 47: goto st126;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
	if ( (*p) == 10 )
		goto tr172;
	goto st127;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
	switch( (*p) ) {
		case 10: goto tr175;
		case 32: goto st129;
		case 35: goto st128;
		case 105: goto st139;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st129;
	goto st126;
st129:
	if ( ++p == pe )
		goto _test_eof129;
case 129:
	switch( (*p) ) {
		case 10: goto tr175;
		case 32: goto st129;
		case 105: goto st130;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st129;
	goto st127;
st130:
	if ( ++p == pe )
		goto _test_eof130;
case 130:
	switch( (*p) ) {
		case 10: goto tr172;
		case 110: goto st131;
	}
	goto st127;
st131:
	if ( ++p == pe )
		goto _test_eof131;
case 131:
	switch( (*p) ) {
		case 10: goto tr172;
		case 99: goto st132;
	}
	goto st127;
st132:
	if ( ++p == pe )
		goto _test_eof132;
case 132:
	switch( (*p) ) {
		case 10: goto tr172;
		case 108: goto st133;
	}
	goto st127;
st133:
	if ( ++p == pe )
		goto _test_eof133;
case 133:
	switch( (*p) ) {
		case 10: goto tr172;
		case 117: goto st134;
	}
	goto st127;
st134:
	if ( ++p == pe )
		goto _test_eof134;
case 134:
	switch( (*p) ) {
		case 10: goto tr172;
		case 100: goto st135;
	}
	goto st127;
st135:
	if ( ++p == pe )
		goto _test_eof135;
case 135:
	switch( (*p) ) {
		case 10: goto tr172;
		case 101: goto st136;
	}
	goto st127;
st136:
	if ( ++p == pe )
		goto _test_eof136;
case 136:
	switch( (*p) ) {
		case 10: goto tr184;
		case 32: goto st136;
		case 34: goto st137;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st136;
	goto st127;
st137:
	if ( ++p == pe )
		goto _test_eof137;
case 137:
	switch( (*p) ) {
		case 10: goto tr187;
		case 34: goto st127;
	}
	goto tr186;
tr186:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st138;
st138:
	if ( ++p == pe )
		goto _test_eof138;
case 138:
#line 2597 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr189;
		case 34: goto tr190;
	}
	goto st138;
tr190:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st187;
st187:
	if ( ++p == pe )
		goto _test_eof187;
case 187:
#line 2615 "wmkdepend.cpp"
	if ( (*p) == 10 )
		goto tr172;
	goto st127;
st139:
	if ( ++p == pe )
		goto _test_eof139;
case 139:
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
		case 110: goto st140;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
st140:
	if ( ++p == pe )
		goto _test_eof140;
case 140:
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
		case 99: goto st141;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
st141:
	if ( ++p == pe )
		goto _test_eof141;
case 141:
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
		case 108: goto st142;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
st142:
	if ( ++p == pe )
		goto _test_eof142;
case 142:
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
		case 117: goto st143;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
st143:
	if ( ++p == pe )
		goto _test_eof143;
case 143:
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
		case 100: goto st144;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
st144:
	if ( ++p == pe )
		goto _test_eof144;
case 144:
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
		case 101: goto st145;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
st145:
	if ( ++p == pe )
		goto _test_eof145;
case 145:
	switch( (*p) ) {
		case 10: goto tr184;
		case 32: goto st136;
		case 34: goto st146;
		case 35: goto st128;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st136;
	goto st126;
st146:
	if ( ++p == pe )
		goto _test_eof146;
case 146:
	switch( (*p) ) {
		case 10: goto tr187;
		case 32: goto tr186;
		case 34: goto st126;
		case 35: goto tr199;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr186;
	goto tr198;
tr198:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st147;
st147:
	if ( ++p == pe )
		goto _test_eof147;
case 147:
#line 2731 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto st138;
		case 34: goto tr201;
		case 35: goto st148;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st138;
	goto st147;
tr201:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st188;
st188:
	if ( ++p == pe )
		goto _test_eof188;
case 188:
#line 2753 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr172;
		case 32: goto st127;
		case 35: goto st128;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st127;
	goto st126;
tr199:
#line 306 "wmkdepend.rl"
	{ tok = p; /* Local token start */ }
	goto st148;
st148:
	if ( ++p == pe )
		goto _test_eof148;
case 148:
#line 2770 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr204;
		case 32: goto st149;
		case 34: goto tr201;
		case 35: goto st148;
		case 105: goto st157;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st149;
	goto st147;
st149:
	if ( ++p == pe )
		goto _test_eof149;
case 149:
	switch( (*p) ) {
		case 10: goto tr204;
		case 32: goto st149;
		case 34: goto tr190;
		case 105: goto st150;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st149;
	goto st138;
st150:
	if ( ++p == pe )
		goto _test_eof150;
case 150:
	switch( (*p) ) {
		case 10: goto tr189;
		case 34: goto tr190;
		case 110: goto st151;
	}
	goto st138;
st151:
	if ( ++p == pe )
		goto _test_eof151;
case 151:
	switch( (*p) ) {
		case 10: goto tr189;
		case 34: goto tr190;
		case 99: goto st152;
	}
	goto st138;
st152:
	if ( ++p == pe )
		goto _test_eof152;
case 152:
	switch( (*p) ) {
		case 10: goto tr189;
		case 34: goto tr190;
		case 108: goto st153;
	}
	goto st138;
st153:
	if ( ++p == pe )
		goto _test_eof153;
case 153:
	switch( (*p) ) {
		case 10: goto tr189;
		case 34: goto tr190;
		case 117: goto st154;
	}
	goto st138;
st154:
	if ( ++p == pe )
		goto _test_eof154;
case 154:
	switch( (*p) ) {
		case 10: goto tr189;
		case 34: goto tr190;
		case 100: goto st155;
	}
	goto st138;
st155:
	if ( ++p == pe )
		goto _test_eof155;
case 155:
	switch( (*p) ) {
		case 10: goto tr189;
		case 34: goto tr190;
		case 101: goto st156;
	}
	goto st138;
st156:
	if ( ++p == pe )
		goto _test_eof156;
case 156:
	switch( (*p) ) {
		case 10: goto tr213;
		case 32: goto st156;
		case 34: goto tr214;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st156;
	goto st138;
tr214:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st189;
st189:
	if ( ++p == pe )
		goto _test_eof189;
case 189:
#line 2878 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr187;
		case 34: goto st127;
	}
	goto tr186;
st157:
	if ( ++p == pe )
		goto _test_eof157;
case 157:
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto st138;
		case 34: goto tr201;
		case 35: goto st148;
		case 110: goto st158;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st138;
	goto st147;
st158:
	if ( ++p == pe )
		goto _test_eof158;
case 158:
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto st138;
		case 34: goto tr201;
		case 35: goto st148;
		case 99: goto st159;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st138;
	goto st147;
st159:
	if ( ++p == pe )
		goto _test_eof159;
case 159:
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto st138;
		case 34: goto tr201;
		case 35: goto st148;
		case 108: goto st160;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st138;
	goto st147;
st160:
	if ( ++p == pe )
		goto _test_eof160;
case 160:
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto st138;
		case 34: goto tr201;
		case 35: goto st148;
		case 117: goto st161;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st138;
	goto st147;
st161:
	if ( ++p == pe )
		goto _test_eof161;
case 161:
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto st138;
		case 34: goto tr201;
		case 35: goto st148;
		case 100: goto st162;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st138;
	goto st147;
st162:
	if ( ++p == pe )
		goto _test_eof162;
case 162:
	switch( (*p) ) {
		case 10: goto tr189;
		case 32: goto st138;
		case 34: goto tr201;
		case 35: goto st148;
		case 101: goto st163;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st138;
	goto st147;
st163:
	if ( ++p == pe )
		goto _test_eof163;
case 163:
	switch( (*p) ) {
		case 10: goto tr213;
		case 32: goto st156;
		case 34: goto tr221;
		case 35: goto st148;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st156;
	goto st147;
tr221:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{ processFile(std::string(tok, (p - tok))); }
#line 313 "wmkdepend.rl"
	{act = 1;}
	goto st190;
st190:
	if ( ++p == pe )
		goto _test_eof190;
case 190:
#line 2993 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr187;
		case 32: goto tr186;
		case 34: goto st126;
		case 35: goto tr199;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr186;
	goto tr198;
st164:
	if ( ++p == pe )
		goto _test_eof164;
case 164:
	switch( (*p) ) {
		case 10: goto tr3;
		case 32: goto st2;
		case 34: goto st57;
		case 35: goto st3;
		case 39: goto st79;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st165:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof165;
case 165:
#line 3023 "wmkdepend.cpp"
	if ( (*p) == 42 )
		goto st166;
	goto st165;
st166:
	if ( ++p == pe )
		goto _test_eof166;
case 166:
	switch( (*p) ) {
		case 42: goto st166;
		case 47: goto tr224;
	}
	goto st165;
tr224:
#line 309 "wmkdepend.rl"
	{ {goto st167;} }
	goto st191;
st191:
	if ( ++p == pe )
		goto _test_eof191;
case 191:
#line 3044 "wmkdepend.cpp"
	goto st0;
st0:
cs = 0;
	goto _out;
	}
	_test_eof167: cs = 167; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof168: cs = 168; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof169: cs = 169; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof170: cs = 170; goto _test_eof; 
	_test_eof171: cs = 171; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof172: cs = 172; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof173: cs = 173; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof174: cs = 174; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof175: cs = 175; goto _test_eof; 
	_test_eof176: cs = 176; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof177: cs = 177; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof178: cs = 178; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof179: cs = 179; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof180: cs = 180; goto _test_eof; 
	_test_eof181: cs = 181; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof182: cs = 182; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof183: cs = 183; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof184: cs = 184; goto _test_eof; 
	_test_eof185: cs = 185; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof186: cs = 186; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof129: cs = 129; goto _test_eof; 
	_test_eof130: cs = 130; goto _test_eof; 
	_test_eof131: cs = 131; goto _test_eof; 
	_test_eof132: cs = 132; goto _test_eof; 
	_test_eof133: cs = 133; goto _test_eof; 
	_test_eof134: cs = 134; goto _test_eof; 
	_test_eof135: cs = 135; goto _test_eof; 
	_test_eof136: cs = 136; goto _test_eof; 
	_test_eof137: cs = 137; goto _test_eof; 
	_test_eof138: cs = 138; goto _test_eof; 
	_test_eof187: cs = 187; goto _test_eof; 
	_test_eof139: cs = 139; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
	_test_eof142: cs = 142; goto _test_eof; 
	_test_eof143: cs = 143; goto _test_eof; 
	_test_eof144: cs = 144; goto _test_eof; 
	_test_eof145: cs = 145; goto _test_eof; 
	_test_eof146: cs = 146; goto _test_eof; 
	_test_eof147: cs = 147; goto _test_eof; 
	_test_eof188: cs = 188; goto _test_eof; 
	_test_eof148: cs = 148; goto _test_eof; 
	_test_eof149: cs = 149; goto _test_eof; 
	_test_eof150: cs = 150; goto _test_eof; 
	_test_eof151: cs = 151; goto _test_eof; 
	_test_eof152: cs = 152; goto _test_eof; 
	_test_eof153: cs = 153; goto _test_eof; 
	_test_eof154: cs = 154; goto _test_eof; 
	_test_eof155: cs = 155; goto _test_eof; 
	_test_eof156: cs = 156; goto _test_eof; 
	_test_eof189: cs = 189; goto _test_eof; 
	_test_eof157: cs = 157; goto _test_eof; 
	_test_eof158: cs = 158; goto _test_eof; 
	_test_eof159: cs = 159; goto _test_eof; 
	_test_eof160: cs = 160; goto _test_eof; 
	_test_eof161: cs = 161; goto _test_eof; 
	_test_eof162: cs = 162; goto _test_eof; 
	_test_eof163: cs = 163; goto _test_eof; 
	_test_eof190: cs = 190; goto _test_eof; 
	_test_eof164: cs = 164; goto _test_eof; 
	_test_eof165: cs = 165; goto _test_eof; 
	_test_eof166: cs = 166; goto _test_eof; 
	_test_eof191: cs = 191; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 1: goto tr0;
	case 2: goto tr0;
	case 3: goto tr0;
	case 4: goto tr0;
	case 168: goto tr0;
	case 5: goto tr0;
	case 6: goto tr0;
	case 7: goto tr0;
	case 8: goto tr0;
	case 9: goto tr0;
	case 10: goto tr0;
	case 11: goto tr0;
	case 12: goto tr0;
	case 13: goto tr0;
	case 14: goto tr0;
	case 15: goto tr0;
	case 16: goto tr0;
	case 17: goto tr0;
	case 18: goto tr0;
	case 19: goto tr0;
	case 20: goto tr0;
	case 21: goto tr0;
	case 169: goto tr0;
	case 22: goto tr0;
	case 23: goto tr0;
	case 170: goto tr0;
	case 171: goto tr0;
	case 24: goto tr0;
	case 25: goto tr0;
	case 26: goto tr0;
	case 27: goto tr0;
	case 28: goto tr0;
	case 29: goto tr0;
	case 30: goto tr0;
	case 31: goto tr0;
	case 32: goto tr0;
	case 172: goto tr0;
	case 33: goto tr0;
	case 34: goto tr0;
	case 173: goto tr0;
	case 35: goto tr0;
	case 36: goto tr0;
	case 37: goto tr0;
	case 38: goto tr0;
	case 39: goto tr0;
	case 40: goto tr0;
	case 41: goto tr0;
	case 42: goto tr0;
	case 174: goto tr227;
	case 43: goto tr0;
	case 44: goto tr0;
	case 45: goto tr0;
	case 46: goto tr0;
	case 47: goto tr0;
	case 48: goto tr0;
	case 49: goto tr0;
	case 175: goto tr0;
	case 176: goto tr0;
	case 50: goto tr0;
	case 51: goto tr0;
	case 52: goto tr0;
	case 53: goto tr0;
	case 54: goto tr0;
	case 55: goto tr0;
	case 56: goto tr0;
	case 177: goto tr0;
	case 178: goto tr228;
	case 60: goto tr82;
	case 61: goto tr82;
	case 79: goto tr0;
	case 80: goto tr0;
	case 81: goto tr0;
	case 179: goto tr0;
	case 82: goto tr0;
	case 83: goto tr0;
	case 84: goto tr0;
	case 85: goto tr0;
	case 86: goto tr0;
	case 87: goto tr0;
	case 88: goto tr0;
	case 89: goto tr0;
	case 90: goto tr0;
	case 91: goto tr0;
	case 92: goto tr0;
	case 93: goto tr0;
	case 94: goto tr0;
	case 180: goto tr227;
	case 181: goto tr229;
	case 95: goto tr0;
	case 182: goto tr228;
	case 96: goto tr82;
	case 97: goto tr82;
	case 98: goto tr0;
	case 99: goto tr0;
	case 100: goto tr0;
	case 101: goto tr0;
	case 102: goto tr0;
	case 103: goto tr0;
	case 104: goto tr0;
	case 105: goto tr0;
	case 106: goto tr0;
	case 107: goto tr0;
	case 183: goto tr227;
	case 108: goto tr0;
	case 109: goto tr0;
	case 110: goto tr0;
	case 111: goto tr0;
	case 112: goto tr0;
	case 113: goto tr0;
	case 114: goto tr0;
	case 115: goto tr0;
	case 116: goto tr0;
	case 184: goto tr227;
	case 185: goto tr229;
	case 117: goto tr0;
	case 118: goto tr0;
	case 119: goto tr0;
	case 120: goto tr0;
	case 121: goto tr0;
	case 122: goto tr0;
	case 123: goto tr0;
	case 124: goto tr0;
	case 186: goto tr227;
	case 126: goto tr0;
	case 127: goto tr0;
	case 128: goto tr0;
	case 129: goto tr0;
	case 130: goto tr0;
	case 131: goto tr0;
	case 132: goto tr0;
	case 133: goto tr0;
	case 134: goto tr0;
	case 135: goto tr0;
	case 136: goto tr0;
	case 137: goto tr0;
	case 138: goto tr0;
	case 187: goto tr227;
	case 139: goto tr0;
	case 140: goto tr0;
	case 141: goto tr0;
	case 142: goto tr0;
	case 143: goto tr0;
	case 144: goto tr0;
	case 145: goto tr0;
	case 146: goto tr0;
	case 147: goto tr0;
	case 188: goto tr227;
	case 148: goto tr0;
	case 149: goto tr0;
	case 150: goto tr0;
	case 151: goto tr0;
	case 152: goto tr0;
	case 153: goto tr0;
	case 154: goto tr0;
	case 155: goto tr0;
	case 156: goto tr0;
	case 189: goto tr227;
	case 157: goto tr0;
	case 158: goto tr0;
	case 159: goto tr0;
	case 160: goto tr0;
	case 161: goto tr0;
	case 162: goto tr0;
	case 163: goto tr0;
	case 190: goto tr227;
	}
	}

	_out: {}
	}

#line 397 "wmkdepend.rl"
       /* ^^^ FSM execution here ^^^ */;

        if (0 == cs)
        {
            // FSM failed before finding a token
            std::cerr
                << EXENAME ": parse error while scanning '"
                << fileName << "'\n";
            break;
        }

        if (ts)
        {
            // Preserve incomplete token
            pending = pe - ts;
            memmove(inbuf, ts, pending);
            te = inbuf + (te - ts);   // token end (after memmove)
            ts = inbuf;               // token start
        }
        else
        {
            pending = 0;
        }
    }
    fclose(infile);
}


int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr
            << EXENAME ": input file not supplied\n";
        return 1;
    }

    unsigned nIncDirs = 0;

    // Prechecks:
    // - help
    // - count number of -I directories
    for (int i = 1; i < argc; ++i)
    {
        if (argv[i][0] != '-') continue;

        switch (argv[i][1])
        {
            case 'h':           // Option: -h, -help
                usage();
                return 0;
                break;

            case 'q':           // Option: -q (quiet)
                optQuiet = true;
                break;

            case 'v':           // Option: -v (verbose)
                optVerbose = true;
                break;

            case 'I':           // Option: -Idir
                ++nIncDirs;
                break;

            // Could check other options, warn about unknown options...
        }
    }

    sourceFile.assign(argv[argc-1]);
    std::string outputFile;

    // Verify that input file has an extension
    {
        auto dot = sourceFile.find_last_of("./");
        if (dot == std::string::npos || sourceFile[dot] != '.')
        {
            std::cerr
                << EXENAME ": no file extension in source file name '"
                << sourceFile << "'\n";

            return 1;
        }
    }

    Files::reset(nIncDirs);

    // Build list of -I directories and add -i ignores
    for (int i = 1; i < argc; ++i)
    {
        const size_t optLen = strlen(argv[i]);

        if (!strncmp(argv[i], "-I", 2))
        {
            // Option: -Idir
            if (optLen > 2)
            {
                Files::dirs.emplace_back(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-i", 2))
        {
            // Option: -iheader
            if (optLen > 2)
            {
                Files::visited.insert(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-e", 2))
        {
            // Option: -eENV
            if (optLen > 2)
            {
                Files::addEnv(argv[i] + 2);
            }
        }
        else if (!strncmp(argv[i], "-o", 2))
        {
            // Option: -oFile */
            if (optLen > 2)
            {
                outputFile.assign(argv[i] + 2);
            }
        }
    }

    if (outputFile.size() && !freopen(outputFile.c_str(), "w", stdout))
    {
        std::cerr
            << EXENAME ": could not open file '"
            << outputFile << "' for output: " << strerror(errno) << "\n";
        return 1;
    }

    fputs("$(OBJECTS_DIR)/", stdout);
    fputs(sourceFile.c_str(), stdout);
    fputs(".dep: \\\n", stdout);

    processFile(sourceFile);

    fputs("\n\n", stdout);
    fflush(stdout);

    return 0;
}


/*****************************************************************************/
