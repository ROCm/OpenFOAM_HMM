
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
        "\nDependency list generator, similar to 'cpp -M'\n\n";
}

// Suppress some error messages
bool optQuiet = false;

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


    //- Add environ replacements
    //
    // Eg,
    //     /openfoam/project/path/directory/xyz
    //  -> $(WM_PROJECT_DIR)/directory/xyz
    void addEnv(std::string key)
    {
        const char *val = ::getenv(key.c_str());

        if (val && *val)
        {
            // "$(ENV)/"  -> "/env/value/"
            std::string orig(val);
            if (orig.back() != '/')
            {
                orig.append("/");
            }

            envlist.emplace_front("$(" + key + ")/", std::move(orig));
        }
    }


    //
    // Open a file for reading and emit its qualified name to stdout.
    // Uses env substitutions at the beginning of the path
    //
    // Eg,
    //     /openfoam/project/path/directory/xyz
    //  -> $(WM_PROJECT_DIR)/directory/xyz
    //
    FILE* fopen_file(const std::string& fileName)
    {
        const auto len = fileName.size();
        const char *fname = fileName.c_str();

        FILE *filePtr = ::fopen(fname, "r");

        if (filePtr)
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
                    fname += entry.len;
                    ::fputs(entry.name.c_str(), stdout);
                    break;
                }
            }

            ::fputs(fname, stdout);
            ::fputs(" \\\n", stdout);
        }
        else if (errno == EMFILE)
        {
            std::cerr
                << EXENAME ": too many open files while opening '"
                << fileName << "'\n"
                << "Please change your open descriptor limit\n";
        }

        return filePtr;
    }


    // Open a not previously visited file for reading, using the include dirs
    // as required.
    //
    // On success, emits the resolved name on stdout
    //
    FILE* open(const std::string& fileName)
    {
        // Bad file name, or already visited
        if (fileName.empty() || visited.find(fileName) != visited.end())
        {
            return nullptr;
        }

        FILE* filePtr = fopen_file(fileName);
        if (!filePtr)
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

                filePtr = fopen_file(fullName);

                if (filePtr)
                {
                    break;
                }
            }
        }

        // Mark as having been visited.
        // This also avoids re-testing and multiple messages.
        visited.insert(fileName);

        // Report failues
        if (!filePtr && !optQuiet)
        {
            std::cerr
                << EXENAME ": could not open file '"
                << fileName << "' for source file '"
                << sourceFile << "'";

            if (dirs.size())
            {
                std::cerr << ": " << strerror(errno);
            }

            std::cerr << '\n' << std::flush;
        }

        return filePtr;
    }

} // end of namespace Files


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine requirements: token start/end, action, code state
// are defined later (prior to use)


#line 290 "wmkdepend.cpp"
static const int wmkdep_start = 37;
static const int wmkdep_error = 0;

static const int wmkdep_en_consume_comment = 35;
static const int wmkdep_en_main = 37;


#line 317 "wmkdepend.rl"



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//
// Open a file and process
//
void processFile(const std::string& fileName)
{
    FILE* infile = Files::open(fileName);
    if (!infile) return;

    // ts, te  = Ragel token start/end points (required naming)
    // act, cs = Ragel action, code state, respectively (required naming)
    char *ts, *te;
    unsigned act, cs;

    
#line 318 "wmkdepend.cpp"
	{
	cs = wmkdep_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 337 "wmkdepend.rl"
    // Token start of include filename (for begInclude, endInclude actions)
    char *ts_inclName = nullptr;

    // Buffering
    char inbuf[READ_BUFLEN];
    size_t pending = 0;

    // Processing loop (as per Ragel pdf example)
    for (bool good = true; good; /*nil*/)
    {
        const size_t avail = READ_BUFLEN - pending;

        if (!avail)
        {
            // We overfilled the buffer while trying to scan a token...
            std::cerr
                << EXENAME ": buffer full while scanning '"
                << fileName << "'\n";
            break;
        }

        // p, pe  = Ragel parsing point and parsing end (required naming)
        // eof    = Ragel EOF point (required naming)

        char *p = inbuf + pending;
        const size_t gcount = ::fread(p, 1, avail, infile);

        char *pe = p + gcount;
        char *eof = nullptr;

        if (!gcount)    // Could also use feof(infile)
        {
            // Tag 'pe' as being the EOF for the FSM as well
            eof = pe;
            good = false;
        }

        
#line 365 "wmkdepend.cpp"
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
	{{p = ((te))-1;} {goto st35;} }
	break;
	default:
	{{p = ((te))-1;}}
	break;
	}
	}
	goto st37;
tr2:
#line 314 "wmkdepend.rl"
	{te = p+1;}
	goto st37;
tr6:
#line 314 "wmkdepend.rl"
	{{p = ((te))-1;}}
	goto st37;
tr19:
#line 297 "wmkdepend.rl"
	{
        // std::cerr << std::string(ts, (te-ts)) << '\n';
        processFile(std::string(ts_inclName, (p - ts_inclName)));
    }
#line 310 "wmkdepend.rl"
	{te = p+1;}
	goto st37;
tr40:
#line 308 "wmkdepend.rl"
	{te = p+1;}
	goto st37;
tr47:
#line 307 "wmkdepend.rl"
	{te = p+1;}
	goto st37;
tr51:
#line 313 "wmkdepend.rl"
	{te = p+1;}
	goto st37;
tr57:
#line 314 "wmkdepend.rl"
	{te = p;p--;}
	goto st37;
st37:
#line 1 "NONE"
	{ts = 0;}
#line 1 "NONE"
	{act = 0;}
	if ( ++p == pe )
		goto _test_eof37;
case 37:
#line 1 "NONE"
	{ts = p;}
#line 429 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr4;
		case 32: goto st2;
		case 34: goto st24;
		case 35: goto st14;
		case 39: goto st28;
		case 47: goto st32;
		case 76: goto st34;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( (*p) == 10 )
		goto tr2;
	goto st1;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	switch( (*p) ) {
		case 10: goto tr4;
		case 32: goto st2;
		case 35: goto st14;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st2;
	goto st1;
tr4:
#line 1 "NONE"
	{te = p+1;}
	goto st38;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
#line 469 "wmkdepend.cpp"
	switch( (*p) ) {
		case 32: goto st3;
		case 35: goto st4;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st3;
	goto tr57;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	switch( (*p) ) {
		case 32: goto st3;
		case 35: goto st4;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st3;
	goto tr6;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 32: goto st4;
		case 105: goto st5;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st4;
	goto tr6;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) == 110 )
		goto st6;
	goto tr6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( (*p) == 99 )
		goto st7;
	goto tr6;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 108 )
		goto st8;
	goto tr6;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 117 )
		goto st9;
	goto tr6;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 100 )
		goto st10;
	goto tr6;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 101 )
		goto st11;
	goto tr6;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	switch( (*p) ) {
		case 32: goto st11;
		case 34: goto st12;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st11;
	goto tr6;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 34: goto tr6;
		case 60: goto tr6;
		case 62: goto tr6;
	}
	goto tr17;
tr17:
#line 293 "wmkdepend.rl"
	{ ts_inclName = p; }
	goto st13;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
#line 570 "wmkdepend.cpp"
	switch( (*p) ) {
		case 34: goto tr19;
		case 60: goto tr6;
		case 62: goto tr6;
	}
	goto st13;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	switch( (*p) ) {
		case 10: goto tr20;
		case 32: goto st14;
		case 105: goto st15;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st14;
	goto st1;
tr20:
#line 1 "NONE"
	{te = p+1;}
	goto st39;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
#line 597 "wmkdepend.cpp"
	switch( (*p) ) {
		case 32: goto st4;
		case 105: goto st5;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st4;
	goto tr57;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	switch( (*p) ) {
		case 10: goto tr2;
		case 110: goto st16;
	}
	goto st1;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	switch( (*p) ) {
		case 10: goto tr2;
		case 99: goto st17;
	}
	goto st1;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 10: goto tr2;
		case 108: goto st18;
	}
	goto st1;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	switch( (*p) ) {
		case 10: goto tr2;
		case 117: goto st19;
	}
	goto st1;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	switch( (*p) ) {
		case 10: goto tr2;
		case 100: goto st20;
	}
	goto st1;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	switch( (*p) ) {
		case 10: goto tr2;
		case 101: goto st21;
	}
	goto st1;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	switch( (*p) ) {
		case 10: goto tr28;
		case 32: goto st21;
		case 34: goto st22;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st21;
	goto st1;
tr28:
#line 1 "NONE"
	{te = p+1;}
	goto st40;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
#line 679 "wmkdepend.cpp"
	switch( (*p) ) {
		case 32: goto st11;
		case 34: goto st12;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st11;
	goto tr57;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	switch( (*p) ) {
		case 10: goto tr31;
		case 34: goto st1;
		case 60: goto st1;
		case 62: goto st1;
	}
	goto tr30;
tr30:
#line 293 "wmkdepend.rl"
	{ ts_inclName = p; }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 706 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr33;
		case 34: goto tr34;
		case 60: goto st1;
		case 62: goto st1;
	}
	goto st23;
tr33:
#line 1 "NONE"
	{te = p+1;}
	goto st41;
tr31:
#line 1 "NONE"
	{te = p+1;}
#line 293 "wmkdepend.rl"
	{ ts_inclName = p; }
	goto st41;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
#line 728 "wmkdepend.cpp"
	switch( (*p) ) {
		case 34: goto tr19;
		case 60: goto tr57;
		case 62: goto tr57;
	}
	goto st13;
tr34:
#line 1 "NONE"
	{te = p+1;}
#line 297 "wmkdepend.rl"
	{
        // std::cerr << std::string(ts, (te-ts)) << '\n';
        processFile(std::string(ts_inclName, (p - ts_inclName)));
    }
#line 310 "wmkdepend.rl"
	{act = 3;}
	goto st42;
tr36:
#line 1 "NONE"
	{te = p+1;}
#line 308 "wmkdepend.rl"
	{act = 2;}
	goto st42;
tr43:
#line 1 "NONE"
	{te = p+1;}
#line 307 "wmkdepend.rl"
	{act = 1;}
	goto st42;
tr49:
#line 1 "NONE"
	{te = p+1;}
#line 312 "wmkdepend.rl"
	{act = 4;}
	goto st42;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
#line 768 "wmkdepend.cpp"
	if ( (*p) == 10 )
		goto tr2;
	goto st1;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 10: goto tr2;
		case 34: goto tr36;
		case 92: goto st25;
	}
	goto st24;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	if ( (*p) == 10 )
		goto tr38;
	goto st24;
tr38:
#line 1 "NONE"
	{te = p+1;}
	goto st43;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
#line 797 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr57;
		case 34: goto tr40;
		case 92: goto st27;
	}
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 10: goto tr6;
		case 34: goto tr40;
		case 92: goto st27;
	}
	goto st26;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	goto st26;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 10: goto tr2;
		case 39: goto tr43;
		case 92: goto st29;
	}
	goto st28;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	if ( (*p) == 10 )
		goto tr45;
	goto st28;
tr45:
#line 1 "NONE"
	{te = p+1;}
	goto st44;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
#line 844 "wmkdepend.cpp"
	switch( (*p) ) {
		case 10: goto tr57;
		case 39: goto tr47;
		case 92: goto st31;
	}
	goto st30;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 10: goto tr6;
		case 39: goto tr47;
		case 92: goto st31;
	}
	goto st30;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	goto st30;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 10: goto tr2;
		case 42: goto tr49;
		case 47: goto st33;
	}
	goto st1;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 10 )
		goto tr51;
	goto st33;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 10: goto tr2;
		case 34: goto st24;
		case 39: goto st28;
	}
	goto st1;
st35:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof35;
case 35:
#line 899 "wmkdepend.cpp"
	if ( (*p) == 42 )
		goto st36;
	goto st35;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 42: goto st36;
		case 47: goto tr54;
	}
	goto st35;
tr54:
#line 302 "wmkdepend.rl"
	{ {goto st37;} }
	goto st45;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
#line 920 "wmkdepend.cpp"
	goto st0;
st0:
cs = 0;
	goto _out;
	}
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
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
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 1: goto tr0;
	case 38: goto tr57;
	case 3: goto tr6;
	case 4: goto tr6;
	case 5: goto tr6;
	case 6: goto tr6;
	case 7: goto tr6;
	case 8: goto tr6;
	case 9: goto tr6;
	case 10: goto tr6;
	case 11: goto tr6;
	case 12: goto tr6;
	case 13: goto tr6;
	case 39: goto tr57;
	case 40: goto tr57;
	case 41: goto tr57;
	case 42: goto tr0;
	case 43: goto tr57;
	case 26: goto tr6;
	case 27: goto tr6;
	case 44: goto tr57;
	case 30: goto tr6;
	case 31: goto tr6;
	}
	}

	_out: {}
	}

#line 376 "wmkdepend.rl"
        if (cs == wmkdep_error)
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
            ::memmove(inbuf, ts, pending);
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
                << EXENAME ": cannot find extension in source file name '"
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


    // Start of output
    if (outputFile.size())
    {
        FILE *reopened = freopen(outputFile.c_str(), "w", stdout);
        if (!reopened)
        {
            std::cerr
                << EXENAME ": could not open file '"
                << outputFile << "' for output: " << strerror(errno)
                << "\n";
            return 1;
        }
    }

    ::fputs("$(OBJECTS_DIR)/", stdout);
    ::fputs(sourceFile.c_str(), stdout);
    ::fputs(".dep: \\\n", stdout);

    processFile(sourceFile);

    ::fputs("\n\n", stdout);
    ::fflush(stdout);

    return 0;
}


/*****************************************************************************/
