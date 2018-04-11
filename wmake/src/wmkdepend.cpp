
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
    A fast dependency list generator that emulates the behaviour and output
    of cpp -M. However, the output contains no duplicates and is thus
    approx. 40% faster than cpp.
    It also handles missing files somewhat more gracefully.

    The algorithm uses a Ragel-generated lexer to scan for includes and
    searches the files found.
    The files are only visited once (their names are hashed),
    which helps make this faster than cpp.

Usage
    wmkdepend [-Idir..] [-iheader...] [-eENV...] [-oFile] [-q] filename

\*---------------------------------------------------------------------------*/
/*
 * With cpp:
 *
 * cpp -x c++ -std=c++11 -nostdinc -nostdinc++
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

// Sizing for read buffer
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
        "  -v        Some verbosity.\n"
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
    FILE* fopen_file(std::string fileName)
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

            std::cerr << "\n" << std::flush;
        }

        return filePtr;
    }

} // end of namespace Files


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine requirements: text start/end, action, code state
// are defined later (prior to use)


#line 283 "wmkdepend.cpp"
static const char _scanInclude_actions[] = {
	0, 1, 0, 1, 2, 1, 3, 1, 
	5, 1, 6, 1, 11, 1, 12, 1, 
	14, 1, 15, 1, 16, 1, 17, 1, 
	18, 2, 1, 13, 2, 3, 4, 2, 
	6, 0, 2, 6, 7, 2, 6, 8, 
	2, 6, 10, 3, 6, 1, 9
};

static const char _scanInclude_key_offsets[] = {
	0, 0, 1, 4, 5, 8, 8, 13, 
	17, 18, 19, 20, 21, 22, 23, 27, 
	28, 29, 31, 33, 35, 37, 39, 41, 
	46, 48, 50, 53, 54, 57, 57, 60, 
	61, 64, 65, 67, 73, 74, 77, 81, 
	85, 86, 89
};

static const char _scanInclude_trans_keys[] = {
	10, 10, 34, 92, 10, 10, 34, 92, 
	10, 32, 105, 9, 13, 32, 105, 9, 
	13, 110, 99, 108, 117, 100, 101, 32, 
	34, 9, 13, 34, 34, 10, 110, 10, 
	99, 10, 108, 10, 117, 10, 100, 10, 
	101, 10, 32, 34, 9, 13, 10, 34, 
	10, 34, 10, 39, 92, 10, 10, 39, 
	92, 10, 42, 47, 10, 10, 34, 39, 
	42, 42, 47, 10, 34, 35, 39, 47, 
	76, 10, 10, 34, 92, 32, 105, 9, 
	13, 32, 34, 9, 13, 34, 10, 39, 
	92, 0
};

static const char _scanInclude_single_lengths[] = {
	0, 1, 3, 1, 3, 0, 3, 2, 
	1, 1, 1, 1, 1, 1, 2, 1, 
	1, 2, 2, 2, 2, 2, 2, 3, 
	2, 2, 3, 1, 3, 0, 3, 1, 
	3, 1, 2, 6, 1, 3, 2, 2, 
	1, 3, 0
};

static const char _scanInclude_range_lengths[] = {
	0, 0, 0, 0, 0, 0, 1, 1, 
	0, 0, 0, 0, 0, 0, 1, 0, 
	0, 0, 0, 0, 0, 0, 0, 1, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 1, 1, 
	0, 0, 0
};

static const unsigned char _scanInclude_index_offsets[] = {
	0, 0, 2, 6, 8, 12, 13, 18, 
	22, 24, 26, 28, 30, 32, 34, 38, 
	40, 42, 45, 48, 51, 54, 57, 60, 
	65, 68, 71, 75, 77, 81, 82, 86, 
	88, 92, 94, 97, 104, 106, 110, 114, 
	118, 120, 124
};

static const char _scanInclude_indicies[] = {
	2, 1, 2, 4, 5, 3, 6, 3, 
	7, 9, 10, 8, 8, 12, 11, 13, 
	11, 1, 14, 15, 14, 7, 16, 7, 
	17, 7, 18, 7, 19, 7, 20, 7, 
	21, 7, 21, 22, 21, 7, 7, 23, 
	25, 24, 2, 26, 1, 2, 27, 1, 
	2, 28, 1, 2, 29, 1, 2, 30, 
	1, 2, 31, 1, 32, 31, 33, 31, 
	1, 35, 1, 34, 37, 38, 36, 2, 
	40, 41, 39, 42, 39, 7, 44, 45, 
	43, 43, 2, 46, 47, 1, 48, 47, 
	2, 3, 39, 1, 50, 49, 50, 51, 
	49, 2, 3, 11, 39, 52, 53, 1, 
	2, 1, 54, 9, 10, 8, 14, 15, 
	14, 54, 21, 22, 21, 54, 25, 24, 
	54, 44, 45, 43, 55, 0
};

static const char _scanInclude_trans_targs[] = {
	35, 1, 35, 2, 36, 3, 37, 35, 
	4, 35, 5, 6, 38, 17, 7, 8, 
	9, 10, 11, 12, 13, 14, 15, 16, 
	16, 35, 18, 19, 20, 21, 22, 23, 
	39, 24, 25, 40, 25, 40, 36, 26, 
	36, 27, 41, 28, 35, 29, 36, 31, 
	35, 33, 34, 42, 30, 32, 35, 0
};

static const char _scanInclude_trans_actions[] = {
	23, 0, 17, 0, 37, 0, 9, 21, 
	0, 13, 0, 0, 9, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 1, 
	0, 25, 0, 0, 0, 0, 0, 0, 
	9, 0, 1, 31, 0, 9, 43, 0, 
	34, 0, 9, 0, 11, 0, 40, 0, 
	15, 0, 0, 3, 0, 0, 19, 0
};

static const char _scanInclude_to_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 5, 0, 28, 0, 0, 0, 0, 
	0, 0, 0
};

static const char _scanInclude_from_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 7, 0, 0, 0, 0, 
	0, 0, 0
};

static const unsigned char _scanInclude_eof_trans[] = {
	0, 1, 0, 0, 8, 8, 0, 8, 
	8, 8, 8, 8, 8, 8, 8, 8, 
	8, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 8, 8, 0, 0, 
	0, 0, 0, 0, 1, 55, 55, 55, 
	55, 55, 0
};

static const int scanInclude_start = 35;
static const int scanInclude_error = 0;

static const int scanInclude_en_consume_comment = 33;
static const int scanInclude_en_main = 35;


#line 313 "wmkdepend.rl"



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//
// Open a file and process
//
void processFile(const std::string& fileName)
{
    FILE* infile = Files::open(fileName);
    if (!infile) return;


    // Ragel text start/end, action, code state
    char *ts = nullptr, *te = nullptr;
    unsigned act, cs;

    // For the include action:
    char *include_start = nullptr;

    
#line 440 "wmkdepend.cpp"
	{
	cs = scanInclude_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 334 "wmkdepend.rl"


    // Buffering
    char buf[READ_BUFLEN];
    size_t bytesPending = 0;

    // Do the first read
    for (bool good = true; good; /*nil*/)
    {
        const size_t avail = READ_BUFLEN - bytesPending;

        if (!avail)
        {
            // We overfilled the buffer while trying to scan a token...
            std::cerr
                << "OUT OF BUFFER SPACE while scanning " << fileName << '\n';
            break;
        }

        char *p = buf + bytesPending;
        const size_t bytesRead = ::fread(p, 1, avail, infile);

        char *pe = p + bytesRead;
        char *eof = nullptr;

        // If we see eof then append the EOF char.
        if (feof(infile))
        {
            eof = pe;
            good = false;
        }

        
#line 482 "wmkdepend.cpp"
	{
	int _klen;
	unsigned int _trans;
	const char *_acts;
	unsigned int _nacts;
	const char *_keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_acts = _scanInclude_actions + _scanInclude_from_state_actions[cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 5:
#line 1 "NONE"
	{ts = p;}
	break;
#line 503 "wmkdepend.cpp"
		}
	}

	_keys = _scanInclude_trans_keys + _scanInclude_key_offsets[cs];
	_trans = _scanInclude_index_offsets[cs];

	_klen = _scanInclude_single_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += (unsigned int)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _scanInclude_range_lengths[cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += (unsigned int)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _scanInclude_indicies[_trans];
_eof_trans:
	cs = _scanInclude_trans_targs[_trans];

	if ( _scanInclude_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _scanInclude_actions + _scanInclude_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 283 "wmkdepend.rl"
	{ include_start = p; }
	break;
	case 1:
#line 286 "wmkdepend.rl"
	{
        if (include_start)
        {
            processFile(std::string(include_start, (p - include_start)));
        }
        include_start = nullptr;
    }
	break;
	case 2:
#line 295 "wmkdepend.rl"
	{ {cs = 35; goto _again;} }
	break;
	case 6:
#line 1 "NONE"
	{te = p+1;}
	break;
	case 7:
#line 303 "wmkdepend.rl"
	{act = 1;}
	break;
	case 8:
#line 304 "wmkdepend.rl"
	{act = 2;}
	break;
	case 9:
#line 306 "wmkdepend.rl"
	{act = 3;}
	break;
	case 10:
#line 308 "wmkdepend.rl"
	{act = 4;}
	break;
	case 11:
#line 303 "wmkdepend.rl"
	{te = p+1;}
	break;
	case 12:
#line 304 "wmkdepend.rl"
	{te = p+1;}
	break;
	case 13:
#line 306 "wmkdepend.rl"
	{te = p+1;}
	break;
	case 14:
#line 309 "wmkdepend.rl"
	{te = p+1;}
	break;
	case 15:
#line 310 "wmkdepend.rl"
	{te = p+1;}
	break;
	case 16:
#line 310 "wmkdepend.rl"
	{te = p;p--;}
	break;
	case 17:
#line 310 "wmkdepend.rl"
	{{p = ((te))-1;}}
	break;
	case 18:
#line 1 "NONE"
	{	switch( act ) {
	case 0:
	{{cs = 0; goto _again;}}
	break;
	case 4:
	{{p = ((te))-1;} {cs = 33; goto _again;} }
	break;
	default:
	{{p = ((te))-1;}}
	break;
	}
	}
	break;
#line 650 "wmkdepend.cpp"
		}
	}

_again:
	_acts = _scanInclude_actions + _scanInclude_to_state_actions[cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 3:
#line 1 "NONE"
	{ts = 0;}
	break;
	case 4:
#line 1 "NONE"
	{act = 0;}
	break;
#line 667 "wmkdepend.cpp"
		}
	}

	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	if ( _scanInclude_eof_trans[cs] > 0 ) {
		_trans = _scanInclude_eof_trans[cs] - 1;
		goto _eof_trans;
	}
	}

	_out: {}
	}

#line 366 "wmkdepend.rl"


        if (cs == scanInclude_error)
        {
            // Machine failed before finding a token
            std::cerr << "PARSE ERROR while scanning " << fileName << '\n';
            break;
        }

        // Now set up the prefix.
        if (ts == nullptr)
        {
            bytesPending = 0;
        }
        else
        {
            // There are data that needs to be shifted over.
            bytesPending = pe - ts;
            ::memmove(buf, ts, bytesPending);
            te -= (ts-buf);
            ts = buf;
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
