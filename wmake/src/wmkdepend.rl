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

%%{
    machine wmkdep;
    write   data nofinal;

    hashInclude = space* '#' space* 'include' space* ; ## include directive
    includeName = [^\"<>]+ ;                    ## Text between "" or <>

    action  begInclude { ts_inclName = fpc; }   ## Remember start position

    ## No 'if (ts_inclName) ...' guard needed (implicit via the FSM)
    action  endInclude
    {
        // std::cerr << std::string(ts, (te-ts)) << '\n';
        processFile(std::string(ts_inclName, (fpc - ts_inclName)));
    }

    consume_comment := any* :>> '*/' @{ fgoto main; };

    main := |*

    # Single and double quoted strings
    ( 'L'? "'" ( [^'\\\n] | /\\./ )* "'") ;     # " swallow
    ( 'L'? '"' ( [^"\\\n] | /\\./ )* '"') ;     # ' swallow

    hashInclude '"' %begInclude includeName %endInclude '"' ;

    '/*' { fgoto consume_comment; };
    '//' [^\n]* '\n' ;
    [^\n]* '\n' ;              # Swallow all other lines

    *|;
}%%


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

    %% write init ;             ## /* FSM initialization here */;

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

        %% write exec ;         ## /* FSM execution here */;

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
