/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "StringStream.H"

#include <cinttypes>
#include <cxxabi.h>
#include <execinfo.h>
#include <dlfcn.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string pOpen(const string& cmd, label line=0)
{
    string res;

    FILE *cmdPipe = popen(cmd.c_str(), "r");
    if (cmdPipe)
    {
        char *buf = nullptr;

        // Read line number of lines
        for (label cnt = 0; cnt <= line; ++cnt)
        {
            size_t linecap = 0;
            ssize_t linelen = ::getline(&buf, &linecap, cmdPipe);

            if (linelen < 0)
            {
                break;
            }

            if (cnt == line)
            {
                res = string(buf);
                // Trim trailing newline
                if (res.size())
                {
                    res.resize(res.size()-1);
                }
                break;
            }
        }

        if (buf != nullptr)
        {
            free(buf);
        }

        pclose(cmdPipe);
    }

    return res;
}


inline word addressToWord(const uintptr_t addr)
{
    OStringStream os;
    #ifdef __APPLE__
    os << "0x" << hex << uint64_t(addr);
    #else
    os << "0x" << hex << addr;
    #endif
    return os.str();
}


inline string& shorterPath(string& s)
{
    s.replace(cwd() + '/', "");
    s.replace(home(), "~");
    return s;
}


void printSourceFileAndLine
(
    Ostream& os,
    const fileName& filename,
    Dl_info *info,
    void *addr
)
{
    uintptr_t address = uintptr_t(addr);
    word myAddress = addressToWord(address);

    // Can use relative addresses for executables and libraries with the
    // Darwin addr2line implementation.
    // On other systems (Linux), only use relative addresses for libraries.

    #ifndef __APPLE__
    if (filename.hasExt("so"))
    #endif
    {
        // Convert address into offset into dynamic library
        uintptr_t offset = uintptr_t(info->dli_fbase);
        intptr_t relativeAddress = address - offset;
        myAddress = addressToWord(relativeAddress);
    }

    if (filename[0] == '/')
    {
        string line = pOpen
        (
            "addr2line -f --demangle=auto --exe "
          + filename
          + " "
          + myAddress,
            1
        );

        if (line.empty())
        {
            os  << " addr2line failed";
        }
        else if (line == "??:0")
        {
            line = filename;
            os  << " in " << shorterPath(line).c_str();
        }
        else
        {
            os  << " at " << shorterPath(line).c_str();
        }
    }
}


fileName absolutePath(const char* fn)
{
    fileName fname(fn);

    if (fname[0] != '/' && fname[0] != '~')
    {
        string tmp = pOpen("which " + fname);

        if (tmp[0] == '/' || tmp[0] == '~')
        {
            fname = tmp;
        }
    }

    return fname;
}


word demangleSymbol(const char* sn)
{
    int st;
    char* cxx_sname = abi::__cxa_demangle
    (
        sn,
        nullptr,
        0,
        &st
    );

    if (st == 0 && cxx_sname)
    {
        word demangled(cxx_sname);
        free(cxx_sname);

        return demangled;
    }

    return sn;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::error::safePrintStack(std::ostream& os)
{
    // Get raw stack symbols
    void *array[100];
    size_t size = backtrace(array, 100);
    char **strings = backtrace_symbols(array, size);

    // See if they contain function between () e.g. "(__libc_start_main+0xd0)"
    // and see if cplus_demangle can make sense of part before +
    for (size_t i = 0; i < size; ++i)
    {
        string msg(strings[i]);
        fileName programFile;
        word address;

        os  << '#' << label(i) << '\t' << msg << std::endl;
    }
}


void Foam::error::printStack(Ostream& os)
{
    // Get raw stack symbols
    const size_t CALLSTACK_SIZE = 128;

    void *callstack[CALLSTACK_SIZE];
    size_t size = backtrace(callstack, CALLSTACK_SIZE);

    Dl_info *info = new Dl_info;

    fileName fname = "???";
    word address;

    for (size_t i=0; i<size; ++i)
    {
        int st = dladdr(callstack[i], info);

        os << '#' << label(i) << "  ";
        if (st != 0 && info->dli_fname != nullptr && info->dli_fname[0] != '\0')
        {
            fname = absolutePath(info->dli_fname);

            os <<
            (
                (info->dli_sname != nullptr)
              ? demangleSymbol(info->dli_sname)
              : "?"
            );
        }
        else
        {
            os << "?";
        }

        printSourceFileAndLine(os, fname, info, callstack[i]);
        os << nl;
    }

    delete info;
}


// ************************************************************************* //
