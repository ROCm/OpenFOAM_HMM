/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "cpuInfo.H"
#include "IFstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope function
template<class T>
inline static void writeEntry
(
    Foam::Ostream& os, const Foam::word& key, const T& value
)
{
    os.writeKeyword(key) << value << Foam::token::END_STATEMENT << '\n';
}


// file-scope function
static bool split(std::string& line, std::string& key, std::string& val)
{
    std::string::size_type sep = line.find(':');

    if (sep == std::string::npos)
    {
        return false;
    }

    std::string::size_type endKey = line.find_last_not_of("\t:", sep);
    std::string::size_type begVal = line.find_first_not_of(" :", sep);

    if (endKey == std::string::npos || begVal == std::string::npos)
    {
        return false;
    }
    ++endKey;

    // replace spaces in key with '_' for ease of use/consistency
    for
    (
        std::string::iterator iter = line.begin();
        iter != line.end();
        ++iter
    )
    {
        if (*iter == ' ')
        {
            *iter = '_';
        }
        else if (*iter == ':')
        {
            break;
        }
    }

    key = line.substr(0, endKey);
    val = line.substr(begVal);

    // std::cerr<<"key=" << key << " val= " << val << '\n';

    return true;
}


// file-scope function - get int
static inline bool getInt(const std::string& str, int& val)
{
    int i;
    if (sscanf(str.c_str(), "%d", &i) == 1)
    {
        val = i;
        return true;
    }
    else
    {
        return false;
    }
}

// file-scope function - get float
static inline bool getFlt(const std::string& str, float& val)
{
    float f;
    if (sscanf(str.c_str(), "%f", &f) == 1)
    {
        val = f;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// parse this type of content:
// ===========================
// processor       : 0
// vendor_id       : GenuineIntel
// cpu family      : 6
// model           : 63
// model name      : Intel(R) Xeon(R) CPU E5-2620 v3 @ 2.40GHz
// stepping        : 2
// microcode       : 0x35
// cpu MHz         : 1200.000
// cache size      : 15360 KB
// physical id     : 0
// siblings        : 12
// core id         : 0
// cpu cores       : 6
// apicid          : 0
// initial apicid  : 0
// fpu             : yes
// fpu_exception   : yes
// cpuid level     : 15
// wp              : yes
// flags           : fpu vme ...
// bugs            :
// bogomips        : 4789.15
// clflush size    : 64
// cache_alignment : 64
// address sizes   : 46 bits physical, 48 bits virtual
// power management:

void Foam::cpuInfo::parse()
{
    int ncpu = 0;

    IFstream is("/proc/cpuinfo");
    while (is.good())
    {
        string line, key, value;
        is.getLine(line);

        if (!split(line, key, value))
        {
            continue;
        }

        if (key == "processor")
        {
            if (ncpu++)
            {
                break; // stop after the first cpu
            }
        }
        else if (key == "vendor_id")   { vendor_id  = value;        }
        else if (key == "model_name")  { model_name = value;        }
        else if (key == "cpu_family")  { getInt(value, cpu_family); }
        else if (key == "model")       { getInt(value, model);      }
        else if (key == "cpu_MHz")     { getFlt(value, cpu_MHz);    }
        else if (key == "cpu_cores")   { getInt(value, cpu_cores);  }
        else if (key == "siblings")    { getInt(value, siblings);   }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuInfo::cpuInfo()
:
    vendor_id(),
    model_name(),
    cpu_family(-1),
    model(-1),
    cpu_MHz(0),
    siblings(0),
    cpu_cores(0)
{
    parse();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cpuInfo::~cpuInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cpuInfo::write(Ostream& os) const
{
    if (!vendor_id.empty())
    {
        writeEntry(os, "vendor_id",     vendor_id);
    }
    if (!model_name.empty())
    {
        writeEntry(os, "model_name",    model_name);
    }
    if (cpu_family != -1)
    {
        writeEntry(os, "cpu_family",    cpu_family);
    }
    if (model != -1)
    {
        writeEntry(os, "model",         model);
    }
    if (cpu_MHz > 0)
    {
        writeEntry(os, "cpu_MHz",       cpu_MHz);
    }
    if (cpu_cores > 0)
    {
        writeEntry(os, "cpu_cores",     cpu_cores);
    }
    if (siblings > 0)
    {
        writeEntry(os, "siblings",      siblings);
    }
}


// ************************************************************************* //
