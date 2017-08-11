/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "memInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"

#include <fstream>
#include <string>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::memInfo::memInfo()
:
    peak_(0),
    size_(0),
    rss_(0)
{
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::memInfo::~memInfo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
//
// Parse the following type of content.
//
// ===========================
// VmPeak:    15920 kB
// VmSize:    15916 kB
// VmLck:         0 kB
// VmPin:         0 kB
// VmHWM:      6972 kB
// VmRSS:      6972 kB
// VmLib:      2208 kB
// VmPTE:        52 kB
// VmPMD:        12 kB
// VmSwap:        0 kB

const Foam::memInfo& Foam::memInfo::update()
{
    // Clear (invalidate) values first
    peak_ = size_ = rss_ = 0;
    std::string line;

    unsigned nKeys = 0;

    std::ifstream is("/proc/" + std::to_string(Foam::pid()) + "/status");
    while (is.good() && nKeys < 3)  // Stop after getting the known keys
    {
        std::getline(is, line);

        const auto keyLen = line.find(':');
        if (keyLen == std::string::npos)
        {
            continue;
        }

        // Value is after the ':', but skip any leading whitespace since
        // strtoi will do it anyhow
        const auto begVal = line.find_first_not_of("\t :", keyLen);
        if (begVal == std::string::npos)
        {
            continue;
        }

        const std::string key = line.substr(0, keyLen);

        if (key == "VmPeak")
        {
            peak_ = std::stoi(line.substr(begVal));
            ++nKeys;
        }
        else if (key == "VmSize")
        {
            size_ = std::stoi(line.substr(begVal));
            ++nKeys;
        }
        else if (key == "VmRSS")
        {
            rss_ = std::stoi(line.substr(begVal));
            ++nKeys;
        }
    }

    return *this;
}


bool Foam::memInfo::valid() const
{
    return peak_ > 0;
}


void Foam::memInfo::write(Ostream& os) const
{
    os.writeEntry("size", size_);
    os.writeEntry("peak", peak_);
    os.writeEntry("rss", rss_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, memInfo& m)
{
    is.readBegin("memInfo");
    is  >> m.peak_ >> m.size_ >> m.rss_;
    is.readEnd("memInfo");

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const memInfo& m)
{
    os  << token::BEGIN_LIST
        << m.peak_ << token::SPACE
        << m.size_ << token::SPACE
        << m.rss_
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
