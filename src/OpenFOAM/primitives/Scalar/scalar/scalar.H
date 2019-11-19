/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

Typedef
    Foam::scalar

Description
    A floating-point number identical to float or double depending on
    whether WM_SP, WM_SPDP or WM_DP is defined.

SourceFiles
    scalar.C

\*---------------------------------------------------------------------------*/

#ifndef scalar_H
#define scalar_H

#include "floatScalar.H"
#include "doubleScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(WM_SP) || defined(WM_SPDP)

// Define scalar as a float

namespace Foam
{
    typedef floatScalar scalar;
    #if defined(WM_SPDP)
    typedef doubleScalar solveScalar;
    #else
    typedef floatScalar solveScalar;
    #endif

    constexpr scalar GREAT = floatScalarGREAT;
    constexpr scalar VGREAT = floatScalarVGREAT;
    constexpr scalar ROOTVGREAT = floatScalarROOTVGREAT;
    constexpr scalar SMALL = floatScalarSMALL;
    constexpr scalar ROOTSMALL = floatScalarROOTSMALL;
    constexpr scalar VSMALL = floatScalarVSMALL;
    constexpr scalar ROOTVSMALL = floatScalarROOTVSMALL;

    inline scalar readScalar(const char* buf)
    {
        return readFloat(buf);
    }

    inline bool readScalar(const char* buf, scalar& val)
    {
        return readFloat(buf, val);
    }

    inline scalar readScalar(const std::string& str)
    {
        return readFloat(str);
    }

    inline bool readScalar(const std::string& str, scalar& val)
    {
        return readFloat(str, val);
    }

    //- Read scalar from stream.
    scalar readScalar(Istream& is);

    //- Read raw scalar from binary stream.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    scalar readRawScalar(Istream& is);

    //- Read raw scalar(s) from binary stream.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    void readRawScalar(Istream& is, scalar* data, size_t nElem = 1);
}

#elif defined(WM_DP)

// Define scalar as a double

namespace Foam
{
    typedef doubleScalar scalar;
    typedef doubleScalar solveScalar;

    constexpr scalar GREAT = doubleScalarGREAT;
    constexpr scalar VGREAT = doubleScalarVGREAT;
    constexpr scalar ROOTVGREAT = doubleScalarROOTVGREAT;
    constexpr scalar SMALL = doubleScalarSMALL;
    constexpr scalar ROOTSMALL = doubleScalarROOTSMALL;
    constexpr scalar VSMALL = doubleScalarVSMALL;
    constexpr scalar ROOTVSMALL = doubleScalarROOTVSMALL;

    inline scalar readScalar(const char* buf)
    {
        return readDouble(buf);
    }

    inline bool readScalar(const char* buf, scalar& val)
    {
        return readDouble(buf, val);
    }

    inline scalar readScalar(const std::string& str)
    {
        return readDouble(str);
    }

    inline bool readScalar(const std::string& str, scalar& val)
    {
        return readDouble(str, val);
    }


    //- Read scalar from stream.
    scalar readScalar(Istream& is);

    //- Read raw scalar from binary stream.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    scalar readRawScalar(Istream& is);

    //- Read raw scalar(s) from binary stream.
    //  \note No internal check for binary vs ascii,
    //        the caller knows what they are doing
    void readRawScalar(Istream& is, scalar* data, size_t nElem = 1);
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Additional transcendental functions and specialisations

namespace Foam
{
    //- Inverse normalized incomplete gamma function
    scalar invIncGamma(const scalar a, const scalar P);

    //- Normalized upper incomplete gamma function
    scalar incGammaRatio_Q(const scalar a, const scalar x);

    //- Normalized lower incomplete gamma function
    scalar incGammaRatio_P(const scalar a, const scalar x);

    //- Upper incomplete gamma function
    scalar incGamma_Q(const scalar a, const scalar x);

    //- Lower incomplete gamma function
    scalar incGamma_P(const scalar a, const scalar x);

    //- Type to use for extended precision
    template<>
    class typeOfSolve<scalar>
    {
    public:

        typedef solveScalar type;
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
