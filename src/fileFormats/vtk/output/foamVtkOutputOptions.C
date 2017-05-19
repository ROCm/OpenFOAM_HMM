/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "foamVtkOutputOptions.H"


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::foamVtkOutput::outputOptions&
Foam::foamVtkOutput::outputOptions::ascii(bool on)
{
    if (on)
    {
        switch (fmtType_)
        {
            case INLINE_BASE64: fmtType_ = INLINE_ASCII;   break;
            case APPEND_BINARY: fmtType_ = APPEND_BASE64;  break;
            case LEGACY_BINARY: fmtType_ = LEGACY_ASCII;   break;
            default: break; // no change
        }
    }
    else
    {
        switch (fmtType_)
        {
            case INLINE_ASCII:  fmtType_ = INLINE_BASE64;  break;
            case APPEND_BASE64: fmtType_ = APPEND_BINARY;  break;
            case LEGACY_ASCII:  fmtType_ = LEGACY_BINARY;  break;
            default: break; // no change
        }
    }

    return *this;
}


Foam::foamVtkOutput::outputOptions&
Foam::foamVtkOutput::outputOptions::append(bool on)
{
    if (on)
    {
        switch (fmtType_)
        {
            case INLINE_ASCII:  fmtType_ = APPEND_BASE64;  break;
            case LEGACY_ASCII:  fmtType_ = APPEND_BASE64;  break;
            case INLINE_BASE64: fmtType_ = APPEND_BINARY;  break;
            case LEGACY_BINARY: fmtType_ = APPEND_BINARY;  break;
            default: break; // no change
        }
    }
    else
    {
        switch (fmtType_)
        {
            case APPEND_BASE64: fmtType_ = INLINE_ASCII;   break;
            case APPEND_BINARY: fmtType_ = INLINE_BASE64;  break;
            default: break; // no change
        }
    }

    return *this;
}


Foam::foamVtkOutput::outputOptions&
Foam::foamVtkOutput::outputOptions::legacy(bool on)
{
    if (on)
    {
        switch (fmtType_)
        {
            case INLINE_ASCII:  fmtType_ = LEGACY_ASCII;   break;
            case INLINE_BASE64: fmtType_ = LEGACY_BINARY;  break;
            case APPEND_BASE64: fmtType_ = LEGACY_ASCII;   break;
            case APPEND_BINARY: fmtType_ = LEGACY_BINARY;  break;
            default: break; // no change
        }
    }
    else
    {
        switch (fmtType_)
        {
            case LEGACY_ASCII:  fmtType_ = INLINE_ASCII;   break;
            case LEGACY_BINARY: fmtType_ = INLINE_BASE64;  break;
            default: break; // no change
        }
    }

    return *this;
}


Foam::string Foam::foamVtkOutput::outputOptions::description() const
{
    switch (fmtType_)
    {
        case INLINE_ASCII:  return "xml ascii";
        case INLINE_BASE64: return "xml base64";
        case APPEND_BASE64: return "xml-append base64";
        case APPEND_BINARY: return "xml-append binary";
        case LEGACY_ASCII:  return "legacy ascii";
        case LEGACY_BINARY: return "legacy binary";
    }

    return "";
}


// ************************************************************************* //
