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

#include "foamVtkAppendBase64Formatter.H"
#include "foamVtkAppendRawFormatter.H"
#include "foamVtkAsciiFormatter.H"
#include "foamVtkBase64Formatter.H"
#include "foamVtkLegacyFormatter.H"

#include "IOstream.H"


// * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * //

Foam::foamVtkOutputOptions::foamVtkOutputOptions()
:
    type_(ASCII),
    precision_(IOstream::defaultPrecision())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::foamVtkFormatter>
Foam::foamVtkOutputOptions::newFormatter(std::ostream& os) const
{
    switch (type_)
    {
        case (LEGACY | BINARY):
            return autoPtr<foamVtkFormatter>
            (
                new foamVtkLegacyFormatter(os)
            );

        case BASE64:  // xml insitu
            return autoPtr<foamVtkFormatter>
            (
                new foamVtkBase64Formatter(os)
            );

        case (APPEND | BASE64):
            return autoPtr<foamVtkFormatter>
            (
                new foamVtkAppendBase64Formatter(os)
            );

        case (APPEND | BINARY):
            return autoPtr<foamVtkFormatter>
            (
                new foamVtkAppendRawFormatter(os)
            );

        default:   // ASCII (legacy or xml) must always work
            return autoPtr<foamVtkFormatter>
            (
                new foamVtkAsciiFormatter(os, precision_)
            );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::foamVtkOutputOptions::ascii(bool on)
{
    if (on)
    {
        // Force ASCII:

        if (type_ & APPEND)
        {
            // Append: ascii = base64 (vs raw binary)
            type_ = (APPEND | BASE64);
        }
        else if (type_ & LEGACY)
        {
            // Legacy: ascii = ascii
            type_ = (LEGACY | ASCII);
        }
        else
        {
            // XML: ascii = ascii
            type_ = ASCII;
        }
    }
    else
    {
        // Non-ASCII:

        if (type_ & APPEND)
        {
            // Append: binary == (raw) binary
            type_ = APPEND | BINARY;
        }
        else if (type_ & LEGACY)
        {
            // Legacy: binary = binary
            type_ = LEGACY | BINARY;
        }
        else
        {
            // XML: binary == (inline) binary == base64
            type_ = BASE64;
        }
    }
}


void Foam::foamVtkOutputOptions::append(bool on)
{
    if (on)
    {
        if (!(type_ & APPEND))
        {
            // XML:    base64 -> raw binary, ascii -> base64
            // Legacy: binary -> raw binary, ascii -> base64
            type_ = APPEND | ((type_ & (BASE64 | BINARY)) ? BINARY : BASE64);
        }
    }
    else if (type_ & APPEND)
    {
        // Only revert back to inline XML base64 versions
        // ASCII needs another step.

        type_ = BASE64;
    }
}


void Foam::foamVtkOutputOptions::legacy(bool on)
{
    if (on)
    {
        if (type_ & APPEND)
        {
            // Append: base64 -> ascii, binary -> binary
            type_ = (LEGACY | ((type_ & BINARY) ? BINARY : ASCII));
        }
        else if (type_ & LEGACY)
        {
            // no-op
        }
        else
        {
            // XML: ascii -> ascii, base64 -> binary
            type_ = (LEGACY | ((type_ & BASE64) ? BINARY : ASCII));
        }
    }
    else if (type_ & LEGACY)
    {
        // Legacy: ascii -> xml ascii, binary -> xml base64
        type_ = (type_ & BINARY) ? BASE64 : ASCII;
    }
}


void Foam::foamVtkOutputOptions::precision(unsigned prec) const
{
    precision_ = prec;
}


Foam::Ostream& Foam::foamVtkOutputOptions::info(Ostream& os) const
{
    os << "type: " << type_;

    switch (type_)
    {
        case (LEGACY | ASCII):
            os << " legacy ascii";
            break;

        case (LEGACY | BINARY):
            os << " legacy binary";
            break;

        case BASE64:
            os << " xml insitu base64";
            break;

        case (APPEND | BASE64):
            os << " xml-append base64";
            break;

        case (APPEND | BINARY):
            os << " xml-append binary";
            break;

        case ASCII:
            os << " xml insitu ascii";
            break;

        case BINARY:
            os << " xml insitu binary - WRONG";
            break;

        default:
            os << " unknown";
            break;
    }

    if (legacy()) os << " legacy";
    if (xml())    os << " xml";
    if (append()) os << " append";
    if (ascii())  os << " ascii";

    return os;
}


// ************************************************************************* //
