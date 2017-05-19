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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::foamVtkOutput::outputOptions::outputOptions()
:
    fmtType_(INLINE_ASCII),
    precision_(IOstream::defaultPrecision())
{}


inline Foam::foamVtkOutput::outputOptions::outputOptions
(
    enum formatType fmtType
)
:
    fmtType_(fmtType),
    precision_(IOstream::defaultPrecision())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::foamVtkOutput::formatter>
Foam::foamVtkOutput::outputOptions::newFormatter(std::ostream& os) const
{
    return foamVtkOutput::newFormatter(os, fmtType_, precision_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::foamVtkOutput::formatType
Foam::foamVtkOutput::outputOptions::format() const
{
    return fmtType_;
}


inline bool Foam::foamVtkOutput::outputOptions::legacy() const
{
    return (fmtType_ & formatType::LEGACY_ASCII & formatType::LEGACY_BINARY);
}


inline bool Foam::foamVtkOutput::outputOptions::xml() const
{
    return !legacy();
}


inline bool Foam::foamVtkOutput::outputOptions::append() const
{
    return (fmtType_ & formatType::APPEND_BASE64 & formatType::APPEND_BINARY);
}


inline bool Foam::foamVtkOutput::outputOptions::insitu() const
{
    return !append();
}


inline bool Foam::foamVtkOutput::outputOptions::ascii() const
{
    return !(fmtType_ & 0xF);
}


inline const Foam::foamVtkOutput::outputOptions&
Foam::foamVtkOutput::outputOptions::precision(unsigned prec) const
{
    precision_ = prec;
    return *this;
}


// ************************************************************************* //
