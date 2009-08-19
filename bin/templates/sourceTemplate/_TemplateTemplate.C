/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "CLASSNAME.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<TemplateClassArgument>
const dataType Foam::CLASSNAME<TemplateArgument>::staticData();


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::CLASSNAME<TemplateArgument>::CLASSNAME()
:
    baseClassName(),
    data_()
{}


template<TemplateClassArgument>
Foam::CLASSNAME<TemplateArgument>::CLASSNAME(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


template<TemplateClassArgument>
Foam::CLASSNAME<TemplateArgument>::CLASSNAME
(
    const CLASSNAME<TemplateArgument>&
)
:
    baseCLASSNAME(),
    data_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::autoPtr<Foam::CLASSNAME<TemplateArgument> >
Foam::CLASSNAME<TemplateArgument>::New()
{
    return autoPtr<CLASSNAME<TemplateArgument> >
    (
        new CLASSNAME<TemplateArgument>
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::CLASSNAME<TemplateArgument>::~CLASSNAME()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<TemplateClassArgument>
void Foam::CLASSNAME<TemplateArgument>::operator=
(
    const CLASSNAME<TemplateArgument>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::CLASSNAME<TemplateArgument>::operator="
            "(const Foam::CLASSNAME<TemplateArgument>&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
