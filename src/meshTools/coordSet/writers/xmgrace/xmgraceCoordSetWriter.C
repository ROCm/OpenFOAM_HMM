/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "xmgraceSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::xmgraceSetWriter<Type>::xmgraceSetWriter()
:
    writer<Type>()
{}


template<class Type>
Foam::xmgraceSetWriter<Type>::xmgraceSetWriter(const dictionary& dict)
:
    writer<Type>(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::xmgraceSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".agr";
}


template<class Type>
void Foam::xmgraceSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "@g0 on" << nl
        << "@with g0" << nl
        << "@    title \"" << points.name() << '"' << nl
        << "@    xaxis label " << '"' << points.axis() << '"' << nl;

    forAll(valueSets, i)
    {
        os  << "@    s" << i << " legend " << '"'
            << valueSetNames[i] << '"' << nl
            << "@target G0.S" << i << nl;

        this->writeTable(points, *valueSets[i], os);

        os  << '&' << nl;
    }
}


template<class Type>
void Foam::xmgraceSetWriter<Type>::write
(
    const bool writeTracks,
    const List<scalarField>& times,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorInFunction
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }
    if (tracks.size() > 0)
    {
        os  << "@g0 on" << nl
            << "@with g0" << nl
            << "@    title \"" << tracks[0].name() << '"' << nl
            << "@    xaxis label " << '"' << tracks[0].axis() << '"' << nl;

        // Data index.
        label sI = 0;

        forAll(tracks, trackI)
        {
            forAll(valueSets, i)
            {
                os  << "@    s" << sI << " legend " << '"'
                    << valueSetNames[i] << "_track" << i << '"' << nl
                    << "@target G0.S" << sI << nl;
                this->writeTable(tracks[trackI], valueSets[i][trackI], os);
                os  << '&' << nl;

                sI++;
            }
        }
    }
}

// ************************************************************************* //
