/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "nastranSetWriter.H"
#include "coordSet.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::nastranSetWriter<Type>::nastranSetWriter()
:
    writer<Type>()
{}


template<class Type>
Foam::nastranSetWriter<Type>::nastranSetWriter(const dictionary& dict)
:
    writer<Type>(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::nastranSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".nas";
}


template<class Type>
void Foam::nastranSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "TITLE=OpenFOAM "
        << this->getBaseName(points, valueSetNames).c_str()
        << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

    forAll(points, pointi)
    {
        fileFormats::NASCore::writeCoord
        (
            os, points[pointi], pointi, fieldFormat::FREE
        );
    }

    if (false)
    {
        // Single track with multiple segments
        const label nEdges = points.size()-1;
        for (label edgei = 0; edgei < nEdges; ++edgei)
        {
            fileFormats::NASCore::writeKeyword
            (
                os,
                "PLOTEL",
                fieldFormat::FREE
            );

            // fieldFormat::SHORT
            //os.setf(std::ios_base::right);
            //os  << setw(8) << edgei+1
            //    << setw(8) << edgei+1
            //    << setw(8) << edgei+2
            //    << nl;
            //os.unsetf(std::ios_base::right);

            // fieldFormat::FREE
            os  << ',' << edgei+1
                << ',' << edgei+1
                << ',' << edgei+2
                << nl;
        }
    }

    os << "ENDDATA" << nl;
}


template<class Type>
void Foam::nastranSetWriter<Type>::write
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
    if (tracks.empty())
    {
        return;
    }

    os  << "TITLE=OpenFOAM "
        << this->getBaseName(tracks[0], valueSetNames).c_str()
        << nl
        << "$" << nl
        << "BEGIN BULK" << nl;

//    label nTracks = tracks.size();
//    label nPoints = 0;
//    forAll(tracks, i)
//    {
//        nPoints += tracks[i].size();
//    }

    label globalPointi = 0;
    for (const coordSet& points : tracks)
    {
        for (const point& p : points)
        {
            fileFormats::NASCore::writeCoord
            (
                os, p, globalPointi, fieldFormat::FREE
            );
            ++globalPointi;
        }
    }

    if (writeTracks)
    {
        // Write ids of track points to file
        label globalEdgei = 0;
        label globalPointi = 0;
        for (const coordSet& points : tracks)
        {
            const label nEdges = points.size()-1;
            for (label edgei = 0; edgei < nEdges; ++edgei)
            {
                fileFormats::NASCore::writeKeyword
                (
                    os,
                    "PLOTEL",
                    fieldFormat::FREE
                );

                // fieldFormat::SHORT
                //os.setf(std::ios_base::right);
                //os  << setw(8) << globalEdgei+1
                //    << setw(8) << globalPointi+1
                //    << setw(8) << globalPointi+2
                //    << nl;
                //os.unsetf(std::ios_base::right);

                // fieldFormat::FREE
                os  << ',' << globalEdgei+1
                    << ',' << globalPointi+1
                    << ',' << globalPointi+2
                    << nl;

                ++globalEdgei;
                ++globalPointi;
            }
        }
    }

    os << "ENDDATA" << nl;
}


// ************************************************************************* //
