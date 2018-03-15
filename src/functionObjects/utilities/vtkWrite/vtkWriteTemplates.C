/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::label
Foam::functionObjects::vtkWrite::countFields() const
{
    return obr_.names<GeoField>(selectFields_).size();
}


template<class GeoField>
Foam::label
Foam::functionObjects::vtkWrite::writeFields
(
    vtk::internalWriter& writer,
    bool verbose
) const
{
    const wordList names = obr_.sortedNames<GeoField>(selectFields_);

    if (verbose && names.size())
    {
        Info<< "    " << GeoField::typeName
            << " " << flatOutput(names) << endl;
    }

    for (const word& fieldName : names)
    {
        writer.write(obr_.lookupObject<GeoField>(fieldName));
    }

    return names.size();
}


template<class GeoField>
Foam::label
Foam::functionObjects::vtkWrite::writeFields
(
    vtk::surfaceMeshWriter& writer,
    bool verbose
) const
{
    const wordList names = obr_.sortedNames<GeoField>(selectFields_);

    if (verbose && names.size())
    {
        Info<< "    " << GeoField::typeName
            << " " << flatOutput(names) << endl;
    }

    for (const word& fieldName : names)
    {
        writer.write(obr_.lookupObject<GeoField>(fieldName));
    }

    return names.size();
}


// ************************************************************************* //
