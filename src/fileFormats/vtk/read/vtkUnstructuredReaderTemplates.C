/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012 OpenFOAM Foundation
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

#include "vtkUnstructuredReader.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "stringIOList.H"
#include "cellModel.H"
#include "vectorIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::vtkUnstructuredReader::readBlock
(
    Istream& inFile,
    const label n,
    List<T>& list
) const
{
    list.setSize(n);
    for (T& val : list)
    {
        inFile >> val;
    }
}


template<class Type>
void Foam::vtkUnstructuredReader::printFieldStats
(
    const objectRegistry& obj
) const
{
    wordList fieldNames(obj.names(Type::typeName));

    if (fieldNames.size())
    {
        Info<< "Read " << fieldNames.size() << " " << Type::typeName
            << " fields:" << nl
            << "Size\tName" << nl
            << "----\t----" << endl;

        for (const word& fieldName : fieldNames)
        {
            Info<< obj.lookupObject<Type>(fieldName).size()
                << "\t" << fieldName
                << endl;
        }
        Info<< endl;
    }
}


// ************************************************************************* //
