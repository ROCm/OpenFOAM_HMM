/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fieldValues::faceSourceDelta::processFields(bool& found)
{
    typedef GeometricField<Type, fvPatchField, volMesh> vf;

    const wordList& fields1 = faceSource1_.fields();

    const dictionary& results1 = faceSource1_.resultDict();
    const dictionary& results2 = faceSource2_.resultDict();

    Type r1(pTraits<Type>::zero);
    Type r2(pTraits<Type>::zero);

    forAll(fields1, i)
    {
        const word& fieldName = fields1[i];
        if (obr_.foundObject<vf>(fieldName) && results2.found(fieldName))
        {
            results1.lookup(fieldName) >> r1;
            results2.lookup(fieldName) >> r2;

            if (log_)
            {
                Info<< "    field: " << fieldName << ", delta: " << r2 - r1
                    << endl;
            }

            if (Pstream::master())
            {
                file()<< tab << r2 - r1;
            }

            found = true;
        }
    }
}


// ************************************************************************* //
