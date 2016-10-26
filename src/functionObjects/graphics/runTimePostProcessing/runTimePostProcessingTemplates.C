/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
void Foam::functionObjects::runTimePostProcessing::readObjects
(
    const dictionary& dict,
    PtrList<Type>& objects
) const
{
    objects.clear();
    forAllConstIter(dictionary, dict, iter)
    {
        if (!iter().isDict())
        {
            FatalIOErrorInFunction(dict)
                << dict.dictName()
                << " objects must be specified in dictionary format"
                << exit(FatalIOError);
        }

        const dictionary& objectDict(iter().dict());
        word objectType = objectDict.lookup("type");

        objects.append
        (
            Type::New(*this, iter().dict(), scene_.colours(), objectType)
        );
    }
}


// ************************************************************************* //
