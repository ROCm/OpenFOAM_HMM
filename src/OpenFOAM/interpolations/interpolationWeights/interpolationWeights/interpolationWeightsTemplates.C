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

#include "interpolationWeights.H"
#include "ListOps.H"
#include "IOobject.H"
#include "HashSet.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class GeoField>
//void interpolationWeights::readFields
//(
//    const word& fieldName,
//    const typename GeoField::Mesh& mesh,
//    const wordList& timeNames,
//    objectRegistry& fieldsCache
//)
//{
//    // Collect all times that are no longer used
//    {
//        HashSet<word> usedTimes(timeNames);
//
//        DynamicList<word> unusedTimes(fieldsCache.size());
//
//        forAllIter(objectRegistry, fieldsCache, timeIter)
//        {
//            const word& tm = timeIter.key();
//            if (!usedTimes.found(tm))
//            {
//                unusedTimes.append(tm);
//            }
//        }
//
//        //Info<< "Unloading times " << unusedTimes << endl;
//
//        forAll(unusedTimes, i)
//        {
//            objectRegistry& timeCache = const_cast<objectRegistry&>
//            (
//                fieldsCache.lookupObject<objectRegistry>(unusedTimes[i])
//            );
//            fieldsCache.checkOut(timeCache);
//        }
//    }
//
//
//    // Load any new fields
//    forAll(timeNames, i)
//    {
//        const word& tm = timeNames[i];
//
//        // Create if not found
//        if (!fieldsCache.found(tm))
//        {
//            //Info<< "Creating registry for time " << tm << endl;
//
//            // Create objectRegistry if not found
//            objectRegistry* timeCachePtr = new objectRegistry
//            (
//                IOobject
//                (
//                    tm,
//                    tm,
//                    fieldsCache,
//                    IOobject::NO_READ,
//                    IOobject::NO_WRITE
//                )
//            );
//            timeCachePtr->store();
//        }
//
//        // Obtain cache for current time
//        const objectRegistry& timeCache =
//            fieldsCache.lookupObject<objectRegistry>
//            (
//                tm
//            );
//
//        // Store field if not found
//        if (!timeCache.found(fieldName))
//        {
//            //Info<< "Loading field " << fieldName
//            //    << " for time " << tm << endl;
//
//            GeoField loadedFld
//            (
//                IOobject
//                (
//                    fieldName,
//                    tm,
//                    mesh.thisDb(),
//                    IOobject::MUST_READ,
//                    IOobject::NO_WRITE,
//                    false
//                ),
//                mesh
//            );
//
//            // Transfer to timeCache (new objectRegistry and store flag)
//            GeoField* fldPtr = new GeoField
//            (
//                IOobject
//                (
//                    fieldName,
//                    tm,
//                    timeCache,
//                    IOobject::NO_READ,
//                    IOobject::NO_WRITE
//                ),
//                loadedFld
//            );
//            fldPtr->store();
//        }
//    }
//}
//
//
//template<class GeoField>
//void interpolationWeights::readFields
//(
//    const word& fieldName,
//    const typename GeoField::Mesh& mesh,
//    const wordList& timeNames,
//    const word& registryName
//)
//{
//    readFields<GeoField>
//    (
//        fieldName,
//        mesh,
//        timeNames,
//        //registry(mesh.thisDb(), registryName)
//        const_cast<objectRegistry&>
//        (
//            mesh.thisDb().subRegistry(registryName, true)
//        )
//    );
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
