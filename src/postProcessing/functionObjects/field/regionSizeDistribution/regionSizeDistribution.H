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

\addtogroup functionObjects
@{

Class
    Foam::regionSizeDistribution

Description
    Looks up a field, interpolates it to the faces and determines a connected
    region from a patch where the field is above a certain value.
    - Writes a field containing all regions starting at given patch
      ('liquid core')
    - All other regions are summed for volume and a histogram is calculated.

SourceFiles
    regionSizeDistribution.C
@}

\*---------------------------------------------------------------------------*/

#ifndef regionSizeDistribution_H
#define regionSizeDistribution_H

#include "pointFieldFwd.H"
#include "writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class regionSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class regionSizeDistribution
{
    // Private data

        //- Name of this set of regionSizeDistribution objects
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Name of field
        word alphaName_;

        //- Patches to walk from
        wordList patchNames_;

        //- Clip value
        scalar threshold_;

        //- Background region volFraction
        scalar volFraction_;

        //- Mumber of bins
        label nBins_;

        //- Output formatter to write
        autoPtr<writer<scalar> > formatterPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionSizeDistribution(const regionSizeDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const regionSizeDistribution&);


public:

    //- Runtime type information
    TypeName("regionSizeDistribution");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        regionSizeDistribution
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~regionSizeDistribution();


    // Member Functions

        //- Return name of the set of regionSizeDistribution
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the regionSizeDistribution data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the regionSizeDistribution and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
