/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

Class
    Foam::cloudInjection

Description
    Cloud injection model

SourceFiles
    cloudInjection.C

\*---------------------------------------------------------------------------*/

#ifndef cloudInjection_H
#define cloudInjection_H

#include "injectionModel.H"
#include "kinematicCloud.H"
#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class cloudInjection Declaration
\*---------------------------------------------------------------------------*/

class cloudInjection
:
    public injectionModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        cloudInjection(const cloudInjection&);

        //- Disallow default bitwise assignment
        void operator=(const cloudInjection&);


protected:

    // Protected data

        //- Name of the cloud owner for newly ejected parcels
        word cloudName_;

        //- Reference to the cloud
        kinematicCloud& cloud_;

        //- Number of particles per parcel
        scalar particlesPerParcel_;

        //- Random number generator
        Random rndGen_;

        //- Parcel size PDF model
        const autoPtr<pdfs::pdf> parcelPDF_;

        //- Diameters of particles to inject into the cloud
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("cloudInjection");


    // Constructors

        //- Construct from surface film model
        cloudInjection(const surfaceFilmModel& owner, const dictionary& dict);


    //- Destructor
    virtual ~cloudInjection();


    // Member Functions

        // Access

            //- Return the cloud name
            inline const word& cloudName() const;

            //- Return a reference to the cloud
            inline const kinematicCloud& cloud() const;


        // Evolution

            //- Correct
            virtual void correct
            (
                scalarField& massToInject,
                scalarField& diameterToInject
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cloudInjectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
