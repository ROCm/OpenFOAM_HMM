/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Class
    Foam::MRFZone

Description
    MRF zone definition based on both cell and face zones and parameters
    obtained from a control dictionary constructed from the given stream.

    The rotation of the MRF region is defined by an origin and axis of
    rotation and an angular speed.

SourceFiles
    MRFZone.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZone_H
#define MRFZone_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class MRFZone Declaration
\*---------------------------------------------------------------------------*/

class MRFZone
{
    // Private data

        const fvMesh& mesh_;

        word name_;

        dictionary dict_;

        label cellZoneID_;
        label faceZoneID_;
        wordList patchNames_;
        labelList patchLabels_;

        dimensionedVector origin_;
        dimensionedVector axis_;
        dimensionedScalar omega_;
        dimensionedVector Omega_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        MRFZone(const MRFZone&);

        //- Disallow default bitwise assignment
        void operator=(const MRFZone&);


public:

    // Constructors

        //- Construct from fvMesh and Istream
        MRFZone(const fvMesh& mesh, Istream& is);

        //- Return clone
        autoPtr<MRFZone> clone() const
        {
            notImplemented("autoPtr<MRFZone> clone() const");
            return autoPtr<MRFZone>(NULL);
        }

        //- Return a pointer to a new MRFZone created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<MRFZone> operator()(Istream& is) const
            {
                return autoPtr<MRFZone>(new MRFZone(mesh_, is));
            }
        };


    // Member Functions

        //- Add the Coriolis force contribution to the momentum equation
        void addCoriolis(fvVectorMatrix& UEqn) const;

        //- Make the given absolute flux relative within the MRF region
        void relativeFlux(surfaceScalarField& phi) const;

        //- Correct the boundary velocity for the roation of the MRF region
        void correctBoundaryVelocity(volVectorField& U) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const MRFZone&)
        {
            notImplemented("Ostream& operator<<(Ostream& os, const MRFZone&)");
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
