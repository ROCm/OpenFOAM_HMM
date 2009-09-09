/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::thermalPorousZone

Description
    Porous zone definition based on cell zones including terms for energy
    equations.

See Also
    porousZone, thermalPorousZones and coordinateSystems

SourceFiles
    thermalPorousZone.C
    thermalPorousZoneTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef thermalPorousZone_H
#define thermalPorousZone_H

#include "porousZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class basicThermo;

/*---------------------------------------------------------------------------*\
                        Class thermalPorousZone Declaration
\*---------------------------------------------------------------------------*/

class thermalPorousZone
:
    public porousZone
{
    // Private data

        //- Temperature in the porous-zone
        dimensionedScalar T_;


        //- Disallow default bitwise copy construct
        thermalPorousZone(const thermalPorousZone&);

        //- Disallow default bitwise assignment
        void operator=(const thermalPorousZone&);


public:

    // Constructors

        //- Construct from components
        thermalPorousZone(const word& name, const fvMesh&, const dictionary&);

        //- Return clone
        autoPtr<thermalPorousZone> clone() const
        {
            notImplemented("autoPtr<thermalPorousZone> clone() const");
            return autoPtr<thermalPorousZone>(NULL);
        }

        //- Return pointer to new thermalPorousZone
        //  created on freestore from Istream
        class iNew
        {
            //- Reference to the finite volume mesh this zone is part of
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<thermalPorousZone> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);

                return autoPtr<thermalPorousZone>
                (
                    new thermalPorousZone(name, mesh_, dict)
                );
            }
        };


    //- Destructor
    virtual ~thermalPorousZone()
    {}


    // Member Functions

        // Access

            //- Return the temperature in the porous-zone
            const dimensionedScalar& T() const
            {
                return T_;
            }

            //- Edit access to the temperature in the porous-zone
            dimensionedScalar& T()
            {
                return T_;
            }

        //- Add the thermal source to the enthalpy equation
        void addEnthalpySource
        (
            const basicThermo& thermo,
            fvScalarMatrix& hEqn
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
//#   include "thermalPorousZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
