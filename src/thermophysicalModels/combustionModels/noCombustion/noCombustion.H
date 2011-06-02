/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    Foam::combustionModel::noCombustion

Description
    Dummy combustion model for 'none' option

SourceFiles
    noCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef noCombustion_H
#define noCombustion_H

#include "combustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                      Class noCombustion Declaration
\*---------------------------------------------------------------------------*/

class noCombustion
:
    public combustionModel
{
    // Private Member Functions

        //- Disallow copy construct
        noCombustion(const noCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const noCombustion&);


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noCombustion
        (
            const dictionary& combustionProperties,
            hsCombustionThermo& thermo,
            const compressible::turbulenceModel& turbulence,
            const surfaceScalarField& phi,
            const volScalarField& rho
        );


    //- Destructor
    virtual ~noCombustion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
