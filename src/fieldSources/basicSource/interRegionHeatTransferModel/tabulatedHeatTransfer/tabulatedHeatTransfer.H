/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::tabulatedHeatTransfer

Description
    Tabulated heat transfer model. The heat exchange area must be
    provided. The 2D table look for heat transfer coefficients uses the
    primary and secondary region velocities

\*---------------------------------------------------------------------------*/

#ifndef tabulatedHeatTransfer_H
#define tabulatedHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class tabulatedHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedHeatTransfer
:
    public interRegionHeatTransferModel
{
private:

    // Private data

        //- 2D look up table
        autoPtr<interpolation2DTable<scalar> > hTable_;

        //- Area of heat exchange
        autoPtr<volScalarField> area_;


public:

    //- Runtime type information
        TypeName("tabulatedHeatTransfer");


    // Constructors

        //- Construct from dictionary
        tabulatedHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~tabulatedHeatTransfer();


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual const tmp<volScalarField> calculateHtc();


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict) ;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
