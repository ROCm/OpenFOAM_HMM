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
    Foam::interRegionHeatTransferModel

Description
    Base class for inter region heat exchange. The derived classes must
    provide the heat transfer coefficient (htc)
    NOTE: mapToMap does to work in paralell

\*---------------------------------------------------------------------------*/

#ifndef interRegionHeatTransferModel_H
#define interRegionHeatTransferModel_H

#include "basicSource.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interRegionHeatTransferModel Declaration
\*---------------------------------------------------------------------------*/

class interRegionHeatTransferModel
:
    public  basicSource
{

private:

    // Private data

        //- Pointer to secondary interRegionHeatTransferModel
        interRegionHeatTransferModel* secIrht_;

        //- First iteration
        bool firstIter_;


    // Private members

        //- Check coupled interRegionHeatTransferModel
        void check();


protected:

    // Protected data

        //- Heat transfer coefficient [W/m2/k] by area/volume [1/m]
        // registered on the master mesh
        volScalarField htc_;


public:

    //- Runtime type information
    TypeName("interRegionHeatTransferModel");


    // Constructors

        //- Construct from dictionary
        interRegionHeatTransferModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~interRegionHeatTransferModel();


    // Member Functions


        // Access

            //- Return the heat transfer coefficient
            const volScalarField& htc() const
            {
                return  htc_;
            }

            //-Source term to fvMatrix<scalar>
            virtual void addSup(fvMatrix<scalar>& eEqn, const label fieldI);

            //- Calculate heat transfer coefficient
            virtual const tmp<volScalarField> calculateHtc() = 0;


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
