/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::fv::effectivenessHeatExchangerSource

Description
    Heat exchanger source model.
    The total heat exchange is given by:

    Qt = e(phi, secondaryMassFlowRate)
        * (secondaryInletT - primaryInletT)
        * phi*Cp;

    where:
        e(phi, secondaryMassFlowRate) is a effectivenes table in function of
        mass net flux coming imto the heat exchanger (phi) and
        secondaryMassFlowRate

    Then the distribution inside the hear exchanger is given by:

    Qcell = Vcell*|Ucell|*(Tcell - Tref)/sum(Vcell*|Ucell|*(Tcell - Tref));

    where:
        Qcell is the source for cell
        Vcell is the volume of the cell
        Ucell is the local cell velocity
        Tcell is the local call temperature
        Tref : min(T) or max(T) in the cell zone depending on the sign of Qt

    Example :

    effectivenessHeatExchangerSource1
    {
        type            effectivenessHeatExchangerSource;
        active          true;
        selectionMode   cellZone;
        cellZone        porosity;

        effectivenessHeatExchangerSourceCoeffs
        {
            fieldNames              (e);
            secondaryMassFlowRate   1.0;
            secondaryInletT         336;
            primaryInletT           293;
            faceZone                facesZoneInletOriented;
            outOfBounds             clamp;
            fileName                "effTable";
        }
    }

    Note:
        the table with name "fileName" should have the same units as
        secondaryMassFlowRate and kg/s for phi
        faceZone is the faces at the inlet of the cellzone, it needs to be
        created with flip map flags. It is used to integrate the net mass flow
        rate into the heat exchanger


SourceFiles
    effectivenessHeatExchangerSource.C

\*---------------------------------------------------------------------------*/

#ifndef effectivenessHeatExchangerSource_H
#define effectivenessHeatExchangerSource_H

#include "fvOption.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class effectivenessHeatExchangerSource Declaration
\*---------------------------------------------------------------------------*/

class effectivenessHeatExchangerSource
:
    public option
{

protected:

    // Protected data

        //- Secondary flow mass rate [Kg/s]
        scalar secondaryMassFlowRate_;

        //- Inlet secondary temperature [K]
        scalar secondaryInletT_;

        //- Primary air temperature at the heat exchanger inlet [K]
        scalar primaryInletT_;

        //- 2D look up table efficiency = F(total primary mass flow rate [Kg/s]
        //, secondary mass flow rate [Kg/s])
        autoPtr<interpolation2DTable<scalar> > eTable_;

        //- Name of velocity field; default = U
        word UName_;

        //- Name of temperature field; default = T
        word TName_;

        //- Name of the flux
        word phiName_;

        //- Name of the faceZone at the heat exchange inlet
        word faceZoneName_;

        //- Id for the face zone
        label zoneID_;

         //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;

        //- List of +1/-1 representing face flip map
        //  (1 use as is, -1 negate)
        labelList faceSign_;

        //- Area of the face zone
        scalar faceZoneArea_;


private:

    // Private Member Functions


        //- Disallow default bitwise copy construct
        effectivenessHeatExchangerSource
        (
            const effectivenessHeatExchangerSource&
        );

        //- Disallow default bitwise assignment
        void operator=(const effectivenessHeatExchangerSource&);

        //- Init heat exchanger source model
        void init();

        //- Calculate total area of faceZone accross procesors
        void calculateTotalArea(scalar& area);


public:

    //- Runtime type information
    TypeName("effectivenessHeatExchangerSource");


    // Constructors

        //- Construct from components
        effectivenessHeatExchangerSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~effectivenessHeatExchangerSource()
    {}


    // Member Functions


        // Public Functions

            //- Source term to fvMatrix<scalar>
            virtual void addSup(fvMatrix<scalar>& eqn, const label fieldI);


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
