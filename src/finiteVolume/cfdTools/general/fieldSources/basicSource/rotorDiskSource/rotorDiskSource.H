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
    Foam::rotorDiskSource

Description
    Cell based momemtum source

    Source approximates the mean effects of rotor forces on a cylindrical
    region within the domain

    Sources described by:

        rotorDiskSourceCoeffs
        {
            fieldName       U;      // name of field on which to apply source
            rhoName         rho;    // density field if compressible case
            nBlades         3;      // number of blades
            tip effect      0.96;   // normalised radius above which lift = 0

            inletFlowType   local;  // inlet flow type specification

            geometryMode    auto;   // geometry specification

            refDirection    (-1 0 0); // reference direction

            flapCoeffs
            {
                beta0           0;  // coning angle [deg]
                beta1           0;  // lateral flapping coeff
                beta2           0;  // longitudinal flapping coeff
            }
            trimCoeffs
            {
                alphac          15; // collective pitch angle [deg]
                A               0;  // lateral cyclic coeff
                B               0;  // longitudinal cyclic coeff
            }
            blade
            {
                ...
            }
            profiles
            {
                ...
            }
        }

    Where:

        geometryMode =
            auto          : determine rototor co-ord system from cells
            specified     : specified co-ord system

        inletFlowType =
            fixed         : specified velocity
            surfaceNormal : specified normal velocity (positive towards rotor)
            local         : use local flow conditions



SourceFiles
    rotorDiskSource.C
    rotorDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef rotorDiskSource_H
#define rotorDiskSource_H

#include "basicSource.H"
#include "cylindricalCS.H"
#include "NamedEnum.H"
#include "bladeModel.H"
#include "profileModelList.H"
#include "volFieldsFwd.H"
#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rotorDiskSource Declaration
\*---------------------------------------------------------------------------*/

class rotorDiskSource
:
    public basicSource
{
public:

    enum geometryModeType
    {
        gmAuto,
        gmSpecified
    };
    static const NamedEnum<geometryModeType, 2> geometryModeTypeNames_;

    enum inletFlowType
    {
        ifFixed,
        ifSurfaceNormal,
        ifLocal
    };
    static const NamedEnum<inletFlowType, 3> inletFlowTypeNames_;


protected:

    // Helper structures to encapsulate flap and trim data

        struct flapData
        {
            scalar beta0;   // coning angle [deg]
            scalar beta1;   // lateral flapping coeff
            scalar beta2;   // longitudinal flapping coeff
        };

        struct trimData
        {
            scalar alphaC;  // collective pitch angle [deg]
            scalar A;       // lateral cyclic coeff
            scalar B;       // longitudinal cyclic coeff
        };


    // Protected data

        //- Name of field on which to apply the source
        word fieldName_;

        //- Name of density field
        word rhoName_;

        //- Rotational speed [rad/s]
        scalar omega_;

        //- Number of blades
        label nBlades_;

        //- Inlet flow type
        inletFlowType inletFlow_;

        //- Inlet velocity for specified iinflow
        vector inletVelocity_;

        //- Tip effect [0-1]
        //  Ratio of blade radius beyond which lift=0
        scalar tipEffect_;

        //- Blade flap coefficients [rad/s]
        flapData flap_;

        //- Blad trim coefficients
        trimData trim_;

        //- Blade data
        bladeModel blade_;

        //- Profile data
        profileModelList profiles_;

        //- Cell centre positions in local rotor frame (Cartesian x, y, z)
        List<point> x_;

        //- Rotation tensor for flap angle
        List<tensor> R_;

        //- Geometric angle of attack [deg]
        List<scalar> alphag_;

        //- Area [m2]
        List<scalar> area_;

        //- Rotor co-ordinate system (r, theta, z)
        cylindricalCS coordSys_;

        //- Maximum radius
        scalar rMax_;

        
    // Protected Member Functions

        //- Check data
        void checkData();

        //- Set the face areas per cell, and optionally correct the rotor axis
        void setFaceArea(vector& axis, const bool correct);

        //- Create the co-ordinate system
        void createCoordinateSystem();

        //- Construct geometry
        void constructGeometry();

        //- Return the inlet flow field
        tmp<vectorField> inflowVelocity(const volVectorField& U) const;

        //- Calculate forces
        template<class RhoType>
        tmp<volVectorField> calculateForces
        (
            const RhoType& rho,
            const vectorField& U,
            const dimensionSet& dims
        );

        //- Helper function to write rotor values
        template<class Type>
        void writeField
        (
            const word& name,
            const List<Type>& values,
            const bool writeNow = false
        ) const;


public:

    //- Runtime type information
    TypeName("rotorDisk");


    // Constructors

        
        //- Construct from components
        rotorDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~rotorDiskSource();


    // Member Functions


        // Checks

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;


        // Source term addition

            //- Source term to fvMatrix<vector>
            virtual void addSup(fvMatrix<vector>& eqn, const label fieldI);


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rotorDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

