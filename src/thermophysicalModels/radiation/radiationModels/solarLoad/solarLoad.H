/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::radiation::solarLoad

Group
    grpRadiationModels

Description

    The solar load radiation model includes Sun primary hits, their
    reflective fluxes and diffusive sky radiative fluxes.

    The primary hit rays are calculated using a face shading algorithm.
    The reflected fluxes are considered diffusive and use a view factors method
    to deposit the energy on "visible" walls. The sky diffusive radiation for
    horizontal and vertical walls is calculated following the Fair Weather
    Conditions Method from the ASHRAE Handbook.

    By default the energy is included in cells adjacent to the patches into
    the energy Equation (wallCoupled = false). On coupled patches the flux is
    by default added to the wall and considered into the solid
    (solidCoupled = true).

    The reflected fluxes uses a grey absoprtion/emission model wich is weighted
    by the spectral distribution. The flag useVFbeamToDiffuse should be
    switched on and the view factors should be calculated using the
    'viewFactorsGen' application.

    The solarLoad model can be used in conjuntion with fvDOM and viewFactor
    radiation models but only using a single band spectrum. On the
    corresponding BC's for these models the flag "solarLoad" must be set to
    true.


SourceFiles
    solarLoad.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_solarLoad_H
#define radiation_solarLoad_H

#include "radiationModel.H"
#include "singleCellFvMesh.H"
#include "scalarListIOList.H"
#include "volFields.H"
#include "faceShading.H"
#include "solarCalculator.H"
#include "IOmapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class solarLoad Declaration
\*---------------------------------------------------------------------------*/

class solarLoad
:
    public radiationModel
{
protected:

    // Static data

        //- Static name for view factor walls
        static const word viewFactorWalls;


private:

    // Private data


        //- Agglomeration List
        labelListIOList finalAgglom_;

        //- Coarse mesh
        autoPtr<singleCellFvMesh> coarseMesh_;

        //- Net radiative heat flux [W/m2]
        volScalarField qr_;

        //- Secondary solar radiative heat flux [W/m2]
        volScalarField qsecondRad_;

        //- Direct hit faces Ids
        autoPtr<faceShading> hitFaces_;

        //- Constant source term
        DimensionedField<scalar, volMesh> Ru_;

        //- Solar calculator
        solarCalculator solarCalc_;

        //- Vertical direction (Default is g vector)
        vector verticalDir_;

        //- Include diffuse reflected heat fluxes from direct heat flux
        bool useVFbeamToDiffuse_;

        //- Selected patches to apply solar load
        labelList includePatches_;

        //- Cached coarse to fine mapping for coarse mesh
        List<labelListList> coarseToFine_;

        //-Number of bands
        label nBands_;

        //- Spectral distribution for the integrated solar heat flux
        scalarList spectralDistribution_;

        //- Map distribute
        autoPtr<IOmapDistribute> map_;

        //- Face-compact map
        labelListIOList visibleFaceFaces_;

        //- Couple solids through mapped boundary patch using qr (default:true)
        bool solidCoupled_;

        //- Couple wall patches using qr (default:false)
        bool wallCoupled_;

        //- Absorptivity list
        List<List<tmp<scalarField>>> absorptivity_;

        //- Update absorptivity
        bool updateAbsorptivity_;

        //- First iteration
        bool firstIter_;

        //- Update Sun position index
        label updateTimeIndex_;


    // Private Member Functions


        //- Initialise
        void initialise(const dictionary&);

        //- Update direct hit faces radiation
        void updateDirectHitRadiation(const labelList&, const labelHashSet&);

        //- Calculate diffusive heat flux
        void calculateQdiff(const labelHashSet&, const labelHashSet&);

        //- Update Sky diffusive radiation
        void updateSkyDiffusiveRadiation
        (
            const labelHashSet&,
            const labelHashSet&
        );

        //- Update hit faces
        bool updateHitFaces();

        //- Update absorptivity
        void updateAbsorptivity(const labelHashSet& includePatches);

        //- No copy construct
        solarLoad(const solarLoad&) = delete;

        //- No copy assignment
        void operator=(const solarLoad&) = delete;


public:

    //- Runtime type information
    TypeName("solarLoad");


    // Constructors

        //- Construct from volScalarField
        solarLoad(const volScalarField& T);

        //- Construct from dictionary and volScalarField
        solarLoad(const dictionary& dict, const volScalarField& T);

        //- Constructor from local components. Does not create a radiationModel.
        //  radWallFieldName is the solar heat field name
        solarLoad
        (
            const dictionary& dict,
            const volScalarField& T,
            const word radWallFieldName
        );


    //- Destructor
    virtual ~solarLoad();


    // Member functions

        // Edit

            //- Solve
            void calculate();

            //- Read radiation properties dictionary
            bool read();

            //- Source term component (for power of T^4)
            virtual tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            virtual tmp<DimensionedField<scalar, volMesh>> Ru() const;


    // Access

        //- Number of bands
        label nBands() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
