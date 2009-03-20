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
    Foam::radiation::fvDOM

Description

    Finite Volume Discrete Ordinary Method. It solves the RTE equation for n
    directions in a participating media. It does not consider scatter.

    Available absorption models:
                                greyMeanAbsoprtionEmission
                                wideBandAbsorptionEmission

    i.e. dictionary
    fvDOMCoeffs
    {
        Nphi    1;          // azimuthal angles in PI/2 on X-Y.(from Y to X)
        Ntheta  2;          // polar angles in P1 (from Z to X-Y plane)
        convergence 1e-4;   //convergence criteria for radiation iteration
    }

    nFlowIterPerRadIter   1; // Number of flow iterations per radiation
                                iteration

    The total number of solid angles is  4 * Nphi * Ntheta.

    In 1D the direction of the rays is X (Nphi and Ntheta are ignored)
    In 2D the direction of the rays is on X-Y plane (only Nphi is considered)
    In 3D (Nphi and Ntheta are considered)

SourceFiles
    fvDOM.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelfvDOM_H
#define radiationModelfvDOM_H

#include "radiativeIntensityRay.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class fvDOM Declaration
\*---------------------------------------------------------------------------*/

class fvDOM
:
    public radiationModel
{

    // Private data

        //- Incident radiation  [W/m2]
        volScalarField G_;

        //- Total Radiative heat flux [W/m2]
        volScalarField Qr_;

        //- Total Absorption coefficient  [1/m]
        volScalarField a_;

        //- wavelength total Absorption coefficient [1/m]
        PtrList<volScalarField> aj_;

        //- Total Emission coefficient [1/m]
        volScalarField e_;

        //- Emission contribution [Kg/m/s^3]
        volScalarField E_;

        //- Number of solid angles in theta
        label Ntheta_;

        //- Number of solid angles in phi
        label Nphi_ ;

        //- Total number of directions
        label Ni_;

        //- Number of wavelength bands
        label lambdaj_;

        //- Black Body
        blackBodyEmission blackBody_;

        //- List of Pointers to RadiativeIntensityRay
        PtrList<radiativeIntensityRay> RadIntRay_;


    // Private member functions

        //- Disallow default bitwise copy construct
        fvDOM(const fvDOM&);

        //- Disallow default bitwise assignment
        void operator=(const fvDOM&);

        //- Update Absorption Coefficients
//        void updateAbsorptionCoeffs(void);

        //- Update Black Body Emissiom
        void updateBlackBodyEmission(void);

public:

    // Static data members

        static label iterRadId;

    //- Runtime type information
        TypeName("fvDOM");


    // Constructors

        //- Construct from components
        fvDOM(const volScalarField& T);


    // Destructor

        ~fvDOM();


    // Member functions

        // Edit

        //- Update radiationSource varible
        void correct();

        //- Read radiationProperties dictionary
        bool read();

        //- Update G and calculate total heat flux on boundary
        void updateG();

        //- Source term component (for power of T^4)
        virtual tmp<volScalarField> Rp() const;

        //- Source term component (constant)
        virtual tmp<DimensionedField<scalar, volMesh> > Ru() const;

        // Access

        //- Intensity Ray on i direction
        inline const radiativeIntensityRay& RadIntRay(label i) const
        {
            return  RadIntRay_[i];
        }

        //- Intensity Ray on i direction and j band-width
        inline const volScalarField& RadIntRayiLambdaj
        (
            const label i,
            const label j
        ) const
        {
            return  RadIntRay_[i].Ilambdaj(j);
        }

        //-Number of angles in theta
        label Ntheta() const
        {
            return Ntheta_;
        }

        //- Number of angles in phi
        label Nphi() const
        {
            return Nphi_;
        }

        //- Number of directions
        label Ni() const
        {
            return Ni_;
        }

        //- Number of wavelengths
        inline const label& lambdaj() const
        {
            return lambdaj_;
        }

        // Const access to a
        inline const volScalarField& a() const
        {
            return a_;
        }

        // Const access to aj
        inline const volScalarField& aj(label i) const
        {
            return aj_[i];
        }

        // Const access to G
        inline const volScalarField& G() const
        {
            return G_;
        }

        //  Const access to Qr
        inline const volScalarField& Qr() const
        {
            return Qr_;
        }

        //  Const access to blavkBody
        virtual const blackBodyEmission& blackBody() const
        {
            return blackBody_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
