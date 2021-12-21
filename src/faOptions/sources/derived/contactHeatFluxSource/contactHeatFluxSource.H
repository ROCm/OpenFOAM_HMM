/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::fa::contactHeatFluxSource

Group
    grpFaOptionsSources

Description
    Applies contact heat flux between specified \c faMesh
    and \c fvMesh within a specified region for compressible flows.

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    contactHeatFluxSource1
    {
        // Mandatory entries (unmodifiable)
        type                contactHeatFluxSource;
        Tprimary            <TprimaryFieldName>;

        // Optional entries (runtime modifiable)
        T                   <Tname>;
        thicknessLayers     (<layer1> <layer2> ... <layerN>);

        // Conditional optional entries (runtime modifiable)

            // when the entry "thicknessLayers" is present
            kappaLayers     (<layer1> <layer2> ... <layerN>);

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type  | Reqd | Dflt
      type      | Type name: contactHeatFluxSource  | word  | yes  | -
      Tprimary  | Name of primary temperature field | word  | yes  | -
      T         | Name of operand temperature field | word  | no   | T
      thicknessLayers | List of thicknesses of layers | scalarList | no | -
      kappaLayers | List of conductivities of layers | scalarList | cndtnl | -
    \endtable

    The inherited entries are elaborated in:
     - \link faOption.H \endlink
     - \link faceSetOption.H \endlink
     - \link temperatureCoupledBase.H \endlink

SourceFiles
    contactHeatFluxSource.C

\*---------------------------------------------------------------------------*/

#ifndef fa_contactHeatFluxSource_H
#define fa_contactHeatFluxSource_H

#include "faOption.H"
#include "Function1.H"
#include "areaFields.H"
#include "faceSetOption.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class contactHeatFluxSource Declaration
\*---------------------------------------------------------------------------*/

class contactHeatFluxSource
:
    public fa::faceSetOption,
    public temperatureCoupledBase
{
    // Private Data

        //- Name of temperature field
        word TName_;

        //- Name of primary temperature field
        word TprimaryName_;

        //- Primary region temperature
        const volScalarField& Tp_;

        //- Temperature - wall [K]
        areaScalarField Tw1_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Total contact resistance
        scalar contactRes_;

        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- Return htc from the primary region
        tmp<areaScalarField> htc() const;


public:

    //- Runtime type information
    TypeName("contactHeatFluxSource");


    // Constructors

        //- Construct from explicit source name and mesh
        contactHeatFluxSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvPatch& patch
        );

        //- No copy construct
        contactHeatFluxSource(const contactHeatFluxSource&) = delete;

        //- No copy assignment
        void operator=(const contactHeatFluxSource&) = delete;


    //- Destructor
    virtual ~contactHeatFluxSource() = default;


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper& mapper
            )
            {
                temperatureCoupledBase::autoMap(mapper);
            }

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>& fld,
                const labelList& map
            )
            {
                temperatureCoupledBase::rmap(fld, map);
            }


        // Evaluation

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const areaScalarField& h,
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
