/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::fv::effectivenessHeatExchangerSource

Group
    grpFvOptionsSources

Description
    Heat exchanger source model for compressible flows, where the heat
    exchanger is modelled as an energy source using a selection of cells.

    The total heat exchange source is given by:
    \f[
        Q_t = e(\phi, \dot{m}_2) (T_2 - T_1) \phi c_p
    \f]

    where:
    \vartable
        Q_t       | Total heat exchange source [J/s]
        e(\phi,\dot{m}_2) | Effectivenes table [-]
        \phi      | Net mass flux entering heat exchanger [kg/s]
        \dot{m}_2 | Secondary flow mass flow rate [kg/s]
        T_1       | Primary flow inlet temperature [K]
        T_2       | Secondary flow inlet temperature [K]
        c_p       | Primary flow specific heat capacity [J/kg/K]
    \endvartable


    The distribution inside the heat exchanger is given by:
    \f[
        Q_c = \frac{V_c |U_c| (T_c - T_{ref})}{\sum(V_c |U_c| (T_c - T_{ref}))}
    \f]

    where:
    \vartable
        Q_c     | Source for cell
        V_c     | Volume of the cell [m3]
        U_c     | Local cell velocity [m/s]
        T_c     | Local cell temperature [K]
        T_{ref} | Min or max(T) in cell zone depending on the sign of Qt [K]
    \endvartable

    Sources applied to either of the below, if exist:
    \verbatim
      e         | Internal energy                            [m2/s2]
      h         | Enthalphy                                  [m2/s2]
    \endverbatim

    Required fields:
    \verbatim
      T         | Temperature                                [K]
      U         | Velocity                                   [m/s]
      phi       | Mass flux                                  [kg/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    effectivenessHeatExchangerSource1
    {
        // Mandatory entries
        type                     effectivenessHeatExchangerSource;
        faceZone                 <faceZoneName>;
        secondaryMassFlowRate    <scalar>;
        secondaryInletT          <scalar>;
        file                     "effectivenessTable";
        outOfBounds              clamp;

        // Optional entries
        U                        <word>;
        T                        <word>;
        phi                      <word>;

        // Conditional optional entries

            // when the total heat exchange is calculated with primary inlet T
            primaryInletT           <scalar>;

            // when the total heat exchange is calculated with a given target
            targetQdot              <scalar>;
            targetQdotCalcInterval  <label>;
            targetQdotRelax         <scalar>;

            // when secondary outlet temperature is requested
            secondaryCp             <Function1<scalar>>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type  | Reqd | Deflt
      type      | Type name: effectivenessHeatExchangerSource <!--
                -->                                 | word  | yes  | -
      secondaryMassFlowRate | Secondary flow mass rate [kg/s] <!--
                -->                                 | scalar | yes | -
      secondaryInletT | Secondary flow inlet temperature [K]  <!--
                -->                                 | scalar | yes | -
      faceZone  | Name of the faceZone at the heat exchanger inlet <!--
                -->                                 | word   | yes | -
      file      | 2D effectiveness table = function of primary    <!--
                --> and secondary mass flow rates [kg/s] | file | yes | -
      primaryInletT | Primary flow temperature at the heat exchanger inlet <!--
                -->                                 | scalar | no  | -
      targetQdot |  Target heat rejection           | scalar | no  | -
      targetQdotCalcInterval | Target heat rejection calculation interval <!--
                -->                                 | label  | no  | -
      targetQdotRelax | Target heat rejection temperature <!--
                --> under-relaxation coefficient    | scalar | no  | -
      secondaryCp | Secondary flow specific heat capacity <!--
                -->                    | Function1\<scalar\> | no  | -
      U         | Name of operand velocity field    | word   | no  | U
      T         | Name of operand temperature field | word   | no  | T
      phi       | Name of operand flux field        | word   | no  | phi
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink
      - \link Function1.H \endlink

    The effectiveness table is described in terms of the primary and secondary
    mass flow rates.  For example, the table:

    \verbatim
                           secondary MFR
                       |  0.1   0.2   0.3
                  -----+-----------------
                  0.02 |   A     B     C
     primary MFR  0.04 |   D     E     F
                  0.06 |   G     H     I
    \endverbatim

    is specified by the following:

    \verbatim
        (
            (
                0.02
                (
                    (0.1    A)
                    (0.2    B)
                    (0.3    C)
                )
            )
            (
                0.04
                (
                    (0.1    D)
                    (0.2    E)
                    (0.3    F)
                )
            )
            (
                0.06
                (
                    (0.1    G)
                    (0.2    H)
                    (0.3    I)
                )
            )
        );
    \endverbatim

Note
  - Primary flow indicates the CFD flow region and
    secondary flow the non-CFD-model region.
  - The table with name \c file should have the same units as the
    secondary mass flow rate and kg/s for \c phi.
  - \c faceZone is the faces at the inlet of the \c cellZone, it needs to be
    created with flip map flags. It is used to integrate the net mass flow
    rate into the heat exchanger.
  - \c primaryInletT sets the primary inlet temperature.  If not set, the
    flux-averaged temperature is used.

SourceFiles
    effectivenessHeatExchangerSource.C

\*---------------------------------------------------------------------------*/

#ifndef fv_effectivenessHeatExchangerSource_H
#define fv_effectivenessHeatExchangerSource_H

#include "cellSetOption.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"
#include "writeFile.H"
#include "Function1.H"

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
    public fv::cellSetOption,
    public functionObjects::writeFile
{
    // Private Data

        //- Flag to use a user-specified primary flow inlet temperature
        bool userPrimaryInletT_;

        //- Flag to use target heat rejection
        bool targetQdotActive_;

        //- Secondary flow specific heat capacity [J/kg/K]
        autoPtr<Function1<scalar>> secondaryCpPtr_;

        //- 2D effectiveness table = function of primary and secondary
        //- mass flow rates [kg/s]
        autoPtr<interpolation2DTable<scalar>> eTable_;

        //- Target heat rejection calculation interval
        label targetQdotCalcInterval_;

        //- Secondary flow mass rate [kg/s]
        scalar secondaryMassFlowRate_;

        //- Secondary flow inlet temperature [K]
        scalar secondaryInletT_;

        //- Primary flow temperature at the heat exchanger inlet [K]
        scalar primaryInletT_;

        //- Target heat rejection
        scalar targetQdot_;

        //- Target heat rejection temperature under-relaxation coefficient
        scalar targetQdotRelax_;

        //- Name of operand velocity field
        word UName_;

        //- Name of operand temperature field
        word TName_;

        //- Name of operand flux field
        word phiName_;

        //- Name of the faceZone at the heat exchanger inlet
        word faceZoneName_;

        //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch IDs per face
        labelList facePatchId_;

        //- List of +1/-1 representing face flip map (1 use as is, -1 negate)
        labelList faceSign_;


    // Private Member Functions

        //- Initialise heat exchanger source model
        void initialise();

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


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

        //- No copy construct
        effectivenessHeatExchangerSource
        (
            const effectivenessHeatExchangerSource&
        ) = delete;

        //- No copy assignment
        void operator=(const effectivenessHeatExchangerSource&) = delete;


    //- Destructor
    virtual ~effectivenessHeatExchangerSource() = default;


    // Member Functions

        //- Add explicit/implicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        )
        {
            NotImplemented;
        }

        //- Add explicit/implicit contribution
        //- to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


        // I-O

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
