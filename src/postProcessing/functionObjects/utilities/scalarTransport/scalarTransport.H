/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::scalarTransport

Group
    grpUtilitiesFunctionObjects

Description
    This function object evolves a passive scalar transport equation.  The
    field is initially zero, to which sources are added.  The field name
    is assigned the name of the function object.  Boundary conditions are
    read from file if the field file is available, or automatically applied
    based on the velocity boundary conditions.

    - the field can be zeroed on start-up using the resetOnStartUp flag
    - to employ the same numerical schemes as the flow velocity, use the
      autoSchemes flag
    - the diffusivity can be set manually using the DT entry, or retrieved
      from the turbulence model (if applicable)

    Example of function object specification to solve a scalar transport
    equation:
    \verbatim
    functions
    {
        scalar1                                     // <--- NOTE: SCALAR NAME
        {
            type            scalarTransport;
            functionObjectLibs ("libutilityFunctionObjects.so");
            resetOnStartUp  no;

            autoSchemes     yes;

            fvOptions
            {
                ...
            }
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: scalarTransport | yes      |
        phiName      | Name of flux field      | no          | phi
        UName        | Name of velocity field  | no          | U
        rhoName      | Name of density field   | no          | rho
        nCorr        | Number of correctors    | no          | 0
        DT           | Diffision coefficient   | no          | auto generated
        resetOnStartUp | Reset source to zero on start-up | no | no
        autoSchemes  | Use U schemes to evolve scalar | no   | yes
        fvOptions    | List of scalar sources  | yes         |
        log          | Log to standard output  | no          | yes
    \endtable

SeeAlso
    fvOption.H
    fvOptionList.H

SourceFiles
    scalarTransport.C
    IOscalarTransport.H


\*---------------------------------------------------------------------------*/

#ifndef scalarTransport_H
#define scalarTransport_H

#include "volFields.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldFwd.H"
#include "fvMatricesFwd.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class scalarTransport Declaration
\*---------------------------------------------------------------------------*/

class scalarTransport
{
    // Private data

        //- Name of this set of scalarTransport objects
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- On/off switch
        bool active_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of velocity field (optional)
        word UName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Diffusion coefficient (optional)
        scalar DT_;

        //- Flag to indicate whether user DT_ is used
        bool userDT_;

        //- Flag to reset scalar field on start-up
        bool resetOnStartUp_;

        //- Number of corrector iterations (optional)
        label nCorr_;

        //- Flag to employ schemes for velocity for scalar transport
        bool autoSchemes_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- Return the boundary types for the scalar field
        wordList boundaryTypes() const;

        //- Return reference to registered transported field
        volScalarField& transportedField();

        //- Return the diffusivity field
        tmp<volScalarField> DT(const surfaceScalarField& phi) const;

        //- Disallow default bitwise copy construct
        scalarTransport(const scalarTransport&);

        //- Disallow default bitwise assignment
        void operator=(const scalarTransport&);


public:

    //- Runtime type information
    TypeName("scalarTransport");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        scalarTransport
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~scalarTransport();


    // Member Functions

        //- Return name of the set of scalarTransport
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the scalarTransport data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the scalarTransport and write
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
