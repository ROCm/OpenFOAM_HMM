/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    Foam::fv::optionAdjoint

Description
    Similar to fv::option but with additional functionality to contribute to
    the sensitivity deriavtives

SourceFiles
    fvOptionAdjoint.C

\*---------------------------------------------------------------------------*/

#ifndef fvOptionAdjoint_H
#define fvOptionAdjoint_H

#include "fvOption.H"
#include "incompressibleAdjointVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class optionAdjoint Declaration
\*---------------------------------------------------------------------------*/

class optionAdjoint
:
    public option
{
public:

    //- Runtime type information
    TypeName("optionAdjoint");


    // Constructors

        //- Construct from components
        optionAdjoint
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<optionAdjoint> clone() const
        {
            notImplemented("autoPtr<optionAdjoint> clone() const");
            return autoPtr<optionAdjoint>(nullptr);
        }


     // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            optionAdjoint,
            dictionary,
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (name, modelType, dict, mesh)
        );


    // Selectors

        //- Return a reference to the selected fvOption model
        static autoPtr<optionAdjoint> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~optionAdjoint() = default;


    // Member Functions
    virtual tmp<volVectorField> dxdbMult(const incompressibleAdjointVars&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
