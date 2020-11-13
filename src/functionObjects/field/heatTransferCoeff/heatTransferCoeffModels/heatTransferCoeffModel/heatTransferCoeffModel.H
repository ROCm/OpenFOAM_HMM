/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

Namespace
    Foam::heatTransferCoeffModels

Description
    A namespace for various heat transfer coefficient model implementations.

Class
    Foam::heatTransferCoeffModel

Description
    An abstract base class for heat transfer coeffcient models.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Mandatory and other optional entries
        ...

        // Mandatory (inherited) entries (runtime modifiable)
        patches     (<patch1> <patch2> ... <patchN>);

        // Optional (inherited) entries (runtime modifiable)
        qr          <qrName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                  | Type | Reqd | Dflt
      patches   | Names of operand patches  | wordRes | yes  | -
      qr        | Name of radiative heat flux  | word | no   | qr
    \endtable

SourceFiles
  - heatTransferCoeffModel.C
  - heatTransferCoeffModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModel_H
#define heatTransferCoeffModel_H

#include "dictionary.H"
#include "HashSet.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class heatTransferCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferCoeffModel
{
protected:

    // Protected Data

        //- Mesh reference
        const fvMesh& mesh_;

        //- Optional list of (wall) patches to process
        labelHashSet patchSet_;

        //- Temperature name
        const word TName_;

        //- Name of radiative heat flux
        word qrName_;


    // Protected Member Functions

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        ) = 0;


public:

    //- Runtime type information
    TypeName("heatTransferCoeffModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferCoeffModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word& TName
            ),
            (dict, mesh, TName)
        );


    // Selectors

        //- Return a reference to the selected heat transfer coefficient model
        static autoPtr<heatTransferCoeffModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );


    // Constructors

        //- Construct from components
        heatTransferCoeffModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );

        //- No copy construct
        heatTransferCoeffModel(const heatTransferCoeffModel&) = delete;

        //- No copy assignment
        void operator=(const heatTransferCoeffModel&) = delete;


    //- Destructor
    virtual ~heatTransferCoeffModel() = default;


    // Member Functions

        //- The mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Wall patches to process
        const labelHashSet& patchSet() const
        {
            return patchSet_;
        }

        //- Temperature name
        const word& TName() const
        {
            return TName_;
        }

        //- Name of radiative heat flux
        const word& qrName() const
        {
            return qrName_;
        }

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        virtual bool calc
        (
            volScalarField& result,
            const FieldField<Field, scalar>& q
        );

        //- Return q boundary fields
        tmp<FieldField<Field, scalar>> q() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
