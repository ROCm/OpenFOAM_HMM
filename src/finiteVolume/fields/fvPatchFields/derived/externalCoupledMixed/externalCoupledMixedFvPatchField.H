/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::externalCoupledMixedFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition provides an interface to an external application.
    Values are transferred as plain text files, comprising of the constituent
    pieces of the `mixed' condition, i.e.

        <value1> <gradient1> <valueFracion1>
        <value2> <gradient2> <valueFracion2>
        <value3> <gradient3> <valueFracion3>
        ...
        <valueN> <gradientN> <valueFracionN>

    At start-up, the boundary creates a lock file, e.g.

        $FOAM_CASE/comms/patchName/OpenFOAM.lock

    ... to signal the external source to wait.  During the boundary condition
    update, boundary values are written to file, e.g.

        $FOAM_CASE/comms/patchName/data.out

    The lock file is then removed, instructing the exterbal source to take
    control of the program execution.  When ready, the external program
    should create the return values, e.g. to file

        $FOAM_CASE/comms/patchName/data.in

    ... and then re-instate the lock file.  The boundary condition will then
    read the return values, and pass program execution back to OpenFOAM.


    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        commsDir     | communications folder   | yes         |
        fileName     | data transfer file name | yes         |
        waitInterval | interval [s] between file checks | no | 1
        timeOut      | time after which error invoked | no   | 100*waitInterval
        calcFrequency | calculation frequency  | no          | 1
        log          | log program control     | no          | no
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            externalCoupled;
        commsDir        "$FOAM_CASE/comms";
        fileName        data;
    }
    \endverbatim

SeeAlso
    mixedFvPatchField

SourceFiles
    externalCoupledMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef externalCoupledMixedFvPatchField_H
#define externalCoupledMixedFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IFstream;

/*---------------------------------------------------------------------------*\
              Class externalCoupledMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class externalCoupledMixedFvPatchField
:
    public mixedFvPatchField<Type>
{

protected:

    // Private data

        //- Path to communications folder
        fileName commsDir_;

        //- Name of data file
        word fName_;

        //- Interval time between checking for return data [s]
        label waitInterval_;

        //- Time out time [s]
        label timeOut_;

        //- Calculation frequency
        label calcFrequency_;

        //- Log flag
        bool log_;


    // Private Member Functions

        //- Return the file path to the base communications folder
        fileName baseDir() const;

        //- Create lock file
        void createLockFile() const;

        //- Remove lock file
        void removeLockFile() const;

        //- Wait and remove lock file
        void writeAndWait(const fileName& transferFile) const;

        //- Initialise input stream for reading
        void initialiseRead(IFstream& is) const;


public:

    //- Runtime type information
    TypeName("externalCoupled");

    //- Name of lock file
    static word lockName;


    // Constructors

        //- Construct from patch and internal field
        externalCoupledMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        externalCoupledMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given externalCoupledMixedFvPatchField
        //  onto a new patch
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new externalCoupledMixedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        externalCoupledMixedFvPatchField
        (
            const externalCoupledMixedFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new externalCoupledMixedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "externalCoupledMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
