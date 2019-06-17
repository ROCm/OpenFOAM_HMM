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
    Foam::fv::optionAdjointList

DescriptionAdjoint
    List of finite volume optionAdjoints

SourceFile
    optionAdjointList.C

\*---------------------------------------------------------------------------*/

#ifndef optionAdjointList_H
#define optionAdjointList_H

#include "PtrList.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "fvOptionAdjoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class optionAdjointList Declaration
\*---------------------------------------------------------------------------*/

class optionAdjointList
:
    public PtrList<optionAdjoint>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        optionAdjointList(const optionAdjointList&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const optionAdjointList&) = delete;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Time index to check that all defined sources have been applied
        label checkTimeIndex_;


    // Protected Member Functions

        //- Return the "optionAdjoints" sub-dictionary if present otherwise
        //- return dict
        const dictionary& optionAdjointsDict(const dictionary& dict) const;

        //- Read optionAdjoints dictionary
        bool readOptionAdjoints(const dictionary& dict);

        //- Check that all sources have been applied
        void checkApplied() const;


public:

    //- Runtime type information
    TypeName("optionAdjointList");


    // Constructors

        //- Construct null
        optionAdjointList(const fvMesh& mesh);

        //- Construct from mesh and dictionary
        optionAdjointList(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~optionAdjointList() = default;


    // Member Functions

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Correct
        template<class Type>
        void correct(GeometricField<Type, fvPatchField, volMesh>& fld);


        // Sources

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& fld
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& fld,
                const word& fieldName
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& fld
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& fld,
                const word& fieldName
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& fld
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& fld,
                const word& fieldName
            );


        // Constraints

            //- Apply constraints to equation
            template<class Type>
            void constrain(fvMatrix<Type>& eqn);


        // I-O

            //- Read dictionary
            virtual bool read(const dictionary& dict);

            //- Write data to Ostream
            virtual bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const optionAdjointList& optionAdjoints
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvOptionAdjointListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
