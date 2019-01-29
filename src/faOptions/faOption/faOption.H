/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::fa::option

Description
    Base abstract class for handling finite area options (i.e. \c faOption).

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    <userDefinedName1>
    {
        // Mandatory entries (unmodifiable)
        type              <faOptionName>;

        // Mandatory entries (runtime modifiable)
        region            <regionName>;

        // Optional entries (unmodifiable/runtime modifiable)
        <faOption>Coeffs
        {
            // subdictionary entries
        }

        // Optional entries (runtime modifiable)
        active            true;
        log               true;
    }
    \endverbatim

    where the entries mean:
    \table
      Property         | Description                   | Type  | Reqd | Dflt
      type             | Name of operand faOption      | word  | yes  | -
      region           | Name of operand region        | word  | yes  | -
      \<faOption\>Coeffs | Dictionary containing settings of <!--
                    --> the selected faOption settings | dictionary | no | -
      active   | Flag to (de)activate faOption         | bool  | no   | true
      log      | Flag to log faOption-related info     | bool  | no   | true
    \endtable

Note
  - Source/sink options are to be added to the right-hand side of equations.

SourceFiles
    faOption.C
    faOptionIO.C

\*---------------------------------------------------------------------------*/

#ifndef faOption_H
#define faOption_H

#include "faMatrices.H"
#include "areaFields.H"
#include "dictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "volSurfaceMapping.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class option Declaration
\*---------------------------------------------------------------------------*/

class option
{
    // Private Member Functions

        //- Construct region mesh and fields
        void constructMeshObjects();


protected:

    // Protected Data

        //- Source name
        const word name_;

        //- Model type
        const word modelType_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the patch
        const fvPatch& patch_;

        //- Top level source dictionary
        dictionary dict_;

        //- Dictionary containing source coefficients
        dictionary coeffs_;

        //- Source active flag
        Switch active_;

        //- Field names to apply source to - populated by derived models
        wordList fieldNames_;

        //- Applied flag list - corresponds to each fieldNames_ entry
        List<bool> applied_;

         //- Region name
        word regionName_;

        //- Pointer to the region mesh database
        autoPtr<faMesh> regionMeshPtr_;

        //-Volume-to surface mapping
        autoPtr<volSurfaceMapping> vsmPtr_;


public:

    //- Runtime type information
    TypeName("option");


     // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            option,
            dictionary,
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const fvPatch& patch
            ),
            (name, modelType, dict, patch)
        );


    // Constructors

        //- Construct from components
        option
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvPatch& patch
        );

        //- Return clone
        autoPtr<option> clone() const
        {
            NotImplemented;
            return nullptr;
        }

        //- Return pointer to new faOption object created
        //- on the freestore from an Istream
        class iNew
        {

            //- Reference to the patch
            const fvPatch& patch_;

            //- Name
            const word& name_;

        public:

            iNew
            (
                const fvPatch& patch,
                const word& name
            )
            :
                patch_(patch),
                name_(name)
            {}

            autoPtr<option> operator()(Istream& is) const
            {
                const dictionary dict(is);

                return autoPtr<option>
                (
                    option::New(name_, dict, patch_)
                );
            }
        };


    // Selectors

        //- Return a reference to the selected faOption model
        static autoPtr<option> New
        (
            const word& name,
            const dictionary& dict,
            const fvPatch& patch
        );


    //- Destructor
    virtual ~option() = default;


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return const access to fvPatch
            inline const fvPatch& patch() const;

            //- Return dictionary
            inline const dictionary& coeffs() const;

            //- Return const access to the source active flag
            inline bool active() const;

            //- Set the applied flag to true for field index fieldi
            inline void setApplied(const label fieldi);

            //- Return the region mesh database
            inline const faMesh& regionMesh() const;

             //- Return volSurfaceMapping
            inline const volSurfaceMapping& vsm() const;

            //- Region name
            inline const word& regionName() const;


        // Edit

            //- Return access to the source active flag
            inline Switch& active();


        // Checks

            //- Is the source active?
            virtual bool isActive();

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;

            //- Check that the source has been applied
            virtual void checkApplied() const;


        // Evaluation

            // Explicit and implicit sources

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Explicit and implicit sources for compressible equations

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Constraints

                virtual void constrain
                (
                    faMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Correction

                virtual void correct(areaScalarField& field);
                virtual void correct(areaVectorField& field);
                virtual void correct(areaSphericalTensorField& field);
                virtual void correct(areaSymmTensorField& field);
                virtual void correct(areaTensorField& field);


        // IO

            //- Write the source header information
            virtual void writeHeader(Ostream&) const;

            //- Write the source footer information
            virtual void writeFooter(Ostream&) const;

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
