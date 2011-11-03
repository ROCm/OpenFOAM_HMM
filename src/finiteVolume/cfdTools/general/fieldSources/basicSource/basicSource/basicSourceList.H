/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::basicSourceList

Description
    List of explict sources

SourceFile
    basicSourceList.C

\*---------------------------------------------------------------------------*/

#ifndef basicSourceList_H
#define basicSourceList_H

#include "PtrList.H"
#include "DimensionedField.H"
#include "basicSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class basicSourceList Declaration
\*---------------------------------------------------------------------------*/

class basicSourceList
:
    public PtrList<basicSource>
{
private:

    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Time index to check that all defined sources have been applied
        label checkTimeIndex_;


    // Private Member Functions

        //- Check that all sources have been applied
        void checkApplied() const;

        //- Disallow default bitwise copy construct
        basicSourceList(const basicSourceList&);

        //- Disallow default bitwise assignment
        void operator=(const basicSourceList&);


public:

    //- Runtime type information
    TypeName("basicSourceList");


    // Constructors

        //- Construct from components with list of field names
        basicSourceList(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~basicSourceList()
    {}


    // Member Functions

        // Evaluation

            //- Apply sources to equation
            template<class Type>
            void apply(fvMatrix<Type>& eqn);

            //- Apply sources to equation with specified name
            template<class Type>
            void apply(fvMatrix<Type>& eqn, const word& name);


        // I-O

            //- Read dictionary
            virtual bool read(const dictionary& dict);

            //- Write data to Istream
            virtual bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const basicSourceList& sources
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicSourceListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
