/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::ensightensightSurfaceReader

Description
    Ensight format surface reader

SourceFiles
    ensightSurfaceReader.C

\*---------------------------------------------------------------------------*/

#ifndef ensightSurfaceReader_H
#define ensightSurfaceReader_H

#include "surfaceReader.H"
#include "ensightReadFile.H"
#include "StringStream.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ensightSurfaceReader Declaration
\*---------------------------------------------------------------------------*/

class ensightSurfaceReader
:
    public surfaceReader
{
protected:

    // Protected Data

        //- Format flag
        IOstream::streamFormat streamFormat_;

        //- Base directory
        fileName baseDir_;

        //- Name of mesh file
        word meshFileName_;

        //- Field names
        List<word> fieldNames_;

        //- Field file names
        List<string> fieldFileNames_;

        //- Number of time steps
        label nTimeSteps_;

        //- Start time index
        label timeStartIndex_;

        //- Time increment
        label timeIncrement_;

        //- Times
        instantList timeValues_;

        //- Pointer to the surface
        autoPtr<meshedSurface> surfPtr_;

        List<Tuple2<string, label>> schema_;


    // Protected Member Functions

        //- Helper function to skip forward n steps in stream
        void skip(const label n, Istream& is) const;

        //- Helper function to read an ascii line from file
        void readLine(IFstream& is, string& buffer) const;

        //- Read and check a section header
        void debugSection(const word& expected, IFstream& is) const;

        //- Read (and throw away) geometry file header
        void readGeometryHeader(ensightReadFile& is) const;

        //- Read the case file
        void readCase(IFstream& is);

        //- Helper function to return Type after skipping n tokens
        template<class Type>
        void readFromLine
        (
            const label nSkip,
            IStringStream& is,
            Type& value
        ) const;

        //- Helper function to return Type after skipping n tokens
        template<class Type>
        void readFromLine
        (
            const label nSkip,
            const string& buffer,
            Type& value
        ) const;

        //- Helper function to return a field
        template<class Type>
        tmp<Field<Type>> readField
        (
            const label timeIndex,
            const label fieldIndex
        ) const;


public:

    //- Runtime type information
    TypeName("ensight");

    // Constructors

        //- Construct from fileName
        ensightSurfaceReader(const fileName& fName);


    //- Destructor
    virtual ~ensightSurfaceReader();


    // Member Functions

        //- Return a reference to the surface geometry
        virtual const meshedSurface& geometry();

        //- Return a list of the available times
        virtual instantList times() const;

        //- Return a list of the available fields at a given time
        virtual wordList fieldNames(const label timeIndex) const;

        //- Return a scalar field at a given time
        virtual tmp<Field<scalar>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const scalar& refValue = pTraits<scalar>::zero
        ) const;

        //- Return a scalar field at a given time
        virtual tmp<Field<vector>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const vector& refValue = pTraits<vector>::zero
        ) const;

        //- Return a sphericalTensor field at a given time
        virtual tmp<Field<sphericalTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const sphericalTensor& reValue = pTraits<sphericalTensor>::zero
        ) const;

        //- Return a symmTensor field at a given time
        virtual tmp<Field<symmTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const symmTensor& reValue = pTraits<symmTensor>::zero
        ) const;

        //- Return a tensor field at a given time
        virtual tmp<Field<tensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const tensor& reValue = pTraits<tensor>::zero
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightSurfaceReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
