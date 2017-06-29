/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::vtk::formatter

Description
    Abstract class for a VTK output stream formatter.

    Includes very simple support for writing XML elements.

SourceFiles
    foamVtkFormatter.C
    foamVtkFormatterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkFormatter_H
#define foamVtkFormatter_H

#include "int.H"
#include "uint64.H"
#include "label.H"
#include "word.H"
#include "UList.H"
#include "LIFOStack.H"
#include "foamVtkCore.H"
#include "foamVtkPTraits.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

class outputOptions;

/*---------------------------------------------------------------------------*\
                          Class formatter Declaration
\*---------------------------------------------------------------------------*/

class formatter
{
    // Private Data

        //- The output stream for the formatter
        std::ostream& os_;

        //- Stack of current XML tags
        LIFOStack<word> xmlTags_;

        //- Tag open/closed/ended state
        mutable bool inTag_;


    // Private Member Functions

        //- Write XML attribute key/value pair
        template<class Type>
        formatter& writeAttribute
        (
            const word& k,
            const Type& v,
            const char quote = '\''
        );


protected:

    // Protected Member Functions

        //- Construct and attach to an output stream
        inline formatter(std::ostream& os);

public:

    // Public typedefs

        //- Use UInt64 for header data
        typedef uint64_t headerType;


    //- Destructor
    virtual ~formatter();


    // Member Functions

      // Access

        //- Access to the underlying output stream
        inline std::ostream& os();

        //- The format-type / output options.
        virtual const outputOptions& opts() const = 0;

        //- Name for the XML output type or the legacy output type.
        virtual const char* name() const = 0;

        //- Name for the XML append encoding
        virtual const char* encoding() const = 0;


        //- Write leading size for binary output
        virtual void writeSize(const uint64_t nBytes) = 0;

        virtual void write(const uint8_t val) = 0;
        virtual void write(const label val)   = 0;
        virtual void write(const float val)   = 0;
        virtual void write(const double val)  = 0;

        //- Flush encoding, write newline etc.
        virtual void flush() = 0;

        //- The encoded length for binary output.
        //  The default is pass-through.
        virtual std::size_t encodedLength(std::size_t n) const;


    // Member Functions

      // Output

        //- Indent according to the currently nested XML tags
        void indent();

        //- Write XML header
        //  \return formatter for chaining
        formatter& xmlHeader();

        //- Write XML comment (at the current indentation level)
        //  \return formatter for chaining
        formatter& xmlComment(const std::string& comment);


        //- Open XML tag
        //  \return formatter for chaining
        formatter& openTag(const word& tagName);

        //- Open XML tag
        //  \return formatter for chaining
        inline formatter& openTag(const vtk::fileTag& tagEnum);

        //- Close XML tag, optional as an empty container.
        //  Always adds a trailing newline.
        //  \return formatter for chaining
        formatter& closeTag(const bool isEmpty = false);

        //- End XML tag, optional with sanity check
        //  Always adds a trailing newline.
        //  \return formatter for chaining
        formatter& endTag(const word& tagName = word::null);

        //- End XML tag with sanity check
        //  Always adds a trailing newline.
        //  \return formatter for chaining
        inline formatter& endTag(const vtk::fileTag& tagEnum);

        //- Write XML tag without any attributes. Combines openTag/closeTag.
        //  \return formatter for chaining
        inline formatter& tag(const word& tagName);

        //- Write XML tag without any attributes. Combines openTag/closeTag.
        //  \return formatter for chaining
        inline formatter& tag(const vtk::fileTag& tagEnum);

        //- Add a "VTKFile" XML tag for contentType, followed by a tag for
        //  the contentType itself. Optionally leave the contentType tag
        //  open for adding additional attributes.
        //  \return formatter for chaining
        formatter& beginVTKFile
        (
            const word& contentType,
            const word& contentVersion,
            const bool leaveOpen = false
        );

        //- Add a "VTKFile" XML tag for contentType, followed by a tag for
        //  the contentType itself. Optionally leave the contentType tag
        //  open for adding additional attributes.
        //  \return formatter for chaining
        inline formatter& beginVTKFile
        (
            const vtk::fileTag& contentType,
            const word& contentVersion,
            const bool leaveOpen = false
        );

        //- Add a "AppendedData" XML tag with the current encoding and output
        //  the requisite '_' prefix.
        //  \return formatter for chaining
        formatter& beginAppendedData();


        //- Open "DataArray" XML tag
        //  \return formatter for chaining
        template<class Type, int nComp=0>
        formatter& openDataArray(const word& dataName);

        //- Open "DataArray" XML tag
        //  \return formatter for chaining
        template<class Type, int nComp=0>
        formatter& openDataArray(const vtk::dataArrayAttr& attrEnum);


        //- Insert a single "PDataArray" XML entry tag.
        //  For some entries, the name is optional.
        //  \return formatter for chaining
        template<class Type, int nComp=0>
        formatter& PDataArray(const word& dataName);


        //- End "DataArray" XML tag
        //  \return formatter for chaining
        inline formatter& endDataArray();

        //- End "AppendedData" XML tag
        //  \return formatter for chaining
        formatter& endAppendedData();

        //- End "VTKFile" XML tag
        //  \return formatter for chaining
        formatter& endVTKFile();


        //- Write XML attribute
        //  \return formatter for chaining
        formatter& xmlAttr
        (
            const word& k,
            const std::string& v,
            const char quote = '\''
        );

        //- Write XML attribute
        //  \return formatter for chaining
        inline formatter& xmlAttr
        (
            const word& k,
            const int32_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        //  \return formatter for chaining
        inline formatter& xmlAttr
        (
            const word& k,
            const int64_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        //  \return formatter for chaining
        inline formatter& xmlAttr
        (
            const word& k,
            const uint64_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        //  \return formatter for chaining
        inline formatter& xmlAttr
        (
            const word& k,
            const scalar v,
            const char quote = '\''
        );

        //- Write XML attribute
        //  \return formatter for chaining
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& attrEnum,
            const int32_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        //  \return formatter for chaining
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& attrEnum,
            const int64_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        //  \return formatter for chaining
        inline formatter& xmlAttr
        (
            const vtk::fileAttr& attrEnum,
            const uint64_t v,
            const char quote = '\''
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkFormatterI.H"

#ifdef NoRepository
    #include "foamVtkFormatterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
