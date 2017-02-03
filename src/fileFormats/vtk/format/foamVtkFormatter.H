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
    foamVtkFormatter

Description
    Abstract class for a VTK output stream formatter.

    Includes very simple support for writing XML tags.

SourceFiles
    foamVtkFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkFormatter_H
#define foamVtkFormatter_H

#include "int.H"
#include "uint64.H"
#include "label.H"
#include "word.H"
#include "UList.H"
#include "LIFOStack.H"
#include "foamVtkPTraits.H"

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class foamVtkFormatter Declaration
\*---------------------------------------------------------------------------*/

class foamVtkFormatter
{
    // Private Data

        //- The output stream for the formatter
        std::ostream& os_;

        //- Stack of current XML tags
        LIFOStack<word> xmlTags_;

        //- Tag open/closed/ended state
        mutable bool inTag_;


        //- Write XML attribute
        template<class Type>
        foamVtkFormatter& xmlAttribute
        (
            const word& k,
            const Type& v,
            const char quote
        );


protected:

    // Protected Member Functions

        //- Construct and attach to an output stream
        foamVtkFormatter(std::ostream& os);

public:

    // Static Data

        //- VTK name for the 'byte_order' attribute
        static const char* const byteOrder;

        //- VTK name for the 'header_type' attribute (UInt64)
        static const char* const headerType;


    //- Destructor
    virtual ~foamVtkFormatter();


    // Member Functions

        //- Access to the underlying output stream
        inline std::ostream& os()
        {
            return os_;
        }


        //- Name for the XML output type.
        //  Possibly the lowercase version of the Legacy output type
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

        //- Indent according to the currently nested XML tags
        void indent();

        //- Write XML header
        foamVtkFormatter& xmlHeader();

        //- Write XML comment (at the current indentation level)
        foamVtkFormatter& comment(const std::string& text);


        //- Open XML tag
        foamVtkFormatter& openTag(const word& tag);

        //- Close XML tag, optional as an empty container.
        //  Always adds a trailing newline.
        foamVtkFormatter& closeTag(const bool isEmpty = false);

        //- End XML tag, optional with sanity check
        //  Always adds a trailing newline.
        foamVtkFormatter& endTag(const word& tag = word::null);

        //- Write XML tag without any attributes. Combines openTag/closeTag.
        foamVtkFormatter& tag(const word& tag);


        //- Open "DataArray" XML tag
        template<class Type, int nComp=0>
        foamVtkFormatter& openDataArray(const word& dataName);


        //- Insert a single "PDataArray" XML entry tag.
        //  For some entries, the name is optional.
        template<class Type, int nComp=0>
        foamVtkFormatter& PDataArray(const word& dataName);


        //- End "DataArray" XML tag
        foamVtkFormatter& endDataArray()
        {
            return endTag("DataArray");
        }


        //- Write XML attribute
        foamVtkFormatter& xmlAttr
        (
            const word& k,
            const std::string& v,
            const char quote = '\''
        );

        //- Write XML attribute
        foamVtkFormatter& xmlAttr
        (
            const word& k,
            const int32_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        foamVtkFormatter& xmlAttr
        (
            const word& k,
            const int64_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        foamVtkFormatter& xmlAttr
        (
            const word& k,
            const uint64_t v,
            const char quote = '\''
        );

        //- Write XML attribute
        foamVtkFormatter& xmlAttr
        (
            const word& k,
            const scalar v,
            const char quote = '\''
        );



    // Member Operators

        //- Write XML attribute
        foamVtkFormatter& operator()(const word& k, const std::string& v);

        //- Write XML attribute
        foamVtkFormatter& operator()(const word& k, const int32_t v);

        //- Write XML attribute
        foamVtkFormatter& operator()(const word& k, const int64_t v);

        //- Write XML attribute
        foamVtkFormatter& operator()(const word& k, const uint64_t v);

        //- Write XML attribute
        foamVtkFormatter& operator()(const word& k, const scalar v);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkFormatterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
