/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::IOstream

Description
    An IOstream is an abstract base class for all input/output systems; be
    they streams, files, token lists etc.

    The basic operations are construct, close, read token, read primitive
    and read binary block.  In addition version control and line number
    counting is incorporated.  Usually one would use the read primitive
    member functions, but if one were reading a stream on unknown data
    sequence one can read token by token, and then analyse.

SourceFiles
    IOstream.C

\*---------------------------------------------------------------------------*/

#ifndef IOstream_H
#define IOstream_H

#include "char.H"
#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include "scalar.H"
#include "fileName.H"
#include "InfoProxy.H"
#include "IOstreamOption.H"

#include <iostream>

using std::ios_base;
using std::istream;
using std::ostream;

using std::cin;
using std::cout;
using std::cerr;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOstream Declaration
\*---------------------------------------------------------------------------*/

class IOstream
:
    public IOstreamOption
{
public:

    // Public Data Types

        //- Enumeration for stream open/closed state
        enum streamAccess : char
        {
            CLOSED = 0,         //!< stream not open
            OPENED              //!< stream is open
        };


    // Public Static Data

        //- Default precision
        static unsigned int precision_;


protected:

    // Protected Data

        //- Name for any generic stream - normally treat as readonly
        static fileName staticName_;

        streamAccess openClosed_;

        ios_base::iostate ioState_;

        //- The label byte-size (could also be stored as byte)
        unsigned short labelByteSize_;

        //- The scalar byte-size (could also be stored as byte)
        unsigned short scalarByteSize_;

        //- The file line
        label lineNumber_;


    // Protected Member Functions

    // Access

        //- Set stream opened
        void setOpened()
        {
            openClosed_ = OPENED;
        }

        //- Set stream closed
        void setClosed()
        {
            openClosed_ = CLOSED;
        }

        //- Set stream state
        void setState(ios_base::iostate state)
        {
            ioState_ = state;
        }

        //- Set stream to be good
        void setGood()
        {
            ioState_ = ios_base::iostate(0);
        }


public:

    // Generated Methods

        //- Copy construct
        IOstream(const IOstream&) = default;

        //- Destructor
        virtual ~IOstream() = default;


    // Constructors

        //- Construct with specified stream option
        explicit IOstream(const IOstreamOption option)
        :
            IOstreamOption(option),
            openClosed_(CLOSED),
            ioState_(ios_base::iostate(0)),
            labelByteSize_(sizeof(label)),
            scalarByteSize_(sizeof(scalar)),
            lineNumber_(0)
        {
            setBad();
        }

        //- Construct with format, version
        IOstream
        (
            streamFormat format,
            versionNumber version,
            compressionType compression=UNCOMPRESSED
        )
        :
            IOstream(IOstreamOption(format, version, compression))
        {}


    // Member Functions

    // Access

        //- Return the name of the stream
        //  Useful for Fstream to return the filename
        virtual const fileName& name() const;

        //- Return non-const access to the name of the stream
        //  Useful to alter the stream name
        virtual fileName& name();


    // Check

        //- Check IOstream status for given operation.
        //  Print IOstream state if error has occurred
        virtual bool check(const char* operation) const;

        //- Check IOstream status for given operation.
        //  Print IOstream state if error has occurred and exit
        void fatalCheck(const char* operation) const;

        //- Return true if stream has been opened
        bool opened() const
        {
            return openClosed_ == OPENED;
        }

        //- Return true if stream is closed
        bool closed() const
        {
            return openClosed_ == CLOSED;
        }

        //- Return true if next operation might succeed
        bool good() const
        {
            return ioState_ == 0;
        }

        //- Return true if end of input seen
        bool eof() const
        {
            return ioState_ & ios_base::eofbit;
        }

        //- Return true if next operation will fail
        bool fail() const
        {
            return ioState_ & (ios_base::badbit | ios_base::failbit);
        }

        //- Return true if stream is corrupted
        bool bad() const
        {
            return ioState_ & ios_base::badbit;
        }

        //- Return true if the stream has not failed
        explicit operator bool() const
        {
            return !fail();
        }

        //- Return true if the stream has failed
        bool operator!() const
        {
            return fail();
        }


    // Element sizes (precision)

        //- The label byte-size associated with the stream
        unsigned labelByteSize() const
        {
            return labelByteSize_;
        }

        //- The scalar byte-size associated with the stream
        unsigned scalarByteSize() const
        {
            return scalarByteSize_;
        }

        //- Set the label byte-size associated with the stream
        void setLabelByteSize(unsigned nbytes)
        {
            labelByteSize_ = nbytes;
        }

        //- Set the scalar byte-size associated with the stream
        void setScalarByteSize(unsigned nbytes)
        {
            scalarByteSize_ = nbytes;
        }


        //- Check if the label byte-size associated with the stream
        //- is the same as the given type
        template<class T = label>
        typename std::enable_if<std::is_integral<T>::value, bool>::type
        checkLabelSize() const
        {
            return labelByteSize_ == sizeof(T);
        }

        //- Check if the scalar byte-size associated with the stream
        //- is the same as the given type
        template<class T = scalar>
        typename std::enable_if<std::is_floating_point<T>::value, bool>::type
        checkScalarSize() const
        {
            return scalarByteSize_ == sizeof(T);
        }


    // Stream State Functions

        //- Const access to the current stream line number
        label lineNumber() const
        {
            return lineNumber_;
        }

        //- Non-const access to the current stream line number
        label& lineNumber()
        {
            return lineNumber_;
        }

        //- Set the stream line number
        //  \return the previous value
        label lineNumber(const label num)
        {
            const label old(lineNumber_);
            lineNumber_ = num;
            return old;
        }

        //- Return flags of stream
        virtual ios_base::fmtflags flags() const = 0;

        //- Return the default precision
        static unsigned int defaultPrecision()
        {
            return precision_;
        }

        //- Reset the default precision
        //  \return the previous value
        static unsigned int defaultPrecision(unsigned int prec)
        {
            unsigned int old(precision_);
            precision_ = prec;
            return old;
        }

        //- Set stream to have reached eof
        void setEof()
        {
            ioState_ |= ios_base::eofbit;
        }

        //- Set stream to have failed
        void setFail()
        {
            ioState_ |= ios_base::failbit;
        }

        //- Set stream to be bad
        void setBad()
        {
            ioState_ |= ios_base::badbit;
        }

        //- Set flags of stream
        virtual ios_base::fmtflags flags(const ios_base::fmtflags f) = 0;

        //- Set flags of stream
        ios_base::fmtflags setf(const ios_base::fmtflags f)
        {
            return flags(flags() | f);
        }

        //- Set flags of given field of stream
        ios_base::fmtflags setf
        (
            const ios_base::fmtflags f,
            const ios_base::fmtflags mask
        )
        {
            return flags((flags() & ~mask) | (f & mask));
        }

        //- Unset flags of stream
        void unsetf(const ios_base::fmtflags f)
        {
            flags(flags() & ~f);
        }


    // Print

        //- Print description of IOstream to Ostream
        virtual void print(Ostream& os) const;

        //- Print information about the stream state bits
        void print(Ostream& os, const int streamState) const;


    // Info

        //- Return info proxy.
        //  Used to print IOstream information to a stream
        InfoProxy<IOstream> info() const
        {
            return *this;
        }
};


// Ostream Operator

template<>
Ostream& operator<<(Ostream& os, const InfoProxy<IOstream>& ip);


// --------------------------------------------------------------------
// ------ Manipulators (not taking arguments)
// --------------------------------------------------------------------

//- An IOstream manipulator
typedef IOstream& (*IOstreamManip)(IOstream&);

//- operator<< handling for manipulators without arguments
inline IOstream& operator<<(IOstream& io, IOstreamManip f)
{
    return f(io);
}


inline IOstream& dec(IOstream& io)
{
    io.setf(ios_base::dec, ios_base::dec|ios_base::hex|ios_base::oct);
    return io;
}

inline IOstream& hex(IOstream& io)
{
    io.setf(ios_base::hex, ios_base::dec|ios_base::hex|ios_base::oct);
    return io;
}

inline IOstream& oct(IOstream& io)
{
    io.setf(ios_base::oct, ios_base::dec|ios_base::hex|ios_base::oct);
    return io;
}

inline IOstream& fixed(IOstream& io)
{
    io.setf(ios_base::fixed, ios_base::floatfield);
    return io;
}

inline IOstream& scientific(IOstream& io)
{
    io.setf(ios_base::scientific, ios_base::floatfield);
    return io;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
