/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::expressions::exprDriver

Description
    Base driver for parsing (field) values.

    Largely based on code and ideas from swak4foam

    Properties
    \table
        Property     | Description                          | Required | Default
        variables    | List of variables for expressions    | no  | ()
        allowShadowing | Allow variables to shadow field names | no  | false
    \endtable

    Debug Properties
    \table
        Property     | Description                          | Required | Default
        debugBaseDriver | Debug level (int) for base driver | no  |
        debugScanner | Add debug for scanner                | no  | false
        debugParser  | Add debug for parser                 | no  | false
    \endtable

SourceFiles
    exprDriverI.H
    exprDriver.C
    exprDriverFields.C
    exprDriverIO.C
    exprDriverTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_exprDriver_H
#define expressions_exprDriver_H

#include "exprString.H"
#include "exprResult.H"
#include "pointField.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class TimeState;

namespace expressions
{

/*---------------------------------------------------------------------------*\
                         Class exprDriver Declaration
\*---------------------------------------------------------------------------*/

class exprDriver
{
public:

    // Data Types

        //- Search/caching controls
        enum searchControls
        {
            NO_SEARCH = 0,
            SEARCH_REGISTRY = 1,   //!< Search registry before disk
            SEARCH_FILES = 2,      //!< Search disk (eg, standalone app)
            CACHE_READ_FIELDS = 4, //!< Cache fields read from disk
            DEFAULT_SEARCH = (SEARCH_REGISTRY)
        };


private:

    // Private Member Functions

        //- Get search/caching controls from dictionary entries
        static int getSearchControls(const dictionary& dict);


protected:

    // Protected Data

    // Stored Data

        //- The dictionary with all input data/specification
        const dictionary& dict_;

        //- The result
        exprResult result_;

        //- Variable definitions, as read from a dictionary
        List<expressions::exprString> variableStrings_;

        //- The variables table
        HashTable<exprResult> variables_;

        //- Special-purpose scalar reference argument
        scalar arg1Value_;

        //- Reference to the time-state
        mutable const TimeState* timeStatePtr_;


    // Controls, tracing etc.

        //- Internal bookkeeping as "look-behind" parsing context
        mutable int stashedTokenId_;

        //- Request debugging for scanner
        bool debugScanner_;

        //- Request debugging for parser
        bool debugParser_;

        //- Allow variable names to mask field names
        bool allowShadowing_;

        //- Use value of previous iteration when oldTime is requested
        bool prevIterIsOldTime_;

        //- Registry/disk/caching control
        searchControls searchCtrl_;


    // Protected Member Functions

        inline bool searchRegistry() const noexcept;
        inline bool searchFiles() const noexcept;
        inline bool cacheReadFields() const noexcept;

        //- Reset the time-state reference
        void resetTimeReference(const TimeState* ts);

        //- Reset the time-state reference
        void resetTimeReference(const TimeState& ts);

        //- Lookup field object.
        //  \return const-ref tmp or invalid if not found
        template<class GeoField>
        static tmp<GeoField>
        cfindFieldObject(const objectRegistry& obr, const word& fldName);

        //- Read an interpolation table
        template<typename TableType>
        static bool readTable
        (
            const word& name,
            const dictionary& dict,
            HashTable<TableType>& tbl,
            bool clear=true
        );

        //- Write an interpolation table
        template<class TableType>
        static void writeTable
        (
            Ostream& os,
            const word& name,
            const HashTable<TableType>& tbl
        );


    // Variables

        //- Non-const access to the named variable (sub-classes only)
        inline virtual exprResult& variable(const word& name);


    // Fields

        //- Fill a random field
        //
        //  \param field the field to populate
        //  \param seed the seed value. If zero or negative, use as an offset
        //      to the current timeIndex (if a time-state is available)
        //  \param gaussian generate a Gaussian distribution
        void fill_random
        (
            scalarField& field,
            label seed = 0,
            const bool gaussian = false
        ) const;

        //- The (global) weighted average of a field, with stabilisation
        template<class Type>
        static Type weightedAverage
        (
            const scalarField& weights,
            const Field<Type>& fld
        );

        //- The (global) weighted sum (integral) of a field
        template<class Type>
        static Type weightedSum
        (
            const scalarField& weights,
            const Field<Type>& fld
        );

        //- Return the location of the min value
        static point getPositionOfMinimum
        (
            const scalarField& vals,
            const pointField& locs
        );

        //- Return the location of the max value
        static point getPositionOfMaximum
        (
            const scalarField& vals,
            const pointField& locs
        );


    // Updating

        //- Update things
        virtual bool update();

        //- Examine current variable values and update stored variables
        virtual void updateSpecialVariables(bool force=false);


    // Results

        //- Get the result from another driver.
        //  Override to allow mapping
        virtual exprResult getRemoteResult(const exprDriver& other) const;


        //- No copy assignment
        void operator=(const exprDriver&) = delete;


public:

    //- Runtime type information
    TypeName("exprDriver");


    // Constructors

        //- Default construct, and default construct with search preferences
        explicit exprDriver
        (
            enum searchControls search = searchControls::DEFAULT_SEARCH,
            const dictionary& dict = dictionary::null,
            const TimeState* ts = nullptr
        );

        //- Copy construct
        exprDriver(const exprDriver& rhs);

        //- Construct from a dictionary
        explicit exprDriver
        (
            const dictionary& dict,
            const TimeState* ts = nullptr
        );


    //- Destructor
    virtual ~exprDriver() = default;


    // Public Member Functions

        //- The underlying field size for the expression
        virtual label size() const
        {
            return 1;
        }

        //- The underlying point field size for the expression
        virtual label pointSize() const
        {
            return 1;
        }

        //- Reference to the current time-state (can be nullptr)
        const TimeState* timeState() const;

        //- The current time value
        virtual scalar timeValue() const;

        //- The current deltaT value
        virtual scalar deltaT() const;


        //- The dictionary with all input data/specification
        const dictionary& dict() const noexcept
        {
            return dict_;
        }

        //- Const access to expression result
        const exprResult& result() const noexcept
        {
            return result_;
        }

        //- Non-const access to expression result
        exprResult& result() noexcept
        {
            return result_;
        }

        //- Clear the result
        void clearResult();

        //- Return the expression result as a tmp field
        template<class Type>
        tmp<Field<Type>> getResult(bool wantPointData=false);

        //- The result type as word - same as result().valueType()
        const word& getResultType() const noexcept
        {
            return result_.valueType();
        }


    // Specials

        //- Get special-purpose scalar reference argument.
        //  Typically available as \c arg() in an expression and
        //  may correspond to table index, time value etc.
        inline scalar argValue() const noexcept;


    // General Controls

        //- Get "look-behind" parsing context (internal bookkeeping)
        inline int stashedTokenId() const noexcept;

        //- Reset "look-behind" parsing context (mutable operation)
        //  \return the previous value
        inline int resetStashedTokenId(int tokenId=0) const noexcept;


        //- Set the scanner/parser debug
        void setDebugging(bool scannerDebug, bool parserDebug);

        //- Set the scanner/parser debug to match the input
        void setDebugging(const exprDriver& rhs);

        //- Toggle CACHE_READ_FIELDS control
        bool setCaching(bool on) noexcept;

        //- Set search behaviour,
        //- with additional CACHE_READ_FIELDS toggle on
        void setSearchBehaviour
        (
            enum searchControls search,
            const bool caching = false
        );

        //- Set search behaviour to be identical to rhs
        void setSearchBehaviour(const exprDriver& rhs);

        //- Read access to scanner debug
        inline bool debugScanner() const noexcept;

        //- Read access to parser debug
        inline bool debugParser() const noexcept;

        bool prevIterIsOldTime() const { return prevIterIsOldTime_; }


    // Variables

        //- Clear temporary variables and resets from expression strings
        virtual void clearVariables();

        //- Set special-purpose scalar reference argument.
        //  Typically available as \c arg() in an expression and
        //  may corrspond to table index, time value etc.
        inline void setArgument(const scalar val) noexcept;

        //- True if named variable exists
        inline virtual bool hasVariable(const word& name) const;

        //- Return const-access to the named variable
        inline virtual const exprResult& variable(const word& name) const;

        //- Add/set string expressions for variables
        //  Can include multiple definitions inline
        void addVariables
        (
            const expressions::exprString& expr,
            bool clear = true  //!< Remove previously defined variables
        );

        //- Add/set string expressions for variables
        //  Can include multiple definitions inline
        void addVariables
        (
            const UList<expressions::exprString>& list,
            bool clear = true  //!< Remove previously defined variables
        );

        //- Add a uniform variable from an outside caller
        template<class T>
        inline void addUniformVariable
        (
            const word& name,
            const T& val
        );


    // Fields

        //- Test existence of a local variable
        template<class T>
        bool isLocalVariable
        (
            const word& name,
            bool wantPointData = false,
            label expectedSize = -1
        ) const;

        //- Retrieve local/global variable as a tmp field
        //
        //  \param name The name of the local/global field
        //  \param expectSize  The size check on the variable, -1 to ignore
        //  \param mandatory A missing variable is Fatal, or return
        //      an invalid tmp
        template<class Type>
        tmp<Field<Type>> getLocalVariable
        (
            const word& name,
            label expectSize,
            const bool mandatory = true
        ) const;


    // Evaluation

        //- Execute the parser.
        //  The return value currently has no meaning.
        virtual unsigned parse
        (
            const std::string& expr,
            size_t pos = 0,
            size_t len = std::string::npos
        ) = 0;

        //- Evaluate the expression and return the field
        template<class Type>
        inline tmp<Field<Type>>
        evaluate
        (
            const expressions::exprString& expr,
            bool wantPointData = false
        );

        //- Evaluate the expression and return a single value
        template<class Type>
        inline Type evaluateUniform
        (
            const expressions::exprString& expr,
            bool wantPointData = false
        );


        //- Evaluate the expression
        //- and save as the specified named variable
        void evaluateVariable
        (
            const word& varName,
            const expressions::exprString& expr
        );

        //- Evaluate an expression on a remote
        //- and save as the specified named variable
        virtual void evaluateVariableRemote
        (
            string remote,
            const word& varName,
            const expressions::exprString& expr
        );


    // Fields

        //- Return a new field with the size()
        template<class Type>
        tmp<Field<Type>>
        newField(const Type& val = pTraits<Type>::zero) const;

        //- Return a new field with the pointSize()
        template<class Type>
        tmp<Field<Type>>
        newPointField(const Type& val = pTraits<Type>::zero) const;


    // Reading

        //- Read an expression string and do substitutions
        static expressions::exprString readExpression
        (
            const word& name,
            const dictionary& dict
        );

        //- Read the list of variable strings
        //  (or initialize with a single string)
        static List<expressions::exprString> readVariableStrings
        (
            const dictionary& dict,
            const word& name = "variables",
            bool mandatory = false
        );

        //- Read an expression string (with the current dictionary)
        //- and do substitutions
        expressions::exprString readExpression(const word& name);


        //- Read variables, tables etc.
        //  Also usable for objects not constructed from a dictionary.
        virtual bool readDict(const dictionary& dict);

        //- Read "variables" and assigns to the list of expression strings
        //  \return the number variable strings read.
        label setVariableStrings
        (
            const dictionary& dict,
            bool mandatory = false
        );


    // Writing

        //- Write "variables"
        Ostream& writeVariableStrings
        (
            Ostream& os,
            const word& keyword = ""
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "exprDriverI.H"

#ifdef NoRepository
    #include "exprDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
