
#line 1 "volumeExprScanner.rl"
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Ragel lexer interface for lemon grammar of volume field expressions

\*---------------------------------------------------------------------------*/

#include "volumeExprScanner.H"
#include "volumeExprDriver.H"
#include "volumeExprLemonParser.h"
#include "volumeExprParser.H"
#include "Enum.H"
#include "macros.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wunused-const-variable"

// Debugging to stderr
#undef  DebugInfo
#define DebugInfo if (debug & 0x2) InfoErr


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

//- Paste token prefix
#define TOKEN_OF(T)         TOK_##T

//- An {int, c_str} enum pairing
#define TOKEN_PAIR(Name,T)  { TOKEN_OF(T), Name }

//- An {int, c_str} enum pairing for field types
#define FIELD_PAIR(Fld,T)  { TOKEN_OF(T), Fld::typeName.c_str() }


// Special handling for these known (stashed) look-back types
static const Enum<int> lookBehindTokenEnums
({
    TOKEN_PAIR("cellSet", CELL_SET),
    TOKEN_PAIR("faceSet", FACE_SET),
    TOKEN_PAIR("pointSet", POINT_SET),
    TOKEN_PAIR("cellZone", CELL_ZONE),
    TOKEN_PAIR("faceZone", FACE_ZONE),
    TOKEN_PAIR("pointZone", POINT_ZONE),
});

#define HAS_LOOKBEHIND_TOKENS


// Special handling of predefined method types. Eg, .x(), .y(), ...
static const Enum<int> fieldMethodEnums
({
    TOKEN_PAIR("x", CMPT_X),
    TOKEN_PAIR("y", CMPT_Y),
    TOKEN_PAIR("z", CMPT_Z),
    TOKEN_PAIR("xx", CMPT_XX),
    TOKEN_PAIR("xy", CMPT_XY),
    TOKEN_PAIR("xz", CMPT_XZ),
    TOKEN_PAIR("yx", CMPT_YX),
    TOKEN_PAIR("yy", CMPT_YY),
    TOKEN_PAIR("yz", CMPT_YZ),
    TOKEN_PAIR("zx", CMPT_ZX),
    TOKEN_PAIR("zy", CMPT_ZY),
    TOKEN_PAIR("zz", CMPT_ZZ),
    TOKEN_PAIR("ii", CMPT_II),
    TOKEN_PAIR("diag", DIAG),   /* tensors only */
    TOKEN_PAIR("T", TRANSPOSE), /* tensors only */
});


// Known field-token types
static const Enum<int>& fieldTokenEnums()
{
    static Enum<int> enums_;

    if (enums_.empty())
    {
        enums_.append
        ({
        #ifdef TOK_SCALAR_ID
            FIELD_PAIR(volScalarField, SCALAR_ID),
            FIELD_PAIR(volVectorField, VECTOR_ID),
            FIELD_PAIR(volTensorField, TENSOR_ID),
            FIELD_PAIR(volSymmTensorField, SYM_TENSOR_ID),
            FIELD_PAIR(volSphericalTensorField, SPH_TENSOR_ID),
        #else
            #error TOK_SCALAR_ID not defined
        #endif
        #ifdef TOK_SSCALAR_ID
            FIELD_PAIR(surfaceScalarField, SSCALAR_ID),
            FIELD_PAIR(surfaceVectorField, SVECTOR_ID),
            FIELD_PAIR(surfaceTensorField, STENSOR_ID),
            FIELD_PAIR(surfaceSymmTensorField, SSYM_TENSOR_ID),
            FIELD_PAIR(surfaceSphericalTensorField, SSPH_TENSOR_ID),
        #else
            #warning TOK_SSCALAR_ID not defined
        #endif
        #ifdef TOK_PSCALAR_ID
            FIELD_PAIR(pointScalarField, PSCALAR_ID),
            FIELD_PAIR(pointVectorField, PVECTOR_ID),
            FIELD_PAIR(pointTensorField, PTENSOR_ID),
            FIELD_PAIR(pointSymmTensorField, PSYM_TENSOR_ID),
            FIELD_PAIR(pointSphericalTensorField, PSPH_TENSOR_ID),
        #else
            #warning TOK_PSCALAR_ID not defined
        #endif
        });
    }

    return enums_;
}


// Simple compile-time function name declarations.
// Useful for handling driver-specific dispatching, or functions that
// are not universally available.
static const Enum<int> funcTokenEnums
({
#ifdef TOK_FLOOR
    TOKEN_PAIR("floor", FLOOR),
    TOKEN_PAIR("ceil", CEIL),
    TOKEN_PAIR("round", ROUND),
#endif
#ifdef TOK_HYPOT
    TOKEN_PAIR("hypot", HYPOT),
#endif

    // Already parsed as function: TOKEN_PAIR("pos", CELL_CENTRE),

    TOKEN_PAIR("point", POINT_EXPR), // Point value
    TOKEN_PAIR("face", FACE_EXPR),   // Face value, Face areaNormal

    TOKEN_PAIR("cellToPoint", CELL_TO_POINT),
    TOKEN_PAIR("cellToFace",  CELL_TO_FACE),
    TOKEN_PAIR("pointToCell", POINT_TO_CELL),

    TOKEN_PAIR("area", FACE_AREA),
    TOKEN_PAIR("vol",  CELL_VOLUME),

    TOKEN_PAIR("fpos", FACE_CENTRE),
    TOKEN_PAIR("pts", POINTS),
});

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Classifying token type based on an identifier name is indeed ugly.
//
// 1)
//   Handle special cases (eg, cellSet,...) first that have been tagged
//   as expected content with the stashed "look-behind" token.
//   Handle not-found errors here directly.
//
// 2)
//   Fallback to determining which field-type (volScalarField etc) the name
//   corresponds to.
//   Handle not-found errors by return -1.
//
static int driverTokenType
(
    const expressions::volumeExpr::parseDriver& driver_,
    const word& ident
)
{
    // Get stashed "look-behind" to decide what type of identifier we expect
    const int lookBehind = driver_.resetStashedTokenId();

    if (lookBehind && lookBehindTokenEnums.found(lookBehind))
    {
        bool good = false;

        switch (lookBehind)
        {
            case TOK_CELL_SET : good = driver_.isCellSet(ident); break;
            case TOK_FACE_SET : good = driver_.isFaceSet(ident); break;
            case TOK_POINT_SET : good = driver_.isPointSet(ident); break;
            case TOK_CELL_ZONE : good = driver_.isCellZone(ident); break;
            case TOK_FACE_ZONE : good = driver_.isFaceZone(ident); break;
            case TOK_POINT_ZONE : good = driver_.isPointZone(ident); break;
        }

        if (good)
        {
            return TOK_IDENTIFIER;
        }

        // Fatal
        driver_.reportFatal
        (
            "Error no " + lookBehindTokenEnums.get(lookBehind) + ": " + ident
        );

        return -2;  // Extra safety
    }

    // Surface variables - distinguish from volume by size
    #ifdef TOK_SSCALAR_ID
    {
        const label len = driver_.mesh().nInternalFaces();

        #undef checkFieldToken
        #define checkFieldToken(TokType, Type)                                \
        if (driver_.isVariable<Type>(ident, false, len))                      \
        {                                                                     \
            return TokType;                                                   \
        }

        checkFieldToken(TOK_SSCALAR_ID, scalar);
        checkFieldToken(TOK_SVECTOR_ID, vector);
        checkFieldToken(TOK_SSYM_TENSOR_ID, symmTensor);
        checkFieldToken(TOK_SSPH_TENSOR_ID, sphericalTensor);
        checkFieldToken(TOK_STENSOR_ID, tensor);
    }
    #endif

    // Point variables
    #ifdef TOK_PSCALAR_ID
    {
        #undef checkFieldToken
        #define checkFieldToken(TokType, Type)                                \
        if (driver_.isVariable<Type>(ident, true))                            \
        {                                                                     \
            return TokType;                                                   \
        }

        checkFieldToken(TOK_PSCALAR_ID, scalar);
        checkFieldToken(TOK_PVECTOR_ID, vector);
        checkFieldToken(TOK_PSPH_TENSOR_ID, sphericalTensor);
        checkFieldToken(TOK_PSYM_TENSOR_ID, symmTensor);
        checkFieldToken(TOK_PTENSOR_ID, tensor);
    }
    #endif

    // Volume variables
    #ifdef TOK_SCALAR_ID
    {
        #undef checkFieldToken
        #define checkFieldToken(TokType, Type)                                \
        if (driver_.isVariable<Type>(ident, false))                           \
        {                                                                     \
            return TokType;                                                   \
        }

        checkFieldToken(TOK_SCALAR_ID, scalar);
        checkFieldToken(TOK_VECTOR_ID, vector);
        checkFieldToken(TOK_SPH_TENSOR_ID, sphericalTensor);
        checkFieldToken(TOK_SYM_TENSOR_ID, symmTensor);
        checkFieldToken(TOK_TENSOR_ID, tensor);
    }
    #endif

    #undef checkFieldToken

    // Check registered fields and/or disk-files
    {
        const word fieldType(driver_.getFieldClassName(ident));

        int tokType = fieldTokenEnums().lookup(fieldType, -1);

        if (tokType > 0)
        {
            return tokType;
        }
    }

    return -1;
}

} // End anonymous namespace


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names

#define EMIT_TOKEN(T)                                                         \
    driver_.parsePosition() = (ts-buf);                                       \
    DebugInfo<< STRINGIFY(T) << " at " << driver_.parsePosition() << nl;      \
    parser_->parse(TOKEN_OF(T), nullptr);                                     \
    driver_.parsePosition() = (p-buf);



#line 320 "volumeExprScanner.cc"
static const int volumeExpr_start = 11;
static const int volumeExpr_first_final = 11;
static const int volumeExpr_error = 0;

static const int volumeExpr_en_main = 11;


#line 460 "volumeExprScanner.rl"



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expressions::volumeExpr::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::expressions::volumeExpr::scanner::dispatch_method
(
    const parseDriver& driver_,
    scanToken& scanTok,
    word&& ident
) const
{
    if (ident[0] == '.')
    {
        ident.erase(0, 1);
    }

    DebugInfo
        << "Method:" << ident
        << " at " << driver_.parsePosition() << nl;

    const int methType = fieldMethodEnums.lookup(ident, -1);

    if (methType > 0)
    {
        // Dispatch '.' and "method" separately
        parser_->parse(TOK_DOT, nullptr);
        parser_->parse(methType, nullptr);

        return true;
    }

    driver_.reportFatal("Unknown method: " + ident);
    return false;
}


bool Foam::expressions::volumeExpr::scanner::dispatch_ident
(
    const parseDriver& driver_,
    scanToken& scanTok,
    word&& ident
) const
{
    int tokType = -1;

    const bool quoted =
    (
        (ident.front() == '"' || ident.front() == '\'')
     && (ident.front() == ident.back())
    );

    if (quoted)
    {
        ident.erase(ident.size()-1);
        ident.erase(0, 1);
    }
    else
    {
        // Check for function name
        tokType = funcTokenEnums.lookup(ident, -1);

        if (tokType > 0)
        {
            DebugInfo
                << "Emit:" << ident << " function:"
                << parser_->tokenName(tokType) << nl;

            parser_->parse(tokType, nullptr);
            return true;
        }

        #ifdef HAS_LOOKBEHIND_TOKENS
        // Specials such "cset" also reset the look-behind
        tokType = lookBehindTokenEnums.lookup(ident, -1);

        if (tokType > 0)
        {
            DebugInfo
                << "Emit:" << ident << " as look-behind:"
                << parser_->tokenName(tokType) << nl;

            driver_.resetStashedTokenId(tokType);
            parser_->parse(tokType, nullptr);
            return true;
        }
        #endif
    }


    // Can also peek at stashed "look-behind"
    // const int lookBehind = driver_.stashedTokenId();

    tokType = driverTokenType(driver_, ident);

    if (tokType > 0)
    {
        DebugInfo
            << "Emit:" << ident << " token:"
            << parser_->tokenName(tokType) << nl;

        scanTok.name = new Foam::word(std::move(ident));
        parser_->parse(tokType, &scanTok);

        return true;
    }


    // Not found? Attempt to strip off '.x' endings etc,
    // but not when quoted

    const auto dot = ident.rfind('.');
    const int methType =
    (
        quoted || dot == std::string::npos
      ? -1
      : fieldMethodEnums.lookup(ident.substr(dot+1), -1)
    );

    if
    (
        methType > 0
     && (tokType = driverTokenType(driver_, ident.substr(0, dot))) > 0
    )
    {
        DebugInfo
            << "Emit:" << ident.substr(0, dot).c_str() << " token:"
            << parser_->tokenName(tokType) << " with "
            << ident.substr(dot).c_str() << " token:"
            << parser_->tokenName(methType) << nl;

        // The field (before the ".")
        ident.erase(dot);

        scanTok.name = new Foam::word(std::move(ident));
        parser_->parse(tokType, &scanTok);

        // Dispatch '.' and "method" separately
        parser_->parse(TOK_DOT, nullptr);
        parser_->parse(methType, nullptr);

        return true;
    }

    driver_.reportFatal
    (
        "Object " + ident + " does not exist or wrong type"
    );

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::volumeExpr::scanner::process
(
    const std::string& str,
    size_t strBeg,
    size_t strLen,
    parseDriver& driver_
)
{
    // Save debug value
    const int oldDebug = debug;

    if (driver_.debugScanner()) { debug |= 0x2; }
    if (driver_.debugParser())  { debug |= 0x4; }

    if (debug & 0x6)
    {
        InfoErr
            << "Begin parse {"
            << str.substr(strBeg, strLen).c_str() << '}' << nl;
    }

    if (!parser_)
    {
        parser_ = new parser();
    }

    driver_.content(str, strBeg, strLen);

    size_t strEnd = str.length();

    if (strBeg > str.length())
    {
        strBeg = str.length();
    }
    else if (strLen != std::string::npos)
    {
        strLen += strBeg;

        if (strLen < str.length())
        {
            strEnd = strLen;
        }
    }


    parser_->start(driver_);

    // Scan token type
    scanToken scanTok;

    // Token start/end (Ragel naming)
    const char* ts;
    const char* te;

    // Local buffer data.
    // - p, pe, eof are Ragel naming
    // - buf is our own naming

    const char* buf = &(str[strBeg]);
    const char* eof = &(str[strEnd]);
    const char* p = buf;
    const char* pe = eof;

    // Initialize FSM variables
    
#line 560 "volumeExprScanner.cc"
	{
	cs = volumeExpr_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 690 "volumeExprScanner.rl"
   /* ^^^ FSM initialization here ^^^ */;

    
#line 572 "volumeExprScanner.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr2:
#line 342 "volumeExprScanner.rl"
	{te = p+1;{
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr4:
#line 342 "volumeExprScanner.rl"
	{te = p+1;{
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr5:
#line 320 "volumeExprScanner.rl"
	{{p = ((te))-1;}{
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        if (readScalar(std::string(ts, te-ts), scanTok.svalue))
        {
            parser_->parse(TOKEN_OF(NUMBER), &scanTok);
        }
        else
        {
            driver_.reportFatal
            (
                "Error parsing number: " + std::string(ts, te-ts)
            );
        }

        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr8:
#line 385 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(EQUAL); }}
	goto st11;
tr9:
#line 439 "volumeExprScanner.rl"
	{{p = ((te))-1;}{ EMIT_TOKEN(TENSOR); }}
	goto st11;
tr11:
#line 447 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(IDENTITY_TENSOR); }}
	goto st11;
tr12:
#line 388 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LOR); }}
	goto st11;
tr16:
#line 370 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PERCENT); }}
	goto st11;
tr19:
#line 371 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LPAREN); }}
	goto st11;
tr20:
#line 372 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(RPAREN); }}
	goto st11;
tr21:
#line 373 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(TIMES); }}
	goto st11;
tr22:
#line 374 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PLUS); }}
	goto st11;
tr23:
#line 376 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COMMA); }}
	goto st11;
tr24:
#line 375 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(MINUS); }}
	goto st11;
tr26:
#line 378 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(DIVIDE); }}
	goto st11;
tr28:
#line 380 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COLON); }}
	goto st11;
tr32:
#line 379 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(QUESTION); }}
	goto st11;
tr35:
#line 391 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(BIT_XOR); }}
	goto st11;
tr52:
#line 364 "volumeExprScanner.rl"
	{te = p;p--;}
	goto st11;
tr53:
#line 369 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(NOT); }}
	goto st11;
tr54:
#line 386 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(NOT_EQUAL); }}
	goto st11;
tr55:
#line 389 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(BIT_AND); }}
	goto st11;
tr56:
#line 387 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LAND); }}
	goto st11;
tr57:
#line 377 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(DOT); }}
	goto st11;
tr60:
#line 320 "volumeExprScanner.rl"
	{te = p;p--;{
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        if (readScalar(std::string(ts, te-ts), scanTok.svalue))
        {
            parser_->parse(TOKEN_OF(NUMBER), &scanTok);
        }
        else
        {
            driver_.reportFatal
            (
                "Error parsing number: " + std::string(ts, te-ts)
            );
        }

        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr62:
#line 348 "volumeExprScanner.rl"
	{te = p;p--;{
        // Tokenized ".method" - dispatch '.' and "method" separately
        driver_.parsePosition() = (ts-buf);
        dispatch_method(driver_, scanTok, word(ts+1, te-ts-1, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr63:
#line 381 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LESS); }}
	goto st11;
tr64:
#line 382 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LESS_EQ); }}
	goto st11;
tr65:
#line 383 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(GREATER); }}
	goto st11;
tr66:
#line 384 "volumeExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(GREATER_EQ); }}
	goto st11;
tr67:
#line 342 "volumeExprScanner.rl"
	{te = p;p--;{
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr69:
#line 1 "NONE"
	{	switch( act ) {
	case 26:
	{{p = ((te))-1;} EMIT_TOKEN(PI); }
	break;
	case 27:
	{{p = ((te))-1;} EMIT_TOKEN(DEG_TO_RAD); }
	break;
	case 28:
	{{p = ((te))-1;} EMIT_TOKEN(RAD_TO_DEG); }
	break;
	case 29:
	{{p = ((te))-1;} EMIT_TOKEN(EXP); }
	break;
	case 31:
	{{p = ((te))-1;} EMIT_TOKEN(LOG10); }
	break;
	case 32:
	{{p = ((te))-1;} EMIT_TOKEN(POW); }
	break;
	case 34:
	{{p = ((te))-1;} EMIT_TOKEN(SQRT); }
	break;
	case 35:
	{{p = ((te))-1;} EMIT_TOKEN(CBRT); }
	break;
	case 39:
	{{p = ((te))-1;} EMIT_TOKEN(ASIN); }
	break;
	case 40:
	{{p = ((te))-1;} EMIT_TOKEN(ACOS); }
	break;
	case 42:
	{{p = ((te))-1;} EMIT_TOKEN(ATAN2); }
	break;
	case 43:
	{{p = ((te))-1;} EMIT_TOKEN(SINH); }
	break;
	case 44:
	{{p = ((te))-1;} EMIT_TOKEN(COSH); }
	break;
	case 45:
	{{p = ((te))-1;} EMIT_TOKEN(TANH); }
	break;
	case 47:
	{{p = ((te))-1;} EMIT_TOKEN(MAGSQR); }
	break;
	case 50:
	{{p = ((te))-1;} EMIT_TOKEN(POS0); }
	break;
	case 51:
	{{p = ((te))-1;} EMIT_TOKEN(NEG0); }
	break;
	case 52:
	{{p = ((te))-1;} EMIT_TOKEN(SIGN); }
	break;
	case 53:
	{{p = ((te))-1;} EMIT_TOKEN(MIN); }
	break;
	case 54:
	{{p = ((te))-1;} EMIT_TOKEN(MAX); }
	break;
	case 55:
	{{p = ((te))-1;} EMIT_TOKEN(AVERAGE); }
	break;
	case 56:
	{{p = ((te))-1;} EMIT_TOKEN(SUM); }
	break;
	case 57:
	{{p = ((te))-1;} EMIT_TOKEN(WEIGHT_AVERAGE); }
	break;
	case 58:
	{{p = ((te))-1;} EMIT_TOKEN(WEIGHT_SUM); }
	break;
	case 59:
	{{p = ((te))-1;} EMIT_TOKEN(RAND); }
	break;
	case 60:
	{{p = ((te))-1;} EMIT_TOKEN(BOOL); }
	break;
	case 61:
	{{p = ((te))-1;} EMIT_TOKEN(VECTOR); }
	break;
	case 63:
	{{p = ((te))-1;} EMIT_TOKEN(SYM_TENSOR); }
	break;
	case 64:
	{{p = ((te))-1;} EMIT_TOKEN(SPH_TENSOR); }
	break;
	case 65:
	{{p = ((te))-1;} EMIT_TOKEN(ZERO); }
	break;
	case 66:
	{{p = ((te))-1;} EMIT_TOKEN(LTRUE); }
	break;
	case 67:
	{{p = ((te))-1;} EMIT_TOKEN(LFALSE); }
	break;
	case 69:
	{{p = ((te))-1;} EMIT_TOKEN(ARG); }
	break;
	case 70:
	{{p = ((te))-1;} EMIT_TOKEN(TIME); }
	break;
	case 71:
	{{p = ((te))-1;} EMIT_TOKEN(DELTA_T); }
	break;
	case 72:
	{{p = ((te))-1;}
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }
	break;
	}
	}
	goto st11;
tr85:
#line 413 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(ATAN); }}
	goto st11;
tr100:
#line 409 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(COS); }}
	goto st11;
tr121:
#line 402 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LOG); }}
	goto st11;
tr128:
#line 418 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(MAG); }}
	goto st11;
tr135:
#line 422 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(NEG); }}
	goto st11;
tr141:
#line 421 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(POS); }}
	goto st11;
tr160:
#line 408 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SIN); }}
	goto st11;
tr176:
#line 405 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SQR); }}
	goto st11;
tr192:
#line 410 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TAN); }}
	goto st11;
tr198:
#line 439 "volumeExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TENSOR); }}
	goto st11;
st11:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 1 "NONE"
	{ts = p;}
#line 925 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 32: goto st12;
		case 33: goto st13;
		case 34: goto st1;
		case 37: goto tr16;
		case 38: goto st14;
		case 39: goto st3;
		case 40: goto tr19;
		case 41: goto tr20;
		case 42: goto tr21;
		case 43: goto tr22;
		case 44: goto tr23;
		case 45: goto tr24;
		case 46: goto st15;
		case 47: goto tr26;
		case 58: goto tr28;
		case 60: goto st20;
		case 61: goto st7;
		case 62: goto st21;
		case 63: goto tr32;
		case 90: goto st24;
		case 94: goto tr35;
		case 95: goto st22;
		case 97: goto st27;
		case 98: goto st41;
		case 99: goto st44;
		case 100: goto st49;
		case 101: goto st59;
		case 102: goto st61;
		case 108: goto st65;
		case 109: goto st69;
		case 110: goto st75;
		case 112: goto st78;
		case 114: goto st81;
		case 115: goto st89;
		case 116: goto st117;
		case 118: goto st129;
		case 119: goto st134;
		case 124: goto st10;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st12;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 89 ) {
			if ( 103 <= (*p) && (*p) <= 122 )
				goto st22;
		} else if ( (*p) >= 65 )
			goto st22;
	} else
		goto tr27;
	goto st0;
st0:
cs = 0;
	goto _out;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 32 )
		goto st12;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st12;
	goto tr52;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 61 )
		goto tr54;
	goto tr53;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( (*p) == 34 )
		goto st0;
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( (*p) == 34 )
		goto tr2;
	goto st2;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 38 )
		goto tr56;
	goto tr55;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( (*p) == 39 )
		goto st0;
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( (*p) == 39 )
		goto tr4;
	goto st4;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr58;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st18;
	} else
		goto st18;
	goto tr57;
tr58:
#line 1 "NONE"
	{te = p+1;}
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 1053 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr58;
	goto tr60;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 43: goto st6;
		case 45: goto st6;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st17;
	goto tr5;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st17;
	goto tr5;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st17;
	goto tr60;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st18;
	} else if ( (*p) >= 65 )
		goto st18;
	goto tr62;
tr27:
#line 1 "NONE"
	{te = p+1;}
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 1104 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr58;
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr27;
	goto tr60;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 61 )
		goto tr64;
	goto tr63;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 61 )
		goto tr8;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 61 )
		goto tr66;
	goto tr65;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
tr68:
#line 1 "NONE"
	{te = p+1;}
#line 342 "volumeExprScanner.rl"
	{act = 72;}
	goto st23;
tr72:
#line 1 "NONE"
	{te = p+1;}
#line 444 "volumeExprScanner.rl"
	{act = 65;}
	goto st23;
tr79:
#line 1 "NONE"
	{te = p+1;}
#line 412 "volumeExprScanner.rl"
	{act = 40;}
	goto st23;
tr80:
#line 1 "NONE"
	{te = p+1;}
#line 448 "volumeExprScanner.rl"
	{act = 69;}
	goto st23;
tr82:
#line 1 "NONE"
	{te = p+1;}
#line 411 "volumeExprScanner.rl"
	{act = 39;}
	goto st23;
tr86:
#line 1 "NONE"
	{te = p+1;}
#line 414 "volumeExprScanner.rl"
	{act = 42;}
	goto st23;
tr91:
#line 1 "NONE"
	{te = p+1;}
#line 430 "volumeExprScanner.rl"
	{act = 55;}
	goto st23;
tr94:
#line 1 "NONE"
	{te = p+1;}
#line 437 "volumeExprScanner.rl"
	{act = 60;}
	goto st23;
tr98:
#line 1 "NONE"
	{te = p+1;}
#line 407 "volumeExprScanner.rl"
	{act = 35;}
	goto st23;
tr101:
#line 1 "NONE"
	{te = p+1;}
#line 416 "volumeExprScanner.rl"
	{act = 44;}
	goto st23;
tr109:
#line 1 "NONE"
	{te = p+1;}
#line 399 "volumeExprScanner.rl"
	{act = 27;}
	goto st23;
tr112:
#line 1 "NONE"
	{te = p+1;}
#line 450 "volumeExprScanner.rl"
	{act = 71;}
	goto st23;
tr114:
#line 1 "NONE"
	{te = p+1;}
#line 401 "volumeExprScanner.rl"
	{act = 29;}
	goto st23;
tr118:
#line 1 "NONE"
	{te = p+1;}
#line 446 "volumeExprScanner.rl"
	{act = 67;}
	goto st23;
tr123:
#line 1 "NONE"
	{te = p+1;}
#line 403 "volumeExprScanner.rl"
	{act = 31;}
	goto st23;
tr127:
#line 1 "NONE"
	{te = p+1;}
#line 429 "volumeExprScanner.rl"
	{act = 54;}
	goto st23;
tr131:
#line 1 "NONE"
	{te = p+1;}
#line 419 "volumeExprScanner.rl"
	{act = 47;}
	goto st23;
tr132:
#line 1 "NONE"
	{te = p+1;}
#line 428 "volumeExprScanner.rl"
	{act = 53;}
	goto st23;
tr136:
#line 1 "NONE"
	{te = p+1;}
#line 424 "volumeExprScanner.rl"
	{act = 51;}
	goto st23;
tr137:
#line 1 "NONE"
	{te = p+1;}
#line 398 "volumeExprScanner.rl"
	{act = 26;}
	goto st23;
tr140:
#line 1 "NONE"
	{te = p+1;}
#line 404 "volumeExprScanner.rl"
	{act = 32;}
	goto st23;
tr142:
#line 1 "NONE"
	{te = p+1;}
#line 423 "volumeExprScanner.rl"
	{act = 50;}
	goto st23;
tr150:
#line 1 "NONE"
	{te = p+1;}
#line 400 "volumeExprScanner.rl"
	{act = 28;}
	goto st23;
tr151:
#line 1 "NONE"
	{te = p+1;}
#line 434 "volumeExprScanner.rl"
	{act = 59;}
	goto st23;
tr159:
#line 1 "NONE"
	{te = p+1;}
#line 425 "volumeExprScanner.rl"
	{act = 52;}
	goto st23;
tr161:
#line 1 "NONE"
	{te = p+1;}
#line 415 "volumeExprScanner.rl"
	{act = 43;}
	goto st23;
tr174:
#line 1 "NONE"
	{te = p+1;}
#line 441 "volumeExprScanner.rl"
	{act = 64;}
	goto st23;
tr177:
#line 1 "NONE"
	{te = p+1;}
#line 406 "volumeExprScanner.rl"
	{act = 34;}
	goto st23;
tr178:
#line 1 "NONE"
	{te = p+1;}
#line 431 "volumeExprScanner.rl"
	{act = 56;}
	goto st23;
tr186:
#line 1 "NONE"
	{te = p+1;}
#line 440 "volumeExprScanner.rl"
	{act = 63;}
	goto st23;
tr193:
#line 1 "NONE"
	{te = p+1;}
#line 417 "volumeExprScanner.rl"
	{act = 45;}
	goto st23;
tr201:
#line 1 "NONE"
	{te = p+1;}
#line 449 "volumeExprScanner.rl"
	{act = 70;}
	goto st23;
tr203:
#line 1 "NONE"
	{te = p+1;}
#line 445 "volumeExprScanner.rl"
	{act = 66;}
	goto st23;
tr208:
#line 1 "NONE"
	{te = p+1;}
#line 438 "volumeExprScanner.rl"
	{act = 61;}
	goto st23;
tr221:
#line 1 "NONE"
	{te = p+1;}
#line 432 "volumeExprScanner.rl"
	{act = 57;}
	goto st23;
tr223:
#line 1 "NONE"
	{te = p+1;}
#line 433 "volumeExprScanner.rl"
	{act = 58;}
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 1371 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr69;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 99: goto st28;
		case 114: goto st30;
		case 115: goto st31;
		case 116: goto st33;
		case 118: goto st36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 115: goto tr79;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto tr80;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 105: goto st32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto tr82;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st34;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto st35;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 46: goto tr68;
		case 50: goto tr86;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr85;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto st38;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st39;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st40;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto tr91;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st42;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st43;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 108: goto tr94;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 98: goto st45;
		case 111: goto st47;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto st46;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 116: goto tr98;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 115: goto st48;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 104: goto tr101;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr100;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st51;
		case 108: goto st56;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 46: goto tr68;
		case 84: goto st52;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 46: goto tr68;
		case 82: goto st54;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st55;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 100: goto tr109;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 116: goto st57;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st58;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 46: goto tr68;
		case 84: goto tr112;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 120: goto st60;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 112: goto tr114;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st62;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 108: goto st63;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 115: goto st64;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto tr118;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st66;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	switch( (*p) ) {
		case 46: goto tr68;
		case 49: goto st68;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr121;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 46: goto tr68;
		case 48: goto tr123;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st70;
		case 105: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st71;
		case 120: goto tr127;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 46: goto tr68;
		case 83: goto st72;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr128;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 113: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto tr131;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto tr132;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st76;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st77;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	switch( (*p) ) {
		case 46: goto tr68;
		case 48: goto tr136;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr135;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 105: goto tr137;
		case 111: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 115: goto st80;
		case 119: goto tr140;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 46: goto tr68;
		case 48: goto tr142;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr141;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st82;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 100: goto st83;
		case 110: goto st88;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	switch( (*p) ) {
		case 46: goto tr68;
		case 84: goto st84;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 46: goto tr68;
		case 68: goto st86;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st87;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto tr150;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 100: goto tr151;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 105: goto st90;
		case 112: goto st93;
		case 113: goto st106;
		case 117: goto st108;
		case 121: goto st109;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st91;
		case 110: goto st92;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto tr159;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 104: goto tr161;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr160;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 104: goto st94;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st95;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto st96;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 105: goto st97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 99: goto st98;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st99;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 108: goto st100;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
	switch( (*p) ) {
		case 46: goto tr68;
		case 84: goto st101;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st102;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto st103;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 115: goto st104;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st105;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto tr174;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto st107;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 116: goto tr177;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr176;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 109: goto tr178;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 109: goto st110;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 109: goto st111;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
	switch( (*p) ) {
		case 46: goto tr68;
		case 84: goto st112;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st113;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto st114;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 115: goto st115;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st116;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto tr186;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st118;
		case 101: goto st120;
		case 105: goto st125;
		case 114: goto st127;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto st119;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 104: goto tr193;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr192;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 110: goto st121;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 115: goto st122;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st123;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto tr197;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
tr197:
#line 1 "NONE"
	{te = p+1;}
	goto st124;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
#line 3212 "volumeExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr68;
		case 58: goto st8;
		case 95: goto tr68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr198;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 58 )
		goto st9;
	goto tr9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 73 )
		goto tr11;
	goto tr9;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 109: goto st126;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto tr201;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 117: goto st128;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto tr203;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st129:
	if ( ++p == pe )
		goto _test_eof129;
case 129:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st130;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st130:
	if ( ++p == pe )
		goto _test_eof130;
case 130:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 99: goto st131;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st131:
	if ( ++p == pe )
		goto _test_eof131;
case 131:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 116: goto st132;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st132:
	if ( ++p == pe )
		goto _test_eof132;
case 132:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 111: goto st133;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st133:
	if ( ++p == pe )
		goto _test_eof133;
case 133:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto tr208;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st134:
	if ( ++p == pe )
		goto _test_eof134;
case 134:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st135;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st135:
	if ( ++p == pe )
		goto _test_eof135;
case 135:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 105: goto st136;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st136:
	if ( ++p == pe )
		goto _test_eof136;
case 136:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st137;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st137:
	if ( ++p == pe )
		goto _test_eof137;
case 137:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 104: goto st138;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st138:
	if ( ++p == pe )
		goto _test_eof138;
case 138:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 116: goto st139;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st139:
	if ( ++p == pe )
		goto _test_eof139;
case 139:
	switch( (*p) ) {
		case 46: goto tr68;
		case 65: goto st140;
		case 83: goto st146;
		case 95: goto tr68;
	}
	if ( (*p) < 66 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st140:
	if ( ++p == pe )
		goto _test_eof140;
case 140:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 118: goto st141;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st141:
	if ( ++p == pe )
		goto _test_eof141;
case 141:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto st142;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st142:
	if ( ++p == pe )
		goto _test_eof142;
case 142:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 114: goto st143;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st143:
	if ( ++p == pe )
		goto _test_eof143;
case 143:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 97: goto st144;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st144:
	if ( ++p == pe )
		goto _test_eof144;
case 144:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 103: goto st145;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st145:
	if ( ++p == pe )
		goto _test_eof145;
case 145:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 101: goto tr221;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st146:
	if ( ++p == pe )
		goto _test_eof146;
case 146:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 117: goto st147;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st147:
	if ( ++p == pe )
		goto _test_eof147;
case 147:
	switch( (*p) ) {
		case 46: goto tr68;
		case 95: goto tr68;
		case 109: goto tr223;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr68;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr68;
	} else
		goto tr68;
	goto tr67;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 124 )
		goto tr12;
	goto st0;
	}
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof129: cs = 129; goto _test_eof; 
	_test_eof130: cs = 130; goto _test_eof; 
	_test_eof131: cs = 131; goto _test_eof; 
	_test_eof132: cs = 132; goto _test_eof; 
	_test_eof133: cs = 133; goto _test_eof; 
	_test_eof134: cs = 134; goto _test_eof; 
	_test_eof135: cs = 135; goto _test_eof; 
	_test_eof136: cs = 136; goto _test_eof; 
	_test_eof137: cs = 137; goto _test_eof; 
	_test_eof138: cs = 138; goto _test_eof; 
	_test_eof139: cs = 139; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
	_test_eof142: cs = 142; goto _test_eof; 
	_test_eof143: cs = 143; goto _test_eof; 
	_test_eof144: cs = 144; goto _test_eof; 
	_test_eof145: cs = 145; goto _test_eof; 
	_test_eof146: cs = 146; goto _test_eof; 
	_test_eof147: cs = 147; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 12: goto tr52;
	case 13: goto tr53;
	case 14: goto tr55;
	case 15: goto tr57;
	case 16: goto tr60;
	case 5: goto tr5;
	case 6: goto tr5;
	case 17: goto tr60;
	case 18: goto tr62;
	case 19: goto tr60;
	case 20: goto tr63;
	case 21: goto tr65;
	case 22: goto tr67;
	case 23: goto tr69;
	case 24: goto tr67;
	case 25: goto tr67;
	case 26: goto tr67;
	case 27: goto tr67;
	case 28: goto tr67;
	case 29: goto tr67;
	case 30: goto tr67;
	case 31: goto tr67;
	case 32: goto tr67;
	case 33: goto tr67;
	case 34: goto tr67;
	case 35: goto tr85;
	case 36: goto tr67;
	case 37: goto tr67;
	case 38: goto tr67;
	case 39: goto tr67;
	case 40: goto tr67;
	case 41: goto tr67;
	case 42: goto tr67;
	case 43: goto tr67;
	case 44: goto tr67;
	case 45: goto tr67;
	case 46: goto tr67;
	case 47: goto tr67;
	case 48: goto tr100;
	case 49: goto tr67;
	case 50: goto tr67;
	case 51: goto tr67;
	case 52: goto tr67;
	case 53: goto tr67;
	case 54: goto tr67;
	case 55: goto tr67;
	case 56: goto tr67;
	case 57: goto tr67;
	case 58: goto tr67;
	case 59: goto tr67;
	case 60: goto tr67;
	case 61: goto tr67;
	case 62: goto tr67;
	case 63: goto tr67;
	case 64: goto tr67;
	case 65: goto tr67;
	case 66: goto tr67;
	case 67: goto tr121;
	case 68: goto tr67;
	case 69: goto tr67;
	case 70: goto tr67;
	case 71: goto tr128;
	case 72: goto tr67;
	case 73: goto tr67;
	case 74: goto tr67;
	case 75: goto tr67;
	case 76: goto tr67;
	case 77: goto tr135;
	case 78: goto tr67;
	case 79: goto tr67;
	case 80: goto tr141;
	case 81: goto tr67;
	case 82: goto tr67;
	case 83: goto tr67;
	case 84: goto tr67;
	case 85: goto tr67;
	case 86: goto tr67;
	case 87: goto tr67;
	case 88: goto tr67;
	case 89: goto tr67;
	case 90: goto tr67;
	case 91: goto tr67;
	case 92: goto tr160;
	case 93: goto tr67;
	case 94: goto tr67;
	case 95: goto tr67;
	case 96: goto tr67;
	case 97: goto tr67;
	case 98: goto tr67;
	case 99: goto tr67;
	case 100: goto tr67;
	case 101: goto tr67;
	case 102: goto tr67;
	case 103: goto tr67;
	case 104: goto tr67;
	case 105: goto tr67;
	case 106: goto tr67;
	case 107: goto tr176;
	case 108: goto tr67;
	case 109: goto tr67;
	case 110: goto tr67;
	case 111: goto tr67;
	case 112: goto tr67;
	case 113: goto tr67;
	case 114: goto tr67;
	case 115: goto tr67;
	case 116: goto tr67;
	case 117: goto tr67;
	case 118: goto tr67;
	case 119: goto tr192;
	case 120: goto tr67;
	case 121: goto tr67;
	case 122: goto tr67;
	case 123: goto tr67;
	case 124: goto tr198;
	case 8: goto tr9;
	case 9: goto tr9;
	case 125: goto tr67;
	case 126: goto tr67;
	case 127: goto tr67;
	case 128: goto tr67;
	case 129: goto tr67;
	case 130: goto tr67;
	case 131: goto tr67;
	case 132: goto tr67;
	case 133: goto tr67;
	case 134: goto tr67;
	case 135: goto tr67;
	case 136: goto tr67;
	case 137: goto tr67;
	case 138: goto tr67;
	case 139: goto tr67;
	case 140: goto tr67;
	case 141: goto tr67;
	case 142: goto tr67;
	case 143: goto tr67;
	case 144: goto tr67;
	case 145: goto tr67;
	case 146: goto tr67;
	case 147: goto tr67;
	}
	}

	_out: {}
	}

#line 692 "volumeExprScanner.rl"
   /* ^^^ FSM execution here ^^^ */;

    if (0 == cs)
    {
        driver_.reportFatal("Parse error while scanning", (p-buf));
    }

    if (p != eof)
    {
        driver_.reportFatal("Parsing failed with remaining content", (p-buf));
    }

    // Terminate parser execution
    parser_->parse(0, nullptr);
    parser_->stop();

    if (debug & 0x6)
    {
        InfoErr<< "Done parse." << nl;
    }

    // Restore debug value
    debug = oldDebug;

    return true;
}


// ************************************************************************* //
