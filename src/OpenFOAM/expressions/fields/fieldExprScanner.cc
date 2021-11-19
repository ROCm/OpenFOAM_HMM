
#line 1 "fieldExprScanner.rl"
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
    Ragel lexer interface for lemon grammar for patch expressions

\*---------------------------------------------------------------------------*/

#include "exprScanToken.H"
#include "fieldExprScanner.H"
#include "fieldExprDriver.H"
#include "fieldExprLemonParser.h"
#include "fieldExprParser.H"
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
#define FIELD_PAIR(Fld,T)   { TOKEN_OF(T), Fld::typeName.c_str() }

// No known look-back types
#undef HAS_LOOKBEHIND_TOKENS

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
    const expressions::fieldExpr::parseDriver& driver_,
    const word& ident
)
{
    #ifdef HAS_LOOKBEHIND_TOKENS
    // Get stashed "look-behind" to decide what type of identifier we expect
    #endif

    // Field variables
    #ifdef TOK_SCALAR_ID
    {
        #undef  doLocalCode
        #define doLocalCode(TokType, Type)                          \
        if (driver_.isLocalVariable<Type>(ident, false))            \
        {                                                           \
            return TokType;                                         \
        }

        doLocalCode(TOK_SCALAR_ID, scalar);
        doLocalCode(TOK_VECTOR_ID, vector);
        doLocalCode(TOK_SYM_TENSOR_ID, symmTensor);
        doLocalCode(TOK_SPH_TENSOR_ID, sphericalTensor);
        doLocalCode(TOK_TENSOR_ID, tensor);
        // Not tested: doLocalCode(TOK_BOOL_ID, bool);
        #undef doLocalCode
    }
    #endif

    return -1;
}

} // End namespace Foam


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names

#define EMIT_TOKEN(T)                                                         \
    driver_.parsePosition() = (ts-buf);                                       \
    DebugInfo<< STRINGIFY(T) << " at " << driver_.parsePosition() << nl;      \
    parser_->parse(TOKEN_OF(T));                                              \
    driver_.parsePosition() = (p-buf);

#define EMIT_VECTOR_TOKEN(X, Y, Z)                                            \
    driver_.parsePosition() = (ts-buf);                                       \
    DebugInfo<< "VECTOR at " << driver_.parsePosition() << nl;                \
    scanToken scanTok;                                                        \
    scanTok.setVector(X,Y,Z);                                                 \
    parser_->parse(TOK_VECTOR_VALUE, scanTok);                                \
    driver_.parsePosition() = (p-buf);



#line 181 "fieldExprScanner.cc"
static const int fieldExpr_start = 13;
static const int fieldExpr_first_final = 13;
static const int fieldExpr_error = 0;

static const int fieldExpr_en_main = 13;


#line 326 "fieldExprScanner.rl"



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expressions::fieldExpr::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::expressions::fieldExpr::scanner::dispatch_method
(
    const parseDriver& driver_,
    word ident
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
        parser_->parse(TOK_DOT);
        parser_->parse(methType);

        return true;
    }

    driver_.reportFatal("Unknown method: " + ident);
    return false;
}


bool Foam::expressions::fieldExpr::scanner::dispatch_ident
(
    const parseDriver& driver_,
    word ident
) const
{
    // Peek at stashed "look-behind". It may influence decisions
    int lookBehindTok = driver_.stashedTokenId();
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

            parser_->parse(tokType);
            return true;
        }

        #ifdef HAS_LOOKBEHIND_TOKENS
        // Specials such "cellSet" etc also reset the look-behind
        tokType = lookBehindTokenEnums.lookup(ident, -1);

        if (tokType > 0)
        {
            DebugInfo
                << "Emit:" << ident << " as look-behind:"
                << parser_->tokenName(tokType) << nl;

            driver_.resetStashedTokenId(tokType);
            parser_->parse(tokType);
            return true;
        }
        #endif
    }

    // Functions: scalar, vector, probably don't need others
    // - "fn:" prefix to avoid any ambiguities
    if (lookBehindTok <= 0 && ident.starts_with("fn:"))
    {
        word funcName(ident.substr(3));  // strip prefix

        do
        {
        }
        while (false);
    }

    if (tokType <= 0)
    {
        tokType = driverTokenType(driver_, ident);
    }

    if (tokType > 0)
    {
        DebugInfo
            << "Emit:" << ident << " token:"
            << parser_->tokenName(tokType) << nl;

        scanToken scanTok;
        scanTok.setWord(ident);
        parser_->parse(tokType, scanTok);

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

        scanToken scanTok;
        scanTok.setWord(ident);
        parser_->parse(tokType, scanTok);

        // Dispatch '.' and "method" separately
        parser_->parse(TOK_DOT);
        parser_->parse(methType);

        return true;
    }

    driver_.reportFatal
    (
        "Object " + ident + " does not exist or wrong type"
    );

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::fieldExpr::scanner::process
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
    
#line 431 "fieldExprScanner.cc"
	{
	cs = fieldExpr_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 566 "fieldExprScanner.rl"
   /* ^^^ FSM initialization here ^^^ */;

    
#line 443 "fieldExprScanner.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr2:
#line 206 "fieldExprScanner.rl"
	{te = p+1;{
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st13;
tr4:
#line 206 "fieldExprScanner.rl"
	{te = p+1;{
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st13;
tr5:
#line 181 "fieldExprScanner.rl"
	{{p = ((te))-1;}{
        // Emit number
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        scanToken scanTok;
        scanTok.setScalar(0);
        if (readScalar(std::string(ts, te-ts), scanTok.scalarValue))
        {
            parser_->parse(TOKEN_OF(NUMBER), scanTok);
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
	goto st13;
tr8:
#line 250 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(EQUAL); }}
	goto st13;
tr9:
#line 302 "fieldExprScanner.rl"
	{{p = ((te))-1;}{ EMIT_TOKEN(TENSOR); }}
	goto st13;
tr11:
#line 313 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(IDENTITY_TENSOR); }}
	goto st13;
tr12:
#line 301 "fieldExprScanner.rl"
	{{p = ((te))-1;}{ EMIT_TOKEN(VECTOR); }}
	goto st13;
tr14:
#line 310 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_VECTOR_TOKEN(1,0,0); }}
	goto st13;
tr15:
#line 311 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_VECTOR_TOKEN(0,1,0); }}
	goto st13;
tr16:
#line 312 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_VECTOR_TOKEN(0,0,1); }}
	goto st13;
tr17:
#line 253 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LOR); }}
	goto st13;
tr21:
#line 235 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PERCENT); }}
	goto st13;
tr24:
#line 236 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LPAREN); }}
	goto st13;
tr25:
#line 237 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(RPAREN); }}
	goto st13;
tr26:
#line 238 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(TIMES); }}
	goto st13;
tr27:
#line 239 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PLUS); }}
	goto st13;
tr28:
#line 241 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COMMA); }}
	goto st13;
tr29:
#line 240 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(MINUS); }}
	goto st13;
tr31:
#line 243 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(DIVIDE); }}
	goto st13;
tr33:
#line 245 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COLON); }}
	goto st13;
tr37:
#line 244 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(QUESTION); }}
	goto st13;
tr40:
#line 256 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(BIT_XOR); }}
	goto st13;
tr56:
#line 229 "fieldExprScanner.rl"
	{te = p;p--;}
	goto st13;
tr57:
#line 234 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LNOT); }}
	goto st13;
tr58:
#line 251 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(NOT_EQUAL); }}
	goto st13;
tr59:
#line 254 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(BIT_AND); }}
	goto st13;
tr60:
#line 252 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LAND); }}
	goto st13;
tr61:
#line 242 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(DOT); }}
	goto st13;
tr64:
#line 181 "fieldExprScanner.rl"
	{te = p;p--;{
        // Emit number
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        scanToken scanTok;
        scanTok.setScalar(0);
        if (readScalar(std::string(ts, te-ts), scanTok.scalarValue))
        {
            parser_->parse(TOKEN_OF(NUMBER), scanTok);
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
	goto st13;
tr66:
#line 213 "fieldExprScanner.rl"
	{te = p;p--;{
        // Tokenized ".method" - dispatch '.' and "method" separately
        driver_.parsePosition() = (ts-buf);
        dispatch_method(driver_, word(ts+1, te-ts-1, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st13;
tr67:
#line 246 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LESS); }}
	goto st13;
tr68:
#line 247 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LESS_EQ); }}
	goto st13;
tr69:
#line 248 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(GREATER); }}
	goto st13;
tr70:
#line 249 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(GREATER_EQ); }}
	goto st13;
tr71:
#line 206 "fieldExprScanner.rl"
	{te = p;p--;{
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st13;
tr73:
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
	{{p = ((te))-1;} EMIT_TOKEN(RAND); }
	break;
	case 58:
	{{p = ((te))-1;} EMIT_TOKEN(BOOL); }
	break;
	case 61:
	{{p = ((te))-1;} EMIT_TOKEN(SYM_TENSOR); }
	break;
	case 62:
	{{p = ((te))-1;} EMIT_TOKEN(SPH_TENSOR); }
	break;
	case 63:
	{{p = ((te))-1;} EMIT_TOKEN(LTRUE); }
	break;
	case 64:
	{{p = ((te))-1;} EMIT_TOKEN(LFALSE); }
	break;
	case 65:
	{{p = ((te))-1;} EMIT_TOKEN(ZERO); }
	break;
	case 70:
	{{p = ((te))-1;} EMIT_TOKEN(ARG); }
	break;
	case 71:
	{{p = ((te))-1;} EMIT_TOKEN(TIME); }
	break;
	case 72:
	{{p = ((te))-1;} EMIT_TOKEN(DELTA_T); }
	break;
	case 73:
	{{p = ((te))-1;}
        // Emit identifier
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }
	break;
	}
	}
	goto st13;
tr89:
#line 278 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(ATAN); }}
	goto st13;
tr104:
#line 274 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(COS); }}
	goto st13;
tr125:
#line 267 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LOG); }}
	goto st13;
tr132:
#line 283 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(MAG); }}
	goto st13;
tr139:
#line 287 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(NEG); }}
	goto st13;
tr145:
#line 286 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(POS); }}
	goto st13;
tr164:
#line 273 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SIN); }}
	goto st13;
tr180:
#line 270 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SQR); }}
	goto st13;
tr196:
#line 275 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TAN); }}
	goto st13;
tr202:
#line 302 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TENSOR); }}
	goto st13;
tr213:
#line 301 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(VECTOR); }}
	goto st13;
st13:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof13;
case 13:
#line 1 "NONE"
	{ts = p;}
#line 817 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 32: goto st14;
		case 33: goto st15;
		case 34: goto st1;
		case 37: goto tr21;
		case 38: goto st16;
		case 39: goto st3;
		case 40: goto tr24;
		case 41: goto tr25;
		case 42: goto tr26;
		case 43: goto tr27;
		case 44: goto tr28;
		case 45: goto tr29;
		case 46: goto st17;
		case 47: goto tr31;
		case 58: goto tr33;
		case 60: goto st22;
		case 61: goto st7;
		case 62: goto st23;
		case 63: goto tr37;
		case 90: goto st26;
		case 94: goto tr40;
		case 95: goto st24;
		case 97: goto st29;
		case 98: goto st43;
		case 99: goto st46;
		case 100: goto st51;
		case 101: goto st61;
		case 102: goto st63;
		case 108: goto st67;
		case 109: goto st71;
		case 110: goto st77;
		case 112: goto st80;
		case 114: goto st83;
		case 115: goto st91;
		case 116: goto st119;
		case 118: goto st131;
		case 124: goto st12;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st14;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 89 ) {
			if ( 103 <= (*p) && (*p) <= 122 )
				goto st24;
		} else if ( (*p) >= 65 )
			goto st24;
	} else
		goto tr32;
	goto st0;
st0:
cs = 0;
	goto _out;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 32 )
		goto st14;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st14;
	goto tr56;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) == 61 )
		goto tr58;
	goto tr57;
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
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	if ( (*p) == 38 )
		goto tr60;
	goto tr59;
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
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr62;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st20;
	} else
		goto st20;
	goto tr61;
tr62:
#line 1 "NONE"
	{te = p+1;}
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
#line 944 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr62;
	goto tr64;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 43: goto st6;
		case 45: goto st6;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st19;
	goto tr5;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st19;
	goto tr5;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st19;
	goto tr64;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st20;
	} else if ( (*p) >= 65 )
		goto st20;
	goto tr66;
tr32:
#line 1 "NONE"
	{te = p+1;}
	goto st21;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
#line 995 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr62;
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr32;
	goto tr64;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 61 )
		goto tr68;
	goto tr67;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 61 )
		goto tr8;
	goto st0;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 61 )
		goto tr70;
	goto tr69;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
tr72:
#line 1 "NONE"
	{te = p+1;}
#line 206 "fieldExprScanner.rl"
	{act = 73;}
	goto st25;
tr76:
#line 1 "NONE"
	{te = p+1;}
#line 309 "fieldExprScanner.rl"
	{act = 65;}
	goto st25;
tr83:
#line 1 "NONE"
	{te = p+1;}
#line 277 "fieldExprScanner.rl"
	{act = 40;}
	goto st25;
tr84:
#line 1 "NONE"
	{te = p+1;}
#line 314 "fieldExprScanner.rl"
	{act = 70;}
	goto st25;
tr86:
#line 1 "NONE"
	{te = p+1;}
#line 276 "fieldExprScanner.rl"
	{act = 39;}
	goto st25;
tr90:
#line 1 "NONE"
	{te = p+1;}
#line 279 "fieldExprScanner.rl"
	{act = 42;}
	goto st25;
tr95:
#line 1 "NONE"
	{te = p+1;}
#line 295 "fieldExprScanner.rl"
	{act = 55;}
	goto st25;
tr98:
#line 1 "NONE"
	{te = p+1;}
#line 300 "fieldExprScanner.rl"
	{act = 58;}
	goto st25;
tr102:
#line 1 "NONE"
	{te = p+1;}
#line 272 "fieldExprScanner.rl"
	{act = 35;}
	goto st25;
tr105:
#line 1 "NONE"
	{te = p+1;}
#line 281 "fieldExprScanner.rl"
	{act = 44;}
	goto st25;
tr113:
#line 1 "NONE"
	{te = p+1;}
#line 264 "fieldExprScanner.rl"
	{act = 27;}
	goto st25;
tr116:
#line 1 "NONE"
	{te = p+1;}
#line 316 "fieldExprScanner.rl"
	{act = 72;}
	goto st25;
tr118:
#line 1 "NONE"
	{te = p+1;}
#line 266 "fieldExprScanner.rl"
	{act = 29;}
	goto st25;
tr122:
#line 1 "NONE"
	{te = p+1;}
#line 308 "fieldExprScanner.rl"
	{act = 64;}
	goto st25;
tr127:
#line 1 "NONE"
	{te = p+1;}
#line 268 "fieldExprScanner.rl"
	{act = 31;}
	goto st25;
tr131:
#line 1 "NONE"
	{te = p+1;}
#line 294 "fieldExprScanner.rl"
	{act = 54;}
	goto st25;
tr135:
#line 1 "NONE"
	{te = p+1;}
#line 284 "fieldExprScanner.rl"
	{act = 47;}
	goto st25;
tr136:
#line 1 "NONE"
	{te = p+1;}
#line 293 "fieldExprScanner.rl"
	{act = 53;}
	goto st25;
tr140:
#line 1 "NONE"
	{te = p+1;}
#line 289 "fieldExprScanner.rl"
	{act = 51;}
	goto st25;
tr141:
#line 1 "NONE"
	{te = p+1;}
#line 263 "fieldExprScanner.rl"
	{act = 26;}
	goto st25;
tr144:
#line 1 "NONE"
	{te = p+1;}
#line 269 "fieldExprScanner.rl"
	{act = 32;}
	goto st25;
tr146:
#line 1 "NONE"
	{te = p+1;}
#line 288 "fieldExprScanner.rl"
	{act = 50;}
	goto st25;
tr154:
#line 1 "NONE"
	{te = p+1;}
#line 265 "fieldExprScanner.rl"
	{act = 28;}
	goto st25;
tr155:
#line 1 "NONE"
	{te = p+1;}
#line 297 "fieldExprScanner.rl"
	{act = 57;}
	goto st25;
tr163:
#line 1 "NONE"
	{te = p+1;}
#line 290 "fieldExprScanner.rl"
	{act = 52;}
	goto st25;
tr165:
#line 1 "NONE"
	{te = p+1;}
#line 280 "fieldExprScanner.rl"
	{act = 43;}
	goto st25;
tr178:
#line 1 "NONE"
	{te = p+1;}
#line 304 "fieldExprScanner.rl"
	{act = 62;}
	goto st25;
tr181:
#line 1 "NONE"
	{te = p+1;}
#line 271 "fieldExprScanner.rl"
	{act = 34;}
	goto st25;
tr182:
#line 1 "NONE"
	{te = p+1;}
#line 296 "fieldExprScanner.rl"
	{act = 56;}
	goto st25;
tr190:
#line 1 "NONE"
	{te = p+1;}
#line 303 "fieldExprScanner.rl"
	{act = 61;}
	goto st25;
tr197:
#line 1 "NONE"
	{te = p+1;}
#line 282 "fieldExprScanner.rl"
	{act = 45;}
	goto st25;
tr205:
#line 1 "NONE"
	{te = p+1;}
#line 315 "fieldExprScanner.rl"
	{act = 71;}
	goto st25;
tr207:
#line 1 "NONE"
	{te = p+1;}
#line 307 "fieldExprScanner.rl"
	{act = 63;}
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 1244 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr73;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st27;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto st28;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto tr76;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 99: goto st30;
		case 114: goto st32;
		case 115: goto st33;
		case 116: goto st35;
		case 118: goto st38;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st31;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 115: goto tr83;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto tr84;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 105: goto st34;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto tr86;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 46: goto tr72;
		case 50: goto tr90;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr89;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st39;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto st40;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st41;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto st42;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto tr95;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st44;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st45;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 108: goto tr98;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 98: goto st47;
		case 111: goto st49;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto st48;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 116: goto tr102;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 115: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 104: goto tr105;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr104;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st52;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto st53;
		case 108: goto st58;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 46: goto tr72;
		case 84: goto st54;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st55;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 46: goto tr72;
		case 82: goto st56;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st57;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 100: goto tr113;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 116: goto st59;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st60;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 46: goto tr72;
		case 84: goto tr116;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 120: goto st62;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 112: goto tr118;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st64;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 108: goto st65;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 115: goto st66;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto tr122;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st68;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto st69;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 46: goto tr72;
		case 49: goto st70;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr125;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 46: goto tr72;
		case 48: goto tr127;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st72;
		case 105: goto st76;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto st73;
		case 120: goto tr131;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 46: goto tr72;
		case 83: goto st74;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr132;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 113: goto st75;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto tr135;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto tr136;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st78;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 46: goto tr72;
		case 48: goto tr140;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr139;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 105: goto tr141;
		case 111: goto st81;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 115: goto st82;
		case 119: goto tr144;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	switch( (*p) ) {
		case 46: goto tr72;
		case 48: goto tr146;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr145;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st84;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 100: goto st85;
		case 110: goto st90;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 46: goto tr72;
		case 84: goto st86;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st87;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	switch( (*p) ) {
		case 46: goto tr72;
		case 68: goto st88;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st89;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto tr154;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 100: goto tr155;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 105: goto st92;
		case 112: goto st95;
		case 113: goto st108;
		case 117: goto st110;
		case 121: goto st111;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 103: goto st93;
		case 110: goto st94;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto tr163;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 104: goto tr165;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr164;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 104: goto st96;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto st98;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 105: goto st99;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 99: goto st100;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st101;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 108: goto st102;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
	switch( (*p) ) {
		case 46: goto tr72;
		case 84: goto st103;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st104;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto st105;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 115: goto st106;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st107;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto tr178;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto st109;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 116: goto tr181;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr180;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 109: goto tr182;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 109: goto st112;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 109: goto st113;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 46: goto tr72;
		case 84: goto st114;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st115;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto st116;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 115: goto st117;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto tr190;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 97: goto st120;
		case 101: goto st122;
		case 105: goto st127;
		case 114: goto st129;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto st121;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 104: goto tr197;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr196;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 110: goto st123;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 115: goto st124;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto tr201;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
tr201:
#line 1 "NONE"
	{te = p+1;}
	goto st126;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
#line 3085 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr72;
		case 58: goto st8;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr202;
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
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 109: goto st128;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto tr205;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st129:
	if ( ++p == pe )
		goto _test_eof129;
case 129:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 117: goto st130;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st130:
	if ( ++p == pe )
		goto _test_eof130;
case 130:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto tr207;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st131:
	if ( ++p == pe )
		goto _test_eof131;
case 131:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 101: goto st132;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st132:
	if ( ++p == pe )
		goto _test_eof132;
case 132:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 99: goto st133;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st133:
	if ( ++p == pe )
		goto _test_eof133;
case 133:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 116: goto st134;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st134:
	if ( ++p == pe )
		goto _test_eof134;
case 134:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 111: goto st135;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
st135:
	if ( ++p == pe )
		goto _test_eof135;
case 135:
	switch( (*p) ) {
		case 46: goto tr72;
		case 95: goto tr72;
		case 114: goto tr212;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr71;
tr212:
#line 1 "NONE"
	{te = p+1;}
	goto st136;
st136:
	if ( ++p == pe )
		goto _test_eof136;
case 136:
#line 3284 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr72;
		case 58: goto st10;
		case 95: goto tr72;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr72;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr72;
	} else
		goto tr72;
	goto tr213;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 58 )
		goto st11;
	goto tr12;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	switch( (*p) ) {
		case 120: goto tr14;
		case 121: goto tr15;
		case 122: goto tr16;
	}
	goto tr12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 124 )
		goto tr17;
	goto st0;
	}
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
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
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
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
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 14: goto tr56;
	case 15: goto tr57;
	case 16: goto tr59;
	case 17: goto tr61;
	case 18: goto tr64;
	case 5: goto tr5;
	case 6: goto tr5;
	case 19: goto tr64;
	case 20: goto tr66;
	case 21: goto tr64;
	case 22: goto tr67;
	case 23: goto tr69;
	case 24: goto tr71;
	case 25: goto tr73;
	case 26: goto tr71;
	case 27: goto tr71;
	case 28: goto tr71;
	case 29: goto tr71;
	case 30: goto tr71;
	case 31: goto tr71;
	case 32: goto tr71;
	case 33: goto tr71;
	case 34: goto tr71;
	case 35: goto tr71;
	case 36: goto tr71;
	case 37: goto tr89;
	case 38: goto tr71;
	case 39: goto tr71;
	case 40: goto tr71;
	case 41: goto tr71;
	case 42: goto tr71;
	case 43: goto tr71;
	case 44: goto tr71;
	case 45: goto tr71;
	case 46: goto tr71;
	case 47: goto tr71;
	case 48: goto tr71;
	case 49: goto tr71;
	case 50: goto tr104;
	case 51: goto tr71;
	case 52: goto tr71;
	case 53: goto tr71;
	case 54: goto tr71;
	case 55: goto tr71;
	case 56: goto tr71;
	case 57: goto tr71;
	case 58: goto tr71;
	case 59: goto tr71;
	case 60: goto tr71;
	case 61: goto tr71;
	case 62: goto tr71;
	case 63: goto tr71;
	case 64: goto tr71;
	case 65: goto tr71;
	case 66: goto tr71;
	case 67: goto tr71;
	case 68: goto tr71;
	case 69: goto tr125;
	case 70: goto tr71;
	case 71: goto tr71;
	case 72: goto tr71;
	case 73: goto tr132;
	case 74: goto tr71;
	case 75: goto tr71;
	case 76: goto tr71;
	case 77: goto tr71;
	case 78: goto tr71;
	case 79: goto tr139;
	case 80: goto tr71;
	case 81: goto tr71;
	case 82: goto tr145;
	case 83: goto tr71;
	case 84: goto tr71;
	case 85: goto tr71;
	case 86: goto tr71;
	case 87: goto tr71;
	case 88: goto tr71;
	case 89: goto tr71;
	case 90: goto tr71;
	case 91: goto tr71;
	case 92: goto tr71;
	case 93: goto tr71;
	case 94: goto tr164;
	case 95: goto tr71;
	case 96: goto tr71;
	case 97: goto tr71;
	case 98: goto tr71;
	case 99: goto tr71;
	case 100: goto tr71;
	case 101: goto tr71;
	case 102: goto tr71;
	case 103: goto tr71;
	case 104: goto tr71;
	case 105: goto tr71;
	case 106: goto tr71;
	case 107: goto tr71;
	case 108: goto tr71;
	case 109: goto tr180;
	case 110: goto tr71;
	case 111: goto tr71;
	case 112: goto tr71;
	case 113: goto tr71;
	case 114: goto tr71;
	case 115: goto tr71;
	case 116: goto tr71;
	case 117: goto tr71;
	case 118: goto tr71;
	case 119: goto tr71;
	case 120: goto tr71;
	case 121: goto tr196;
	case 122: goto tr71;
	case 123: goto tr71;
	case 124: goto tr71;
	case 125: goto tr71;
	case 126: goto tr202;
	case 8: goto tr9;
	case 9: goto tr9;
	case 127: goto tr71;
	case 128: goto tr71;
	case 129: goto tr71;
	case 130: goto tr71;
	case 131: goto tr71;
	case 132: goto tr71;
	case 133: goto tr71;
	case 134: goto tr71;
	case 135: goto tr71;
	case 136: goto tr213;
	case 10: goto tr12;
	case 11: goto tr12;
	}
	}

	_out: {}
	}

#line 568 "fieldExprScanner.rl"
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
    parser_->parse(0);
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
