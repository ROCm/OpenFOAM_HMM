
#line 1 "STLAsciiParseRagel.rl"
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

Description
    Ragel-based parsing of STL ASCII format.
    The goto-based finite state machine (FSM) is generated with

        ragel -G2 -o STLAsciiParseRagel.C STLAsciiParseRagel.rl

\*---------------------------------------------------------------------------*/

#include "STLAsciiParse.H"
#include "STLReader.H"
#include "OSspecific.H"

// https://en.wikipedia.org/wiki/STL_%28file_format%29#ASCII_STL
//
// Format
//
// solid [name]
//
// * where name is an optional string.
// * The file continues with any number of triangles,
//   each represented as follows:
//
// [color ...]
// facet normal ni nj nk
//     outer loop
//         vertex v1x v1y v1z
//         vertex v2x v2y v2z
//         vertex v3x v3y v3z
//     endloop
// endfacet
//
// * where each n or v is a floating-point number.
// * The file concludes with
//
// endsolid [name]

// We take some parsing shortcuts.
// - Ignore 'color' lines
// - Only look for initial 'facet '. Ignore 'normal ...'
// - Ignore name for 'endsolid'
//

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names

// Define the machine actions

#line 92 "STLAsciiParseRagel.rl"




#line 133 "STLAsciiParseRagel.rl"



//
// FSM globals
//


#line 90 "STLAsciiParseRagel.C"
static const int stlAscii_start = 1;
static const int stlAscii_error = 0;

static const int stlAscii_en_main = 1;


#line 141 "STLAsciiParseRagel.rl"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Foam
{
namespace Detail
{

//- A lexer for parsing STL ASCII files.
//  Returns DynamicList(s) of points and facets (zoneIds).
//  The facets are within a solid/endsolid grouping
class STLAsciiParseRagel
:
    public Detail::STLAsciiParse
{
    // Private Data
    word  startError_;


public:

    //- From input stream and the approximate number of vertices in the STL
    STLAsciiParseRagel(const label approxNpoints)
    :
        Detail::STLAsciiParse(approxNpoints)
    {}

    //- Execute lexer
    void execute(std::istream& is);
};

} // end of namespace Detail
} // end of namespace Foam


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Length of the input read buffer
#define INBUFLEN 16384

void Foam::Detail::STLAsciiParseRagel::execute(std::istream& is)
{
    if (!is)
    {
        return;
    }

    // cs = code state
    int cs;

    
#line 150 "STLAsciiParseRagel.C"
	{
	cs = stlAscii_start;
	}

#line 192 "STLAsciiParseRagel.rl"
   /* ^^^ FSM initialization here ^^^ */;

    // Local token start
    char *tok = nullptr;

    // Buffering
    char inbuf[INBUFLEN];
    std::streamsize pending = 0;

    // Line-oriented processing loop (as per Ragel pdf example)

    while (is)
    {
        if (pending >= INBUFLEN)
        {
            // We overfilled the buffer while trying to scan a token...
            FatalErrorInFunction
                << "buffer full while scanning near line " << lineNum_ << nl;
            break;
        }

        char *data = inbuf + pending;   // current data buffer
        const std::streamsize buflen = INBUFLEN - pending; // space in buffer

        is.read(data, buflen);
        const std::streamsize gcount = is.gcount();
        if (!gcount)
        {
            break;
        }

        // p,pe = Ragel parsing point, parsing end (default naming)
        // eof  = Ragel EOF point (default naming)

        char *p = inbuf;
        const char *pe = data + gcount;
        const char *eof = nullptr;
        if (!is)
        {
            eof = pe;   // Tag 'pe' as being the EOF for the FSM as well
        }

        // Line-oriented: search backwards to find last newline
        {
            --pe;
            while (*pe != '\n' && pe >= inbuf)
            {
                --pe;
            }
            ++pe;
        }

        
#line 209 "STLAsciiParseRagel.C"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	switch( (*p) ) {
		case 32: goto st1;
		case 83: goto st2;
		case 115: goto st25;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st1;
	goto st0;
st0:
cs = 0;
	goto _out;
tr177:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 238 "STLAsciiParseRagel.C"
	if ( (*p) == 79 )
		goto st3;
	goto st0;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( (*p) == 76 )
		goto st4;
	goto st0;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( (*p) == 73 )
		goto st5;
	goto st0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( (*p) == 68 )
		goto st6;
	goto st0;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	switch( (*p) ) {
		case 9: goto tr8;
		case 10: goto tr9;
		case 32: goto tr8;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto tr8;
	goto st0;
tr8:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
#line 78 "STLAsciiParseRagel.rl"
	{ beginSolid(word::validate(tok, p)); }
	goto st7;
tr11:
#line 78 "STLAsciiParseRagel.rl"
	{ beginSolid(word::validate(tok, p)); }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 289 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 9: goto tr11;
		case 10: goto tr12;
		case 32: goto tr11;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto tr11;
	goto st7;
tr9:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
#line 78 "STLAsciiParseRagel.rl"
	{ beginSolid(word::validate(tok, p)); }
	goto st8;
tr12:
#line 78 "STLAsciiParseRagel.rl"
	{ beginSolid(word::validate(tok, p)); }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 312 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto tr13;
		case 67: goto tr14;
		case 69: goto tr15;
		case 70: goto tr16;
		case 99: goto tr17;
		case 101: goto tr18;
		case 102: goto tr19;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr13;
	goto st0;
tr13:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 333 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto st9;
		case 67: goto st10;
		case 69: goto st17;
		case 70: goto st29;
		case 99: goto st136;
		case 101: goto st125;
		case 102: goto st132;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st9;
	goto st0;
tr14:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 354 "STLAsciiParseRagel.C"
	if ( (*p) == 79 )
		goto st11;
	goto st0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( (*p) == 76 )
		goto st12;
	goto st0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 79 )
		goto st13;
	goto st0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 82 )
		goto st14;
	goto st0;
tr92:
#line 80 "STLAsciiParseRagel.rl"
	{ endFacet(); }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 387 "STLAsciiParseRagel.C"
	if ( (*p) == 10 )
		goto st15;
	goto st14;
tr93:
#line 80 "STLAsciiParseRagel.rl"
	{ endFacet(); }
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
#line 399 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto tr32;
		case 69: goto tr15;
		case 70: goto tr16;
		case 101: goto tr18;
		case 102: goto tr19;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr32;
	goto st0;
tr32:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 418 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto st16;
		case 69: goto st17;
		case 70: goto st29;
		case 101: goto st125;
		case 102: goto st132;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st16;
	goto st0;
tr15:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
#line 437 "STLAsciiParseRagel.C"
	if ( (*p) == 78 )
		goto st18;
	goto st0;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( (*p) == 68 )
		goto st19;
	goto st0;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( (*p) == 83 )
		goto st20;
	goto st0;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 79 )
		goto st21;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 76 )
		goto st22;
	goto st0;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 73 )
		goto st23;
	goto st0;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 68 )
		goto st24;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 10 )
		goto st140;
	goto st24;
st140:
	if ( ++p == pe )
		goto _test_eof140;
case 140:
	switch( (*p) ) {
		case 32: goto tr176;
		case 83: goto tr177;
		case 115: goto tr178;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr176;
	goto st0;
tr176:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st141;
st141:
	if ( ++p == pe )
		goto _test_eof141;
case 141:
#line 510 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto st141;
		case 83: goto st2;
		case 115: goto st25;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st141;
	goto st0;
tr178:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 527 "STLAsciiParseRagel.C"
	if ( (*p) == 111 )
		goto st26;
	goto st0;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( (*p) == 108 )
		goto st27;
	goto st0;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	if ( (*p) == 105 )
		goto st28;
	goto st0;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	if ( (*p) == 100 )
		goto st6;
	goto st0;
tr16:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st29;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
#line 560 "STLAsciiParseRagel.C"
	if ( (*p) == 65 )
		goto st30;
	goto st0;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( (*p) == 67 )
		goto st31;
	goto st0;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	if ( (*p) == 69 )
		goto st32;
	goto st0;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	if ( (*p) == 84 )
		goto st33;
	goto st0;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 9: goto st34;
		case 32: goto st34;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st34;
	goto st0;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	if ( (*p) == 10 )
		goto tr51;
	goto tr50;
tr50:
#line 79 "STLAsciiParseRagel.rl"
	{ beginFacet(); }
	goto st35;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
#line 611 "STLAsciiParseRagel.C"
	if ( (*p) == 10 )
		goto st36;
	goto st35;
tr51:
#line 79 "STLAsciiParseRagel.rl"
	{ beginFacet(); }
	goto st36;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
#line 623 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto tr54;
		case 79: goto tr55;
		case 111: goto tr56;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr54;
	goto st0;
tr54:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st37;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
#line 640 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto st37;
		case 111: goto st38;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st37;
	goto st0;
tr56:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st38;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
#line 656 "STLAsciiParseRagel.C"
	if ( (*p) == 117 )
		goto st39;
	goto st0;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	if ( (*p) == 116 )
		goto st40;
	goto st0;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	if ( (*p) == 101 )
		goto st41;
	goto st0;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	if ( (*p) == 114 )
		goto st42;
	goto st0;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 9: goto st43;
		case 32: goto st43;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st43;
	goto st0;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 9: goto st43;
		case 32: goto st43;
		case 108: goto st44;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st43;
	goto st0;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	if ( (*p) == 111 )
		goto st45;
	goto st0;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	if ( (*p) == 111 )
		goto st46;
	goto st0;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	if ( (*p) == 112 )
		goto st47;
	goto st0;
tr135:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st47;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
#line 733 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto st47;
		case 69: goto st48;
		case 86: goto st72;
		case 101: goto st92;
		case 118: goto st98;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st47;
	goto st0;
tr136:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st48;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
#line 752 "STLAsciiParseRagel.C"
	if ( (*p) == 78 )
		goto st49;
	goto st0;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	if ( (*p) == 68 )
		goto st50;
	goto st0;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	if ( (*p) == 76 )
		goto st51;
	goto st0;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	if ( (*p) == 79 )
		goto st52;
	goto st0;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	if ( (*p) == 79 )
		goto st53;
	goto st0;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	if ( (*p) == 80 )
		goto st54;
	goto st0;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	if ( (*p) == 10 )
		goto st55;
	goto st54;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 32: goto tr79;
		case 69: goto tr80;
		case 101: goto tr81;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr79;
	goto st0;
tr79:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st56;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
#line 818 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto st56;
		case 69: goto st57;
		case 101: goto st65;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st56;
	goto st0;
tr80:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st57;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
#line 835 "STLAsciiParseRagel.C"
	if ( (*p) == 78 )
		goto st58;
	goto st0;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	if ( (*p) == 68 )
		goto st59;
	goto st0;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	if ( (*p) == 70 )
		goto st60;
	goto st0;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	if ( (*p) == 65 )
		goto st61;
	goto st0;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	if ( (*p) == 67 )
		goto st62;
	goto st0;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	if ( (*p) == 69 )
		goto st63;
	goto st0;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	if ( (*p) == 84 )
		goto st64;
	goto st0;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	if ( (*p) == 10 )
		goto tr93;
	goto tr92;
tr81:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
#line 896 "STLAsciiParseRagel.C"
	if ( (*p) == 110 )
		goto st66;
	goto st0;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	if ( (*p) == 100 )
		goto st67;
	goto st0;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	if ( (*p) == 102 )
		goto st68;
	goto st0;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	if ( (*p) == 97 )
		goto st69;
	goto st0;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	if ( (*p) == 99 )
		goto st70;
	goto st0;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	if ( (*p) == 101 )
		goto st71;
	goto st0;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	if ( (*p) == 116 )
		goto st64;
	goto st0;
tr137:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st72;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
#line 950 "STLAsciiParseRagel.C"
	if ( (*p) == 69 )
		goto st73;
	goto st0;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	if ( (*p) == 82 )
		goto st74;
	goto st0;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	if ( (*p) == 84 )
		goto st75;
	goto st0;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	if ( (*p) == 69 )
		goto st76;
	goto st0;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	if ( (*p) == 88 )
		goto st77;
	goto st0;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	switch( (*p) ) {
		case 9: goto st78;
		case 32: goto st78;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st78;
	goto st0;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 9: goto st78;
		case 32: goto st78;
		case 43: goto tr106;
		case 45: goto tr106;
		case 46: goto tr107;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr108;
	} else if ( (*p) >= 12 )
		goto st78;
	goto st0;
tr106:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st79;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
#line 1018 "STLAsciiParseRagel.C"
	if ( (*p) == 46 )
		goto st80;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st114;
	goto st0;
tr107:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st80;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
#line 1032 "STLAsciiParseRagel.C"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st81;
	goto st0;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	switch( (*p) ) {
		case 9: goto tr112;
		case 32: goto tr112;
		case 69: goto st111;
		case 101: goto st111;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st81;
	} else if ( (*p) >= 12 )
		goto tr112;
	goto st0;
tr112:
#line 84 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st82;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
#line 1066 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 9: goto st82;
		case 32: goto st82;
		case 43: goto tr115;
		case 45: goto tr115;
		case 46: goto tr116;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr117;
	} else if ( (*p) >= 12 )
		goto st82;
	goto st0;
tr115:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st83;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
#line 1088 "STLAsciiParseRagel.C"
	if ( (*p) == 46 )
		goto st84;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st110;
	goto st0;
tr116:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st84;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
#line 1102 "STLAsciiParseRagel.C"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st85;
	goto st0;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 9: goto tr121;
		case 32: goto tr121;
		case 69: goto st107;
		case 101: goto st107;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st85;
	} else if ( (*p) >= 12 )
		goto tr121;
	goto st0;
tr121:
#line 84 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st86;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
#line 1136 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 9: goto st86;
		case 32: goto st86;
		case 43: goto tr124;
		case 45: goto tr124;
		case 46: goto tr125;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr126;
	} else if ( (*p) >= 12 )
		goto st86;
	goto st0;
tr124:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st87;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
#line 1158 "STLAsciiParseRagel.C"
	if ( (*p) == 46 )
		goto st88;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st106;
	goto st0;
tr125:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st88;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
#line 1172 "STLAsciiParseRagel.C"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st89;
	goto st0;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 9: goto tr130;
		case 10: goto tr131;
		case 32: goto tr130;
		case 69: goto st103;
		case 101: goto st103;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st89;
	} else if ( (*p) >= 12 )
		goto tr130;
	goto st0;
tr130:
#line 84 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st90;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
#line 1207 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 9: goto st90;
		case 10: goto st91;
		case 32: goto st90;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st90;
	goto st0;
tr131:
#line 84 "STLAsciiParseRagel.rl"
	{
        const char saveC = *p;
        *p = '\0';  // Make nul-terminated

        addVertexComponent(tok);
        *p = saveC; // Restore previous character
    }
	goto st91;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
#line 1230 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 32: goto tr135;
		case 69: goto tr136;
		case 86: goto tr137;
		case 101: goto tr138;
		case 118: goto tr139;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto tr135;
	goto st0;
tr138:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st92;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
#line 1249 "STLAsciiParseRagel.C"
	if ( (*p) == 110 )
		goto st93;
	goto st0;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	if ( (*p) == 100 )
		goto st94;
	goto st0;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	if ( (*p) == 108 )
		goto st95;
	goto st0;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
	if ( (*p) == 111 )
		goto st96;
	goto st0;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	if ( (*p) == 111 )
		goto st97;
	goto st0;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	if ( (*p) == 112 )
		goto st54;
	goto st0;
tr139:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st98;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
#line 1296 "STLAsciiParseRagel.C"
	if ( (*p) == 101 )
		goto st99;
	goto st0;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
	if ( (*p) == 114 )
		goto st100;
	goto st0;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
	if ( (*p) == 116 )
		goto st101;
	goto st0;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	if ( (*p) == 101 )
		goto st102;
	goto st0;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
	if ( (*p) == 120 )
		goto st77;
	goto st0;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	switch( (*p) ) {
		case 43: goto st104;
		case 45: goto st104;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st105;
	goto st0;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st105;
	goto st0;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	switch( (*p) ) {
		case 9: goto tr130;
		case 10: goto tr131;
		case 32: goto tr130;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st105;
	} else if ( (*p) >= 12 )
		goto tr130;
	goto st0;
tr126:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st106;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
#line 1369 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 9: goto tr130;
		case 10: goto tr131;
		case 32: goto tr130;
		case 46: goto st89;
		case 69: goto st103;
		case 101: goto st103;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st106;
	} else if ( (*p) >= 12 )
		goto tr130;
	goto st0;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
	switch( (*p) ) {
		case 43: goto st108;
		case 45: goto st108;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st109;
	goto st0;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st109;
	goto st0;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
	switch( (*p) ) {
		case 9: goto tr121;
		case 32: goto tr121;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st109;
	} else if ( (*p) >= 12 )
		goto tr121;
	goto st0;
tr117:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st110;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
#line 1424 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 9: goto tr121;
		case 32: goto tr121;
		case 46: goto st85;
		case 69: goto st107;
		case 101: goto st107;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st110;
	} else if ( (*p) >= 12 )
		goto tr121;
	goto st0;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
	switch( (*p) ) {
		case 43: goto st112;
		case 45: goto st112;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st113;
	goto st0;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st113;
	goto st0;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 9: goto tr112;
		case 32: goto tr112;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st113;
	} else if ( (*p) >= 12 )
		goto tr112;
	goto st0;
tr108:
#line 75 "STLAsciiParseRagel.rl"
	{ tok = p; /* Local token start */ }
	goto st114;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
#line 1478 "STLAsciiParseRagel.C"
	switch( (*p) ) {
		case 9: goto tr112;
		case 32: goto tr112;
		case 46: goto st81;
		case 69: goto st111;
		case 101: goto st111;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st114;
	} else if ( (*p) >= 12 )
		goto tr112;
	goto st0;
tr55:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st115;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
#line 1500 "STLAsciiParseRagel.C"
	if ( (*p) == 85 )
		goto st116;
	goto st0;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
	if ( (*p) == 84 )
		goto st117;
	goto st0;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
	if ( (*p) == 69 )
		goto st118;
	goto st0;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
	if ( (*p) == 82 )
		goto st119;
	goto st0;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	switch( (*p) ) {
		case 9: goto st120;
		case 32: goto st120;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st120;
	goto st0;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	switch( (*p) ) {
		case 9: goto st120;
		case 32: goto st120;
		case 76: goto st121;
	}
	if ( 12 <= (*p) && (*p) <= 13 )
		goto st120;
	goto st0;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
	if ( (*p) == 79 )
		goto st122;
	goto st0;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	if ( (*p) == 79 )
		goto st123;
	goto st0;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	if ( (*p) == 80 )
		goto st124;
	goto st0;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
	if ( (*p) == 10 )
		goto st91;
	goto st124;
tr18:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st125;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
#line 1584 "STLAsciiParseRagel.C"
	if ( (*p) == 110 )
		goto st126;
	goto st0;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
	if ( (*p) == 100 )
		goto st127;
	goto st0;
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
	if ( (*p) == 115 )
		goto st128;
	goto st0;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
	if ( (*p) == 111 )
		goto st129;
	goto st0;
st129:
	if ( ++p == pe )
		goto _test_eof129;
case 129:
	if ( (*p) == 108 )
		goto st130;
	goto st0;
st130:
	if ( ++p == pe )
		goto _test_eof130;
case 130:
	if ( (*p) == 105 )
		goto st131;
	goto st0;
st131:
	if ( ++p == pe )
		goto _test_eof131;
case 131:
	if ( (*p) == 100 )
		goto st24;
	goto st0;
tr19:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st132;
st132:
	if ( ++p == pe )
		goto _test_eof132;
case 132:
#line 1638 "STLAsciiParseRagel.C"
	if ( (*p) == 97 )
		goto st133;
	goto st0;
st133:
	if ( ++p == pe )
		goto _test_eof133;
case 133:
	if ( (*p) == 99 )
		goto st134;
	goto st0;
st134:
	if ( ++p == pe )
		goto _test_eof134;
case 134:
	if ( (*p) == 101 )
		goto st135;
	goto st0;
st135:
	if ( ++p == pe )
		goto _test_eof135;
case 135:
	if ( (*p) == 116 )
		goto st33;
	goto st0;
tr17:
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	goto st136;
st136:
	if ( ++p == pe )
		goto _test_eof136;
case 136:
#line 1671 "STLAsciiParseRagel.C"
	if ( (*p) == 111 )
		goto st137;
	goto st0;
st137:
	if ( ++p == pe )
		goto _test_eof137;
case 137:
	if ( (*p) == 108 )
		goto st138;
	goto st0;
st138:
	if ( ++p == pe )
		goto _test_eof138;
case 138:
	if ( (*p) == 111 )
		goto st139;
	goto st0;
st139:
	if ( ++p == pe )
		goto _test_eof139;
case 139:
	if ( (*p) == 114 )
		goto st14;
	goto st0;
	}
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
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

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 140: 
#line 76 "STLAsciiParseRagel.rl"
	{ ++lineNum_; }
	break;
#line 1847 "STLAsciiParseRagel.C"
	}
	}

	_out: {}
	}

#line 244 "STLAsciiParseRagel.rl"
       /* ^^^ FSM execution here ^^^ */;

        if (0 == cs)
        {
            // FSM failed before finding a token
            FatalErrorInFunction
                << "parse error while scanning near line " << lineNum_ << nl;

            if (p)
            {
                std::string::size_type errLen = (pe - p);
                if (errLen > 80)
                {
                    errLen = 80;
                }

                FatalErrorInFunction
                    << "context: " << std::string(p, errLen) << nl;
            }
            break;
        }

        // How much still in the bufffer?
        pending = data + gcount - pe;

        if (pending)
        {
            memmove(inbuf, pe, pending);
        }

        if (gcount < buflen)
        {
            break; // done
        }
    }
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//
// Member Function
//
bool Foam::fileFormats::STLReader::readAsciiRagel
(
    const fileName& filename
)
{
    IFstream is(filename);
    if (!is)
    {
        FatalErrorInFunction
            << "file " << filename << " not found"
            << exit(FatalError);
    }

    // Create with approx number of vertices in the STL (from file size)
    Detail::STLAsciiParseRagel lexer(Foam::fileSize(filename)/400);
    lexer.execute(is.stdStream());

    transfer(lexer);

    return true;
}


// ************************************************************************* //
