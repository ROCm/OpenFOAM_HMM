divert(-1)dnl
#-----------------------------------*- m4 -*-----------------------------------
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | www.openfoam.com
#      \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2019 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GNU General Public
#     License GPL-3.0 or later <https://www.gnu.org/licenses/gpl-3.0>
#
# Description
#     Defines standard operator precedence macro for lemon grammar.
#     Follows https://en.cppreference.com/w/cpp/language/operator_precedence
#
#------------------------------------------------------------------------------

define([operator_precedence],
[// [https://en.cppreference.com/w/cpp/language/operator_precedence]
%right QUESTION COLON .                 // 16: right-to-left
%left LOR  .                            // 15:
%left LAND .                            // 14:
// %left BIT_OR  .                        // 13 (unused)
%left BIT_XOR .                         // 12
%left BIT_AND .                         // 11
%left EQUAL NOT_EQUAL .                 // 10
%left LESS_EQ GREATER_EQ LESS GREATER . // 9
%left PLUS MINUS .                      // 6
%left TIMES DIVIDE PERCENT .            // 5
%right NEGATE NOT .                     // 3: right-to-left
%left DOT  .                            // 2: (method)]
)

#------------------------------------------------------------------------------
divert(0)dnl
