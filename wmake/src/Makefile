#----------------------------*- makefile-gmake -*------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     Makefile
#
# Description
#     Makefile used to compile wmake utilities.
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# The Makefile use a POSIX shell
#------------------------------------------------------------------------------

SHELL      = /bin/sh

#------------------------------------------------------------------------------
# Unset suffices list (suffix rules are not used)
#------------------------------------------------------------------------------

.SUFFIXES:


#------------------------------------------------------------------------------
# Set compilation and dependency building rules
#------------------------------------------------------------------------------

# Locally set optimized compilation
WM_COMPILE_OPTION = Opt

GENERAL_RULES   = $(WM_DIR)/rules/General
include $(GENERAL_RULES)/general


#------------------------------------------------------------------------------
# Targets
#------------------------------------------------------------------------------

.PHONY: all clean

all:  $(WMAKE_BIN)/dirToString $(WMAKE_BIN)/wmkdep
	@echo built wmake-bin for $(WM_ARCH)$(WM_COMPILER)

clean:
	@echo clean wmake-bin for $(WM_ARCH)$(WM_COMPILER)
	@rm -rf $(WMAKE_BIN) 2>/dev/null
	@rmdir $(shell dirname $(WMAKE_BIN)) 2>/dev/null || true

$(WMAKE_BIN)/dirToString: dirToString.c
	@mkdir -p $(WMAKE_BIN)
	$(call QUIET_MESSAGE,compile,$(<F))
	$E $(cc) $(cFLAGS) $(<F) -o $@

$(WMAKE_BIN)/wmkdep: wmkdep.l
	@mkdir -p $(WMAKE_BIN)
	$(call QUIET_MESSAGE,flex,$(<F))
	$E flex -o $@.c $(<F) && $(cc) $(cFLAGS) $@.c -o $@
	@rm -f $@.c 2>/dev/null

#------------------------------------------------------------------------------
