#/bin/sh
cd ${0%/*} || exit 1    # run from this directory

# this will have to do until we make a makefile rule

Coco \
    -frames $WM_THIRD_PARTY_DIR/coco-r \
    SimpleCalc.atg

