#!/bin/bash
laplacianFoam -case ./right -world RIGHT 2>&1 | tee log.run_right
read dummy
