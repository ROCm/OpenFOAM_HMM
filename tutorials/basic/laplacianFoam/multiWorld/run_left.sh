#!/bin/bash
laplacianFoam -case left -world LEFT 2>&1 | tee log.run_left
read dummy
