blockMesh

setSet -batch createObstacle.setSet

subsetMesh c0 -patch walls

#mv 0 0-orig
mv 0.001 0

mv 0/polyMesh/* constant/polyMesh/
rmdir 0/polyMesh

setFields
