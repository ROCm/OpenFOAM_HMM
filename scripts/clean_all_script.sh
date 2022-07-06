find build/*/applications build/*/src -name "*.o" | xargs rm
find build/*/applications build/*/src -name "*.dep" | xargs rm
find build/*/applications build/*/src -name "*.so" | xargs rm
