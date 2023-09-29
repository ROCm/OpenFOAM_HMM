find build/*/applications build/*/src -name "*.o" | xargs rm
find build/*/applications build/*/src -name "*.dep" | xargs rm
find build/*/applications build/*/src -name "*.so" | xargs rm
find platforms/linux64ClangDPInt32Opt/lib/ -name "*.so" | xargs rm

