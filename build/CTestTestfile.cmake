# CMake generated Testfile for 
# Source directory: /data/smew1/speidel/genomics/relate_aDNA
# Build directory: /data/smew1/speidel/genomics/relate_aDNA/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTest "/data/smew1/speidel/genomics/relate_aDNA/bin/Tests")
set_tests_properties(UnitTest PROPERTIES  _BACKTRACE_TRIPLES "/data/smew1/speidel/genomics/relate_aDNA/CMakeLists.txt;39;add_test;/data/smew1/speidel/genomics/relate_aDNA/CMakeLists.txt;0;")
subdirs("include/aDNA")
