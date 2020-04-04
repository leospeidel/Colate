add_executable(RelateaDNA RelateaDNA.cpp aDNA_EM.hpp aDNA_EM.cpp)
#target_link_libraries(RelateaDNA ${RELATE_LIBS} OpenMP::OpenMP_CXX)
target_link_libraries(RelateaDNA ${RELATE_LIBS} /usr/local/Cellar/gperftools/2.7/lib/libprofiler.dylib)
