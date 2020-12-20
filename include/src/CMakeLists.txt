project (relate_lib)

add_library(relateStatic STATIC cxxopts data sample anc mutations tree_sequence)
target_link_libraries(relateStatic gzstreamStatic tskitStatic kastoreStatic)
#add_library(relateShared SHARED cxxopts data sample anc mutations tree_sequence)
#target_link_libraries(relateShared gzstreamShared tskitShared kastoreShared)
