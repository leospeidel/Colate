project                     (Relate)
cmake_minimum_required      (VERSION 3.1)

if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
endif()
set (CMAKE_CXX_STANDARD 11)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} ) 
endif( ZLIB_FOUND )

#find_package( OpenMP REQUIRED )
#if ( OpenMP_FOUND )
#  include_directories( ${OpenMP_INCLUDE_DIRS} ) 
#endif( OpenMP_FOUND )

# Set dir
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(HTSLIB_DIR ${CMAKE_SOURCE_DIR}/include/vcf)
set(ADNA_DIR ${CMAKE_SOURCE_DIR}/include/coal)
set(FTREE_DIR ${CMAKE_SOURCE_DIR}/include/ftree)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/include/test)

set(RELATE_LIBS relateShared gzstreamShared kastoreShared tskitShared)
link_directories(${CMAKE_SOURCE_DIR}/../relate_lib/bin)
link_directories(${HTSLIB_DIR})
set(HTSLIB hts)

# Add subdirectories
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/../relate_lib/include/src)
include_directories(${CMAKE_SOURCE_DIR}/../relate_lib/include/src/tskit)
include_directories(${CMAKE_SOURCE_DIR}/../relate_lib/include/src/gzstream)
include_directories(${CMAKE_SOURCE_DIR}/include/vcf)
include_directories(${CMAKE_SOURCE_DIR}/include/coal)
include_directories(${CMAKE_SOURCE_DIR}/include/ftree)

add_subdirectory(${HTSLIB_DIR})
add_subdirectory(${ADNA_DIR})
add_subdirectory(${FTREE_DIR})

# tests
add_executable(Tests ${TEST_DIR}/Tests.cpp ${ADNA_DIR}/coal_EM.hpp ${ADNA_DIR}/coal_EM.cpp)
target_link_libraries(Tests ${RELATE_LIBS} UnitTest)

enable_testing()
add_test(UnitTest ${CMAKE_BINARY_DIR}/Tests)
add_custom_command(TARGET Tests COMMENT "Run tests" POST_BUILD WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND Tests)
