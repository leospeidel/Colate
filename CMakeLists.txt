project                     (Relate)
cmake_minimum_required      (VERSION 3.1)

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
endif()
set (CMAKE_CXX_STANDARD 11)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
  	message(STATUS "zlib Found: ${ZLIB_VERSION}")
	  message(STATUS "zlib Include: ${ZLIB_INCLUDE_DIR}")
	  message(STATUS "zlib Libraries: ${ZLIB_LIBRARIES}")
    include_directories( ${ZLIB_INCLUDE_DIRS} ) 
endif( ZLIB_FOUND )

set( OPENSSL_ROOT_DIR "/usr/local/opt/openssl/")
find_package(OpenSSL REQUIRED)
if ( OPENSSL_FOUND )
    message(STATUS "OpenSSL Found: ${OPENSSL_VERSION}")
    message(STATUS "OpenSSL Include: ${OPENSSL_INCLUDE_DIR}")
    message(STATUS "OpenSSL Libraries: ${OPENSSL_LIBRARIES}")
		include_directories( ${OPENSSL_INCLUDE_DIRS} )
endif()

#find_package( OpenMP REQUIRED )
#if ( OpenMP_FOUND )
#  include_directories( ${OpenMP_INCLUDE_DIRS} ) 
#endif( OpenMP_FOUND )

# Set dir
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(HTSLIB_DIR ${CMAKE_SOURCE_DIR}/include/vcf)
set(ADNA_DIR ${CMAKE_SOURCE_DIR}/include/coal)
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

add_subdirectory(${HTSLIB_DIR})
add_subdirectory(${ADNA_DIR})

# tests
add_executable(Tests ${TEST_DIR}/Tests.cpp ${ADNA_DIR}/coal_EM.hpp ${ADNA_DIR}/coal_EM.cpp ${ADNA_DIR}/coal_EM_old.hpp ${ADNA_DIR}/coal_EM_old.cpp)
target_link_libraries(Tests ${RELATE_LIBS} UnitTest)

enable_testing()
add_test(UnitTest ${CMAKE_BINARY_DIR}/Tests)
add_custom_command(TARGET Tests COMMENT "Run tests" POST_BUILD WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND Tests)
