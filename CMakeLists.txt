project                     (Colate)
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

# Search OpenSSL

find_package(OpenSSL REQUIRED) 	
if ( OPENSSL_FOUND )
if(APPLE)
	include_directories(/usr/local/Cellar/openssl@1.1/1.1.1i/include)
	list(APPEND LIB_LIST /usr/local/Cellar/openssl@1.1/1.1.1i/lib/libssl.dylib)
	list(APPEND LIB_LIST /usr/local/Cellar/openssl@1.1/1.1.1i/lib/libcrypto.dylib)
	message(STATUS "OpenSSL Version: ${OPENSSL_VERSION} ${OPENSSL_INCLUDE_DIR} ${OPENSSL_LIBRARIES}")
else()
	include_directories(${OPENSSL_INCLUDE_DIR})
	list(APPEND LIB_LIST ${OPENSSL_LIBRARIES})
	message(STATUS "OpenSSL Version: ${OPENSSL_VERSION} ${OPENSSL_INCLUDE_DIR} ${OPENSSL_LIBRARIES}")
endif()
endif()

add_definitions(-static)
add_definitions(-static-libgcc)
set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static -static-libgcc")

# Set dir
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/src)
set(GZSTREAM_DIR ${CMAKE_SOURCE_DIR}/include/src/gzstream)
set(TSKIT_DIR ${CMAKE_SOURCE_DIR}/include/src/tskit)
set(HTSLIB_DIR ${CMAKE_SOURCE_DIR}/include/vcf)
set(COLATE_DIR ${CMAKE_SOURCE_DIR}/include/coal)
set(TEST_DIR ${CMAKE_SOURCE_DIR}/include/test)

link_directories(${HTSLIB_DIR})
set(HTSLIB hts)


# Add subdirectories
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${GZSTREAM_DIR})
include_directories(${TEST_DIR})
include_directories(${TSKIT_DIR})
include_directories(${HTSLIB_DIR})
include_directories(${COLATE_DIR})

add_subdirectory(${PROJECT_SOURCE_DIR})
add_subdirectory(${GZSTREAM_DIR})
add_subdirectory(${TSKIT_DIR})
add_subdirectory(${TEST_DIR})
add_subdirectory(${HTSLIB_DIR})
add_subdirectory(${COLATE_DIR})

# tests
add_executable(Tests ${TEST_DIR}/Tests.cpp ${COLATE_DIR}/coal_EM.hpp ${COLATE_DIR}/coal_EM.cpp ${COLATE_DIR}/coal_EM_old.hpp ${COLATE_DIR}/coal_EM_old.cpp)
target_link_libraries(Tests relateStatic UnitTest)

enable_testing()
add_test(UnitTest ${CMAKE_BINARY_DIR}/Tests)
add_custom_command(TARGET Tests COMMENT "Run tests" POST_BUILD WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND Tests)
