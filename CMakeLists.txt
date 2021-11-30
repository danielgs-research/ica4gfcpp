cmake_minimum_required(VERSION 3.16)
project(ICALIB)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package (Threads REQUIRED)

find_library(GMP gmp)



file(GLOB NTL_HEADER ABSOLUTE "vendor/ntl/include/*.h")
# file(GLOB NTL_SRC ABSOLUTE "vendor/ntl/src/*.cpp")
get_filename_component(NTL_SRC "libntl.a" ABSOLUTE BASE_DIR /usr/local/lib)
set(ENV{NTL_INCLUDE} "/usr/local/include/NTL")

# find_library(NTL 
# NAMES libntl ntl ntllib 
# HINTS /usr/local)
add_executable(ICALIB test.cpp )

message(${NTL_SRC})
target_link_libraries(ICALIB ${NTL_SRC} Threads::Threads gmp) #libgmp)

# message("NTL INCLUDES")
# message(${NTL_SRC})
# add_library(NTL ${NTL_HEADERS})
# add_library(NTL ${NTL_SRC} ${NTL_HEADER})
# target_include_directories(ICALIB PRIVATE $ENV{NTL_INCLUDE} )
#target_link_libraries(ICALIB )
#target_link_libraries(ICALIB PUBLIC ${NTL_SRC} )
# target_include_directories(ICALIB PRIVATE $ENV{NTL_INCLUDE})