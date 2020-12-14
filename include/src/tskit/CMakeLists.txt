project (tskit)

add_library(kastoreStatic STATIC kastore)
add_library(kastoreShared SHARED kastore)

add_library(tskitStatic STATIC kastore core tables trees genotypes convert stats tskit)
add_library(tskitShared SHARED kastore core tables trees genotypes convert stats tskit)

