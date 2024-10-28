add_executable(Colate Colate.cpp coal_tree.hpp coal_tree.cpp coal_EM.hpp coal_EM.cpp coal_EM_old.hpp coal_EM_old.cpp)
target_link_libraries(Colate ${HTSLIB})

add_executable(CoalRate CoalRate.cpp coal_tree.hpp coal_tree.cpp coal_EM.hpp coal_EM.cpp coal_EM_old.hpp coal_EM_old.cpp)
target_link_libraries(CoalRate ${HTSLIB})
