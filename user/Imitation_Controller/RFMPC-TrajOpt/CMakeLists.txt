include_directories(SYSTEM "../../../third-party/qpOASES/include")
include_directories(SYSTEM "../../../third-party/qpSWIFT/include")

file(GLOB_RECURSE sources "*.cpp")
add_library(rfmpc-trajopt SHARED ${sources})
target_link_libraries(rfmpc-trajopt qpSWIFT-static qpOASES)



