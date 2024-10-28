include_directories(${PROJECT_SOURCE_DIR})
include_directories("../../robot/include")
include_directories("../../common/include/")
include_directories("../../common/include/Utilities")
include_directories("../../lcm-types/cpp")

add_executable(follow_trajectory_ctrl Follow_Trajectory_Controller.cpp main.cpp)
target_link_libraries(follow_trajectory_ctrl robot biomimetics dynacore_param_handler)