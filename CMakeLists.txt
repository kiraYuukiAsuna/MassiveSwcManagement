cmake_minimum_required(VERSION 3.28)
project(MassiveSwcManagement)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 23)

add_executable(MassiveSwcManagement
        main.cpp
        Source/MOctree/BitsetUtil.cppm
        Source/MOctree/MOctree.cppm
)

target_sources(MassiveSwcManagement PRIVATE
        FILE_SET CXX_MODULES FILES
        Source/STLModules/std.ixx
        Source/STLModules/std.compat.ixx
        Source/MOctree/BitsetUtil.cppm
        Source/MOctree/MOctree.cppm
)

