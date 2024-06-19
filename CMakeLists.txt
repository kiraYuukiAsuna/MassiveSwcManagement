cmake_minimum_required(VERSION 3.28)
project(MassiveSwcManagement)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 23)

add_executable(MassiveSwcManagement
        main.cpp
)

target_sources(MassiveSwcManagement PRIVATE
        FILE_SET CXX_MODULES FILES
        Source/STLModules/std.ixx
        Source/STLModules/std.compat.ixx
        Source/ZOrderTree/BitsetUtil.cppm
        Source/ZOrderTree/ZOrderTree.cppm
        Source/ZOrderTree/TypeAlias.cppm
)

