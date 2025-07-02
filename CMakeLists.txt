cmake_minimum_required(VERSION 3.30)
project(Goldwasser_Micali_Multithread)

set(CMAKE_CXX_STANDARD 20)


add_executable(Goldwasser_Micali_Multithread test.cpp
        Goldwasser_Micali_Multithread.cpp
        Goldwasser_Micali_Multithread.h
        base64.cpp
        base64.h
)

set_target_properties(Goldwasser_Micali_Multithread  PROPERTIES OUTPUT_NAME "test")