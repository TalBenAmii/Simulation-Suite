cmake_minimum_required(VERSION 3.10)
project(simulation-suite VERSION 0.0.1 LANGUAGES C CXX)
include(CTest)

add_executable(branch-predictor bp_main.c bp.cpp bp_api.h)
