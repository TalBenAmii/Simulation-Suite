cmake_minimum_required(VERSION 3.10)
project(simulation-suite VERSION 0.0.1 LANGUAGES C)
include(CTest)

add_executable(raw-dependency-analyzer dflow_calc.c dflow_calc.h dflow_main.c)
