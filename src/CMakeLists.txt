set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# ===== External Dependency: GoogleTest
if (NOT TARGET gtest)
    configure_file(CMakeLists.txt.in.googletest googletest-download/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            RESULT_VARIABLE result
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
            ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
            EXCLUDE_FROM_ALL)
    if (CMAKE_VERSION VERSION_LESS 2.8.11)
        include_directories("${gtest_SOURCE_DIR}/include")
    endif()
endif()

set(TCROSSER_HEADERS include/tcrosser.h)

add_library(tcrosser INTERFACE)
target_include_directories(
        tcrosser
        INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/src/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(tcrosser INTERFACE cxx_std_17)


install(TARGETS tcrosser
        EXPORT tcrosser_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

add_executable(tcrosser_test tcrosser_test.cpp)
add_dependencies(tcrosser_test tcrosser)
target_compile_features(tcrosser_test PRIVATE cxx_std_17)
target_link_libraries(tcrosser_test tcrosser gtest gtest_main)
add_test(NAME tcrosser_test COMMAND tcrosser_test)

if(MSVC)
    set(MY_CXX_DEBUG_OPTIONS ${MY_CXX_DEBUG_OPTIONS} /W4)
    set(MY_CXX_RELEASE_OPTIONS ${MY_CXX_RELEASE_OPTIONS} /W4)
else()
    set(MY_CXX_DEBUG_OPTIONS ${MY_CXX_DEBUG_OPTIONS}  -Wall -Wextra -pedantic)
    set(MY_CXX_RELEASE_OPTIONS ${MY_CXX_RELEASE_OPTIONS} -Wall -Wextra -pedantic)
endif()

target_compile_options(tcrosser_test PUBLIC "$<$<CONFIG:DEBUG>:${MY_CXX_DEBUG_OPTIONS}>")
target_compile_options(tcrosser_test PUBLIC "$<$<CONFIG:RELEASE>:${MY_CXX_RELEASE_OPTIONS}>")

add_executable(tcrosser_main tcrosser_main.cpp)
add_dependencies(tcrosser_main tcrosser)
target_compile_features(tcrosser_main PRIVATE cxx_std_17)
