option(TEST_SOLUTION "Build solution" OFF)

function(add_hse_executable NAME)
    add_executable(${NAME} ${ARGN})
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endfunction()

function(add_catch TARGET)
    add_hse_executable(${TARGET} ${ARGN})
    target_link_libraries(${TARGET} contrib_catch_main)

    if (TEST_SOLUTION)
        add_custom_target(
                run_${TARGET}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                DEPENDS ${TARGET}
                COMMAND ${CMAKE_BINARY_DIR}/${TARGET})

        add_dependencies(test-all run_${TARGET})
    endif ()
endfunction()

add_custom_target(test-all)
