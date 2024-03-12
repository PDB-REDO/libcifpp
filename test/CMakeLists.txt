# We're using the older version 2 of Catch2

find_package(Catch2 QUIET)

if(NOT Catch2_FOUND)
	FetchContent_Declare(
		Catch2
		GIT_REPOSITORY https://github.com/catchorg/Catch2.git
		GIT_TAG v2.13.9)

	FetchContent_MakeAvailable(Catch2)

	set(Catch2_VERSION "2.13.9")
endif()

list(
	APPEND
	CIFPP_tests
	unit-v2
	unit-3d
	format
	model
	rename-compound
	sugar
	spinner
	# reconstruction
	validate-pdbx)

add_library(test-main OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/test-main.cpp")

target_link_libraries(test-main cifpp::cifpp Catch2::Catch2)

if(${Catch2_VERSION} VERSION_GREATER_EQUAL 3.0.0)
	target_compile_definitions(test-main PUBLIC CATCH22=0)
else()
	target_compile_definitions(test-main PUBLIC CATCH22=1)
endif()

foreach(CIFPP_TEST IN LISTS CIFPP_tests)
	set(CIFPP_TEST "${CIFPP_TEST}-test")
	set(CIFPP_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${CIFPP_TEST}.cpp")

	add_executable(
		${CIFPP_TEST} ${CIFPP_TEST_SOURCE} $<TARGET_OBJECTS:test-main>)

	if(${Catch2_VERSION} VERSION_GREATER_EQUAL 3.0.0)
		target_compile_definitions(${CIFPP_TEST} PUBLIC CATCH22=0)
	else()
		target_compile_definitions(${CIFPP_TEST} PUBLIC CATCH22=1)
	endif()

	target_link_libraries(${CIFPP_TEST} PRIVATE Threads::Threads cifpp::cifpp
		Catch2::Catch2)
	target_include_directories(${CIFPP_TEST} PRIVATE "${EIGEN_INCLUDE_DIR}")

	if(MSVC)
		# Specify unwind semantics so that MSVC knowns how to handle exceptions
		target_compile_options(${CIFPP_TEST} PRIVATE /EHsc)
	endif()

	add_custom_target(
		"run-${CIFPP_TEST}"
		DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch ${CIFPP_TEST})

	add_custom_command(
		OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch
		COMMAND $<TARGET_FILE:${CIFPP_TEST}> --data-dir
		${CMAKE_CURRENT_SOURCE_DIR})

	add_test(NAME ${CIFPP_TEST} COMMAND $<TARGET_FILE:${CIFPP_TEST}> --data-dir
		${CMAKE_CURRENT_SOURCE_DIR})
endforeach()