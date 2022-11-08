cmake_minimum_required(VERSION 3.16..3.19)

function(add_git_submodule dir)
	# add a Git submodule directory to CMake, assuming the
	# Git submodule directory is a CMake project.
	#
	# Usage: in CMakeLists.txt
	#
	# include(AddGitSubmodule.cmake)
	# add_git_submodule(mysubmod_dir)
	find_package(Git QUIET)

	if(NOT EXISTS "${PROJECT_SOURCE_DIR}/${dir}/CMakeLists.txt")
		if(NOT (GIT_FOUND AND IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/.git"))
			message(FATAL_ERROR "${CMAKE_CURRENT_SOURCE_DIR} is not a git repository and the submodule ${dir} is not complete. Cannot continue.")
		else()
			if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
				execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- ${dir}
					WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
					COMMAND_ERROR_IS_FATAL ANY)
			else()
				execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -- ${dir}
					WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
			endif()
		endif()
	endif()

	set(ENABLE_TESTING OFF)

	add_subdirectory(${dir} ${ARGV})
endfunction(add_git_submodule)