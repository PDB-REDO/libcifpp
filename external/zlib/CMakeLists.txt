cmake_minimum_required(VERSION 3.16)

project(private-static-zlib LANGUAGES C)

include(FetchContent)

FetchContent_Declare(
	zlib
	GIT_REPOSITORY https://github.com/madler/zlib
	GIT_TAG v1.3
)

FetchContent_Populate(zlib)

add_library(zlibstatic STATIC
	${zlib_SOURCE_DIR}/adler32.c
	${zlib_SOURCE_DIR}/crc32.c
	${zlib_SOURCE_DIR}/deflate.c
	${zlib_SOURCE_DIR}/infback.c
	${zlib_SOURCE_DIR}/inffast.c
	${zlib_SOURCE_DIR}/inflate.c
	${zlib_SOURCE_DIR}/inftrees.c
	${zlib_SOURCE_DIR}/trees.c
	${zlib_SOURCE_DIR}/zutil.c
	${zlib_SOURCE_DIR}/compress.c
	${zlib_SOURCE_DIR}/uncompr.c
	${zlib_SOURCE_DIR}/gzclose.c
	${zlib_SOURCE_DIR}/gzlib.c
	${zlib_SOURCE_DIR}/gzread.c
	${zlib_SOURCE_DIR}/gzwrite.c)

add_library(ZLIB::ZLIB ALIAS zlibstatic)
