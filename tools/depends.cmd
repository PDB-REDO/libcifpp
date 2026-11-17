@ECHO OFF

@REM  SET ZLIB_VERSION=1.3
@REM  SET LZMA_VERSION=5.6.3

@REM  IF NOT EXIST build_ci\libs (
@REM    MKDIR build_ci\libs
@REM  )
@REM  CD build_ci\libs

@REM  @REM Install ZLib
@REM  IF NOT EXIST zlib-%ZLIB_VERSION%.zip (
@REM    ECHO Downloading https://github.com/libarchive/zlib/archive/v%ZLIB_VERSION%.zip
@REM    curl -L -o zlib-%ZLIB_VERSION%.zip https://github.com/libarchive/zlib/archive/v%ZLIB_VERSION%.zip || EXIT /b 1
@REM  )
@REM  IF NOT EXIST zlib-%ZLIB_VERSION% (
@REM    ECHO Unpacking zlib-%ZLIB_VERSION%.zip
@REM    C:\windows\system32\tar.exe -x -f zlib-%ZLIB_VERSION%.zip || EXIT /b 1
@REM  )
@REM  CD zlib-%ZLIB_VERSION%
@REM  cmake -B build || EXIT /b 1
@REM  cmake --build build --target ALL_BUILD --config Release || EXIT /b 1
@REM  cmake --build build --target RUN_TESTS --config Release || EXIT /b 1
@REM  cmake --build build --target INSTALL --config Release || EXIT /b 1
@REM  CD ..

@REM  @REM Install LibLZMA
@REM  IF NOT EXIST xz-%LZMA_VERSION%.zip (
@REM    ECHO Downloading https://github.com/libarchive/xz/archive/refs/tags/v%LZMA_VERSION%.zip
@REM    curl -L -o xz-%LZMA_VERSION%.zip https://github.com/libarchive/xz/archive/refs/tags/v%LZMA_VERSION%.zip || exit /b 1
@REM  )
@REM  IF NOT EXIST xz-%LZMA_VERSION% (
@REM    ECHO Unpacking xz-%LZMA_VERSION%.zip
@REM    C:\windows\system32\tar.exe -x -f xz-%LZMA_VERSION%.zip || EXIT /b 1
@REM  )
@REM  CD xz-%LZMA_VERSION%
@REM  cmake -B build || EXIT /b 1
@REM  cmake --build build --target ALL_BUILD --config Release || EXIT /b 1
@REM  cmake --build build --target RUN_TESTS --config Release || EXIT /b 1
@REM  cmake --build build --target INSTALL --config Release || EXIT /b 1
@REM  CD ..

vcpkg install

@EXIT /b 0
