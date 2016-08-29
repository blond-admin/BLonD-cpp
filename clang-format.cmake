# additional target to perform clang-format run, requires clang-format
# get all project files

if(${WITH_FORMAT})
    file(GLOB_RECURSE ALL_SOURCE_FILES
            src/*.cpp
            include/*.h
            demos/*.cpp
            unit-tests/*.cpp
            benchmarks/*.cpp)

    find_program (CLANG_FMT_CMD
            NAMES "clang-format"
            )
    if(NOT CLANG_FMT_CMD)
        message(ERROR "clang-format not found!")
    else()
        add_custom_target( clangformat
            COMMAND ${CLANG_FMT_CMD} -style="{BasedOnStyle: LLVM, IndentWidth: 4, Language: Cpp, DerivePointerAlignment: false, PointerAlignment: Left, AllowShortIfStatementsOnASingleLine: false, AllowShortLoopsOnASingleLine: false, NamespaceIndentation: All}" -i ${ALL_SOURCE_FILES}
                )
    endif()
endif()