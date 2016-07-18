# additional target to perform clang-format run, requires clang-format
# get all project files

if(${WITH_FORMAT})

file(GLOB_RECURSE ALL_SOURCE_FILES src/*.cpp include/*.h)

find_program (CLANG_FMT_CMD
        NAMES "clang-format"
        )
add_custom_target( clangformat
        COMMAND ${CLANG_FMT_CMD} -style="{BasedOnStyle: LLVM, IndentWidth: 4, Language: Cpp, DerivePointerAlignment: false, PointerAlignment: Left, AllowShortIfStatementsOnASingleLine: false, AllowShortLoopsOnASingleLine: false, NamespaceIndentation: All}" -i ${ALL_SOURCE_FILES}
        )
endif()