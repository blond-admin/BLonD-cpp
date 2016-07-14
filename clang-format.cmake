# additional target to perform clang-format run, requires clang-format
# get all project files

file(GLOB_RECURSE ALL_SOURCE_FILES *.cpp *.h)
foreach (SOURCE_FILE ${ALL_SOURCE_FILES})
    string(FIND ${SOURCE_FILE} ${PROJECT_TRDPARTY_DIR} PROJECT_TRDPARTY_DIR_FOUND)
    if (NOT ${PROJECT_TRDPARTY_DIR_FOUND} EQUAL -1)
        list(REMOVE_ITEM ALL_SOURCE_FILES ${SOURCE_FILE})
    endif()
endforeach()

set(CLANG_FMT_CMD "clang-format" )
add_custom_target( clangformat
        COMMAND ${CLANG_FMT_CMD}
        -style=LLVM
        -style="{BasedOnStyle: LLVM,
                IndentWidth: 4,
                Language: Cpp,
                DerivePointerAlignment: false,
                PointerAlignment: Left,
                AllowShortIfStatementsOnASingleLine: false,
                AllowShortLoopsOnASingleLine: false,
                NamespaceIndentation: All
               }"
        -i ${ALL_SOURCE_FILES}
        )
