function(GetnerateLibName arg)
    set(arg "${CMAKE_STATIC_LIBRARY_PREFIX}${arg}${CMAKE_STATIC_LIBRARY_SUFFIX}" PARENT_SCOPE)
endfunction()