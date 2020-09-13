find_path(LEMON_DIRECTORY
          NAMES include/lemon/concepts/graph.h
          PATHS "C:/Program Files (x86)/LEMON" "C:/Program Files/LEMON" "$ENV{HOME}/lemon-1.3.1/build"
         )

MESSAGE("LEMON: ${LEMON_DIRECTORY}")

if(LEMON_DIRECTORY)

    set(LEMON_FOUND TRUE)

    set(LEMON_INCLUDE_DIRECTORIES "${LEMON_DIRECTORY}/include")
    set(LEMON_LINK_DIRECTORIES "${LEMON_DIRECTORY}/lib")

    if(WIN32)
        if (CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(LEMON_LINK_LIBRARIES "${LEMON_DIRECTORY}/lib/lemon_debug.lib")
        else()
            set(LEMON_LINK_LIBRARIES "${LEMON_DIRECTORY}/lib/lemon.lib")
        endif()
    elseif(APPLE)
        set(LEMON_LINK_LIBRARIES "${LEMON_DIRECTORY}/lemon/libemon.a")
    else()
        set(LEMON_LINK_LIBRARIES "")
    endif()

else()
    set(LEMON_FOUND FALSE)
    MESSAGE("LEMON: directory not found, but may still work")
endif()