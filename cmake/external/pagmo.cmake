if (NOT __PAGMO_INCLUDED) # guard against multiple includes
    set(__PAGMO_INCLUDED TRUE)

    # use the system-wide smart-math if present
    find_package(pagmo)
    if (PAGMO_FOUND)
        set(PAGMO_EXTERNAL FALSE)
    else ()
        set(PAGMO_FOUND FALSE)
        set(PAGMO_EXTERNAL TRUE)
    endif ()

endif ()