#set the minimum cmake version required
CMAKE_MINIMUM_REQUIRED (VERSION 2.8.11)

#set the name of the project
PROJECT (mandelbrot)

#####################################################################
#enable warning for MSVC and gcc

if (MSVC)
	# Force to always compile with W4
	if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
		string (REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	else ()
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall")
	endif ()
elseif (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

#####################################################################
#enable c++ 11 for MSVC and gcc

include (CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG ("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG ("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else ()
	message (STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif ()

#####################################################################
#enable openmp for MSVC and gcc

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#####################################################################
set (MANDELBROT_MAIN_SOURCES
	main.c
	mandelbrot.c
	mandelbrot.h
	mutil.c
	mutil.h
	mymacro.h
	nrdef.h
	nrutil.c
	nrutil.h
	pi.c
	pi.h
	simd_macro.h
	vnrdef.h
	vnrutil.c
	vnrutil.h
)

add_executable (mandelbrot ${MANDELBROT_MAIN_SOURCES})

# if(MSVC)
	# set (CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
# endif()