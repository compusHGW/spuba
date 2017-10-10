CMAKE_MINIMUM_REQUIRED (VERSION 3.0)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

configure_file( jsoncpp-download/CMakeLists.txt.in jsoncpp-download/CMakeLists.txt)
execute_process( COMMAND ${CMAKE_COMMAND} -G ${CMAKE_GENERATOR} . WORKING_DIRECTORY jsoncpp-download)
execute_process( COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY jsoncpp-download )

SET(CMAKE_BUILD_TYPE release)

SET( SOURCES
		main.cpp
        ParticleSource.cpp
		sphere.cpp
	    tube.cpp
	    yield.cpp
	    cone.cpp
		output.cpp
	    bins.cpp
		intersect_objects.cpp
        circle.hpp
		input_json.cpp
		intersection.cpp
		Simulation.cpp Simulation.h
		Cylinder.hpp
		histogram.cpp
		Thruster.hpp
        ParticleSource.h
		bgeom.hpp
		Chamber.hpp
		Geometry.h
		tube.cpp tube.hpp
		sphere.cpp sphere.hpp
		utility.hpp
		particle3d.hpp
		shared.hpp
		const.h
		)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++17 -DBOOST_DISABLE_ASSERT -pthread" )
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

link_directories( depos.x ${JSONCPP_LIBRARY_DIRS} )
ADD_EXECUTABLE ( depos.x ${SOURCES} )


target_link_libraries( depos.x ${CMAKE_BINARY_DIR}/jsoncpp-build/src/lib_json/libjsoncpp.a)
target_include_directories( depos.x PRIVATE ${CMAKE_BINARY_DIR}/jsoncpp-src/include )

