cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

#set(CMAKE_PREFIX_PATH $ENV{HOME}/lib CMAKE_PREFIX_PATH)
#find_package(Acts REQUIRED COMPONENTS Core IdentificationPlugin TGeoPlugin DD4hepPlugin PATHS /home/whit/lib/cmake NO_DEFAULT_PATH)
find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )

#-----------------------------------------------------------------------------------
set(a_lib_name topside2)
dd4hep_configure_output()
#dd4hep_package (${a_lib_name} MAJOR 0 MINOR 0 PATCH 1
#    USES         [ROOT   REQUIRED COMPONENTS Geom GenVector] 
#                 [DD4hep REQUIRED COMPONENTS DDCore DDRec]
#    OPTIONAL     XERCESC
#    INCLUDE_DIRS include
#    )

dd4hep_add_plugin(${a_lib_name} SOURCES src/*.cpp
        USES ActsCore ActsPluginDD4hep
  )
target_link_libraries(${a_lib_name}
  PUBLIC DD4hep::DDCore  DD4hep::DDRec
  )

#-----------------------------------------------------------------------------------
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/${a_lib_name}
    FILES_MATCHING PATTERN "*.xml"
    )
