cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
PROJECT(reference_detector
  VERSION 0.0.1
  LANGUAGES CXX
  DESCRIPTION "A template dd4hep+acts detector"
  )

find_package( DD4hep REQUIRED COMPONENTS DDCore DDG4 )

#find_package(Acts REQUIRED COMPONENTS Core IdentificationPlugin TGeoPlugin DD4hepPlugin PATHS /home/whit/lib/cmake NO_DEFAULT_PATH)
#find_package(Acts REQUIRED COMPONENTS Core PluginIdentification PluginTGeo PluginDD4hep )

#-----------------------------------------------------------------------------------
set(a_lib_name reference_detector)

dd4hep_configure_output()

dd4hep_add_plugin(${a_lib_name} SOURCES src/*.cpp
          USES ActsCore ActsPluginDD4hep
  )
target_link_libraries(${a_lib_name}
  PUBLIC DD4hep::DDCore  DD4hep::DDRec
  )

#-----------------------------------------------------------------------------------
# Install the detector description files.
install(DIRECTORY compact/
    DESTINATION share/${PROJECT_NAME}/${a_lib_name}
    FILES_MATCHING PATTERN "*.xml"
    )
