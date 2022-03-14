# CMake 3.8 is required for C++ 17 support
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(DYJets)

# Use FindROOT from the repo, not the system one -- it doesn't have all the features we need
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Require C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
find_package(Boost 1.41.0 REQUIRED COMPONENTS filesystem iostreams program_options system)
find_package(ROOT 6.12 REQUIRED COMPONENTS Hist Physics TreePlayer Unfold)
find_package(Threads REQUIRED)
FIND_PACKAGE(yaml-cpp 0.5 REQUIRED)

# Directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

###################################################################################################
# RooUnfold
###################################################################################################

# Generate dictionary
#include(GNUInstallDirs) # Otherwise ROOT_GENERATE_DICTIONARY fails
#include(${ROOT_USE_FILE})
#ROOT_GENERATE_DICTIONARY(
#    "${CMAKE_BINARY_DIR}/CMakeFiles/RooUnfoldDict"
#    "RooUnfold/src/*.h"
#    LINKDEF "RooUnfold/src/RooUnfold_LinkDef.h"
#)

file(GLOB RooUnfold_src "RooUnfold/src/*.cxx")
add_library(RooUnfold ${RooUnfold_src})  #"${CMAKE_BINARY_DIR}/CMakeFiles/RooUnfoldDict.cxx")

target_include_directories(RooUnfold PUBLIC "RooUnfold/src")
target_include_directories(RooUnfold SYSTEM PUBLIC ROOT::Hist ROOT::Unfold)

###################################################################################################
# Python tools
###################################################################################################

file(INSTALL Main/send-to-cluster DESTINATION bin USE_SOURCE_PERMISSIONS)
file(INSTALL Main/task-check      DESTINATION bin USE_SOURCE_PERMISSIONS)

###################################################################################################
# Configuration
###################################################################################################

configure_file(Includes/cmake_config.h.in Includes/cmake_config.h)

###################################################################################################
# Utilities library
###################################################################################################

add_library(DYJets-util
    Sources/chains.cc
    Sources/event_counter.cc
    Sources/functions.cc # Legacy
    Sources/histo_set.cc
    Sources/histo_set2D.cc
    Sources/job.cc
    Sources/logging.cc
    Sources/options.cc
    Sources/style_list.cc
    Sources/tables.cc
    Sources/timer.cc
)

target_compile_definitions(DYJets-util PUBLIC DYJETS_NEW_API)

# Required headers
target_include_directories(DYJets-util PUBLIC Includes ${CMAKE_BINARY_DIR}/Includes)

# Dependencies
target_link_libraries(DYJets-util PUBLIC Boost::filesystem)
target_link_libraries(DYJets-util PUBLIC Boost::program_options)
target_link_libraries(DYJets-util PUBLIC Boost::system)
target_link_libraries(DYJets-util PUBLIC ROOT::Physics)
target_link_libraries(DYJets-util PUBLIC ROOT::TreePlayer)
target_link_libraries(DYJets-util PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(DYJets-util PUBLIC yaml-cpp)

###################################################################################################
# Data library
###################################################################################################

add_library(DYJets-data
    Sources/catalog.cc
    Sources/mc_group.cc
    Sources/sample.cc
)

# Preprocessor macro to disable incompatible parts of the old code
target_compile_definitions(DYJets-data PUBLIC DYJETS_NEW_API)

# Required headers
target_include_directories(DYJets-data PUBLIC Includes)

# Dependencies
target_link_libraries(DYJets-data PUBLIC DYJets-util)
target_link_libraries(DYJets-data PUBLIC ROOT::Hist)

###################################################################################################
# Physics library
###################################################################################################

add_library(DYJets-physics
    Sources/BTagCalibrationStandalone.cc
    Sources/boson_jets_analyzer.cc
    Sources/btagger.cc
    Sources/electrons.cc
    Sources/genleps.cc
    Sources/jets.cc
    Sources/JetResolution.cc
    Sources/JetResolutionObject.cc
    Sources/muons.cc
    Sources/pileup.cc
    Sources/RoccoR.cc
    Sources/standalone_LumiReWeighting.cc
    Sources/triggers.cc
    Sources/reweighing.cc
    Sources/weights.cc
    Sources/zfinder.cc

)

# Preprocessor macro to disable incompatible parts of the old code
target_compile_definitions(DYJets-physics PUBLIC DYJETS_NEW_API)

# Dependencies
target_link_libraries(DYJets-physics PUBLIC DYJets-util)

###################################################################################################
# Comparison library
###################################################################################################

add_library(DYJets-compare
    Sources/comparison_entry.cc
    Sources/compare_builder_base.cc
    Sources/reco_compare_builder.cc
    Sources/signal_bg_compare_builder.cc
    Sources/ss_yield_compare_builder.cc
)

# Preprocessor macro to disable incompatible parts of the old code
target_compile_definitions(DYJets-compare PUBLIC DYJETS_NEW_API)

# Dependencies
target_link_libraries(DYJets-compare PUBLIC DYJets-data)

###################################################################################################
# Higgs targets
###################################################################################################

add_executable(higgs-compare Main/higgs-compare.cc)
target_link_libraries(higgs-compare PRIVATE DYJets-compare)

add_executable(higgs-job-list Main/higgs-job-list.cc)
target_link_libraries(higgs-job-list PRIVATE DYJets-data)

add_executable(higgs-loop Main/higgs-loop.cc Sources/higgs_analyzer.cc)
target_link_libraries(higgs-loop PRIVATE DYJets-data)
target_link_libraries(higgs-loop PRIVATE DYJets-physics)

###################################################################################################
# Drell-Yan targets
###################################################################################################

add_executable(dyjets-compare Main/dyjets-compare.cc)
target_link_libraries(dyjets-compare PRIVATE DYJets-compare)

add_executable(dyjets-matrix Main/dyjets-matrix.cc)
target_link_libraries(dyjets-matrix PRIVATE DYJets-compare)

add_executable(dyjets-job-list Main/dyjets-job-list.cc)
target_link_libraries(dyjets-job-list PRIVATE DYJets-data)

add_executable(dyjets-loop Main/dyjets-loop.cc Sources/dyjets_analyzer.cc)
target_link_libraries(dyjets-loop PRIVATE DYJets-data)
target_link_libraries(dyjets-loop PRIVATE DYJets-physics)

#add_executable(dyjets-sample-weights Main/dyjets-sample-weights.cc)
#target_link_libraries(dyjets-sample-weights PRIVATE DYJets-data)
#target_link_libraries(dyjets-sample-weights PRIVATE DYJets-physics)

###################################################################################################
# Legacy library
###################################################################################################

#add_library(DYJets-legacy
#    Sources/ArgParser.cc
#    Sources/BLUEMeth.cc
#    Sources/BTagCalibrationStandalone.cc
#    Sources/Combination.cc
#    Sources/Config.cc
#    Sources/ConfigVJets.cc
#    Sources/fixYscale.cc
#    Sources/functions.cc
#    Sources/getFilesAndHistogramsZJets.cc
#    Sources/HistoSetZJets.cc
#    Sources/JetCorrectionUncertainty.cc
#    Sources/JetCorrectorParameters.cc
#    Sources/JetResolution.cc
#    Sources/JetResolutionObject.cc
#    Sources/printTable.cc
#    Sources/PlotSettings.cc
#    Sources/RecoComparison.cc
#    Sources/ResultTables.cc
#    Sources/RoccoR.cc
#    Sources/rootFunctions.cc
#    Sources/SectionedConfig.cc
#    Sources/SimpleJetCorrectionUncertainty.cc
#    Sources/TextTools.cc
#    Sources/Uncertainties.cc
#    Sources/variablesOfInterestZJets.cc

#    # Some sources shared with the new code
#    Sources/catalog.cc
#    Sources/logging.cc
#    Sources/standalone_LumiReWeighting.cc
#    Sources/tables.cc
#)

## Required headers
#target_include_directories(DYJets-legacy PUBLIC Includes ${CMAKE_BINARY_DIR}/Includes)
#target_include_directories(DYJets-legacy SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
#target_include_directories(DYJets-legacy SYSTEM PUBLIC ${ROOT_INCLUDE_DIRS})

## Dependencies
#target_link_libraries(DYJets-legacy PUBLIC ${ROOT_LIBRARIES})
#target_link_libraries(DYJets-legacy PUBLIC ${CMAKE_THREAD_LIBS_INIT})
#target_link_libraries(DYJets-legacy PUBLIC RooUnfold)
#target_link_libraries(DYJets-legacy PUBLIC yaml-cpp)

####################################################################################################
## Legacy targets
####################################################################################################

#add_executable(plotCrossSections Main/plotCrossSections.cc)
#target_link_libraries(plotCrossSections PRIVATE DYJets-legacy)

#add_executable(plotRespMatrix Sources/plotRespMatrix.cc)
#target_link_libraries(plotRespMatrix PRIVATE DYJets-legacy)

#add_executable(runIntegralCheck Main/runIntegralCheck.cc)
#target_link_libraries(runIntegralCheck PRIVATE DYJets-legacy)

#add_executable(runRecoComparison Main/runRecoComparison.cc)
#target_link_libraries(runRecoComparison PRIVATE DYJets-legacy)

#add_executable(runUnfoldingZJets Main/runUnfoldingZJets.cc Sources/UnfoldingZJets.cc)
#target_link_libraries(runUnfoldingZJets PRIVATE DYJets-legacy)

#add_executable(runZJets_newformat Main/runZJets_newformat.cc Sources/ZJets_newformat.cc)
#target_compile_options(runZJets_newformat PRIVATE -fno-var-tracking-assignments) # Unfortunate :'(
#target_link_libraries(runZJets_newformat PRIVATE DYJets-legacy)

#add_executable(testDecompVar Sources/testDecompVar.cc)
#target_link_libraries(testDecompVar PRIVATE DYJets-legacy)

