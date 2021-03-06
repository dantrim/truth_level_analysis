atlas_subdir( truth_level_analysis )

find_package(lwtnn)

atlas_depends_on_subdirs(
    PUBLIC
    Event/xAOD/xAODBase
    Control/xAODRootAccess
    Event/xAOD/xAODEventInfo
    Event/xAOD/xAODTruth
    Event/xAOD/xAODMissingET
    Event/xAOD/xAODJet
    PhysicsAnalysis/MCTruthClassifier
    Tools/PathResolver
)

find_package( ROOT REQUIRED COMPONENTS RIO Hist Tree Core Gpad Graf Graf3d )

set(_common
    src/*.cxx
    src/utility/*.cxx
    INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} truth_level_analysis truth_level_analysis/utility ${LWTNN_INCLUDE_DIRS}
    LINK_LIBRARIES ${ROOT_LIBRARIES} ${LWTNN_LIBRARIES}
    xAODRootAccess
    xAODJet
    xAODTruth
    xAODEventInfo
    xAODMissingET
    PathResolver
)


function(MakeExec filename)
    set(execname)
    get_filename_component(execname ${filename} NAME_WE)
    atlas_add_executable( ${execname} "util/${execname}.cxx" ${_common} )
endfunction(MakeExec)

file(GLOB files "util/*.cxx")
foreach(file ${files})
    MakeExec(${file})
endforeach()

# locate data/ dir to be picked up by PathResolver
atlas_install_data( data/* )
