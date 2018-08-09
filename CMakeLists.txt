atlas_subdir( truth_level_analysis )

atlas_depends_on_subdirs(
    PUBLIC
    Event/xAOD/xAODBase
    Control/xAODRootAccess
    Event/xAOD/xAODEventInfo
    Event/xAOD/xAODTruth
    Event/xAOD/xAODMissingET
    Event/xAOD/xAODJet
    PhysicsAnalysis/MCTruthClassifier
)

find_package( ROOT REQUIRED COMPONENTS RIO Hist Tree Core Gpad Graf Graf3d )

set(_common
    src/*.cxx
    src/utility/*.cxx
    INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} truth_level_analysis truth_level_analysis/utility
    LINK_LIBRARIES ${ROOT_LIBRARIES}
    xAODRootAccess
    xAODJet
    xAODTruth
    xAODEventInfo
    xAODMissingET
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
