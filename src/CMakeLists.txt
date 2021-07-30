#CMake file for Adding executables

function( ADD_RBM_TESTCASE arg1)

    add_executable(${arg1} ${arg1}.cpp)
    target_link_libraries(${arg1} MODES ${MANDATORY_LIBRARIES})

endfunction()

#ADD_RBM_TESTCASE(AdaptiveMODES_Offline)
#ADD_RBM_TESTCASE(Adaptive_MODESConservative)
#ADD_RBM_TESTCASE(AdaptiveMODES_Online)
ADD_RBM_TESTCASE(AdaptiveMODES_DTR_Offline)
ADD_RBM_TESTCASE(AdaptiveMODES_DTR_Online)
ADD_RBM_TESTCASE(Nmodes_Sensitivity_DTR)
ADD_RBM_TESTCASE(Nmodes_Sensitivity_Direct)
#ADD_RBM_TESTCASE(SingleMODES)
ADD_RBM_TESTCASE(POD_Rec)
#ADD_RBM_TESTCASE(ExtractPointEvoTime)


if(ENABLE_PAGMO)
    ADD_RBM_TESTCASE(Adaptive_Sampling)
    ADD_RBM_TESTCASE(AdaptSample_ResEval)
endif()

