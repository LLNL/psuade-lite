add_executable (mcmc-sim simulator.c)
SET_TARGET_PROPERTIES(mcmc-sim
  PROPERTIES OUTPUT_NAME simulator)
file(COPY ${PROJECT_SOURCE_DIR}/Examples/MCMCTest/mcmcFile DESTINATION ${PROJECT_BINARY_DIR}/Examples/MCMCTest/ )


if (USE_MARS)
  add_test (MCMCTest ${PROJECT_SOURCE_DIR}/scripts/psuade_2test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/MCMCTest/psuadeRS.in ${PROJECT_SOURCE_DIR}/Examples/MCMCTest/psScript.in ${PROJECT_SOURCE_DIR}/Examples/MCMCTest/psScript.out)

  install(TARGETS mcmc-sim DESTINATION Examples/MCMCTest)
  install(FILES mcmcFile DESTINATION Examples/MCMCTest)
  if(UNIX)
    install(FILES psuadeRS.in DESTINATION Examples/MCMCTest)
  else(UNIX)
    install(FILES psuadeRSForWin.in DESTINATION Examples/MCMCTest)
  endif(UNIX)

endif (USE_MARS)