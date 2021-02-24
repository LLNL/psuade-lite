add_executable (optrosenbrockm-sim simulator.c)
SET_TARGET_PROPERTIES(optrosenbrockm-sim
  PROPERTIES OUTPUT_NAME simulator)
target_link_libraries(optrosenbrockm-sim m)

if (USE_BOBYQA)
  add_test (OptRosenbrockM ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/OptRosenbrockM/psuadeBobyqa.in ${PROJECT_SOURCE_DIR}/Examples/OptRosenbrockM/psScript.out)

  install(TARGETS optrosenbrockm-sim DESTINATION Examples/OptRosenbrockM)
  if(UNIX)
    install(FILES psuadeBobyqa.in DESTINATION Examples/OptRosenbrockM)
  else(UNIX)
    install(FILES psuadeBobyqaForWin.in DESTINATION Examples/OptRosenbrockM)
  endif(UNIX)

endif (USE_BOBYQA)