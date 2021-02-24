add_executable (morris20-sim simulator.c)
SET_TARGET_PROPERTIES(morris20-sim
  PROPERTIES OUTPUT_NAME simulator)

install(TARGETS morris20-sim DESTINATION Examples/Morris20)

add_test (Morris20MOAT ${PROJECT_SOURCE_DIR}/scripts/psuade_2test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/Morris20/psuadeMOAT.in ${PROJECT_SOURCE_DIR}/Examples/Morris20/psMOATScript.in ${PROJECT_SOURCE_DIR}/Examples/Morris20/psMOATScript.out)

if(UNIX)
  install(FILES psuadeMOAT.in DESTINATION Examples/Morris20)
else(UNIX)
  install(FILES psuadeMOATForWin.in DESTINATION Examples/Morris20)
endif(UNIX)

if (USE_MARS)
add_test (Morris20LH ${PROJECT_SOURCE_DIR}/scripts/psuade_2test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/Morris20/psuadeLH.in ${PROJECT_SOURCE_DIR}/Examples/Morris20/psLHScript.in ${PROJECT_SOURCE_DIR}/Examples/Morris20/psLHScript.out)

  if(UNIX)
    install(FILES psuadeLH.in DESTINATION Examples/Morris20)
  else(UNIX)
    install(FILES psuadeLHForWin.in DESTINATION Examples/Morris20)
  endif(UNIX)

endif(USE_MARS)