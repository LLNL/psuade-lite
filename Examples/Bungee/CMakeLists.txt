add_executable (bungee-sim simulator.c)
SET_TARGET_PROPERTIES(bungee-sim
  PROPERTIES OUTPUT_NAME simulator)

if (USE_MARS)
  add_test (Bungee ${PROJECT_SOURCE_DIR}/scripts/psuade_2test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/Bungee/psuade.in ${PROJECT_SOURCE_DIR}/Examples/Bungee/psScript.in ${PROJECT_SOURCE_DIR}/Examples/Bungee/psScript.out)
  install(TARGETS bungee-sim DESTINATION Examples/Bungee)

  if(UNIX)
    install(FILES psuade.in DESTINATION Examples/Bungee)
  else(UNIX)
    install(FILES psuadeForWin.in DESTINATION Examples/Bungee)
  endif(UNIX)


endif (USE_MARS)