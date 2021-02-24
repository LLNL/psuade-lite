file(COPY psuadeARSM1.in ValidationData DESTINATION ${PROJECT_BINARY_DIR}/Examples/ARSM )

add_executable (ARSM-sim simulator.c)
target_link_libraries(ARSM-sim m)
SET_TARGET_PROPERTIES(ARSM-sim
  PROPERTIES OUTPUT_NAME simulator)

install(TARGETS ARSM-sim DESTINATION Examples/ARSM)

add_test (ARSM1 ${PROJECT_SOURCE_DIR}/scripts/psuade_1test.sh
  ${PROJECT_BINARY_DIR}/bin/psuade
  ${PROJECT_SOURCE_DIR}/Examples/ARSM/psScript1.in
  ${PROJECT_SOURCE_DIR}/Examples/ARSM/psScript1.out 
  psuadeMetisInfo)

install(FILES ValidationData psuadeARSM1.in DESTINATION Examples/ARSM)
