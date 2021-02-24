add_executable (pdf-sim simulator.c)
SET_TARGET_PROPERTIES(pdf-sim
  PROPERTIES OUTPUT_NAME simulator)

install(TARGETS pdf-sim DESTINATION Examples/PDFTest)

add_test (pdf1 ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psuade1.in  ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psScript1.out)
add_test (pdf2 ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psuade2.in  ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psScript2.out)
add_test (pdf3 ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psuade3.in  ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psScript3.out)
add_test (pdf4 ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psuade4.in  ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psScript4.out)
add_test (pdf5 ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psuade5.in  ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psScript5.out)
add_test (pdf6 ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psuade6.in  ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psScript6.out)
add_test (pdf7 ${PROJECT_SOURCE_DIR}/scripts/psuade_0test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psuade7.in  ${PROJECT_SOURCE_DIR}/Examples/PDFTest/psScript7.out)

install(FILES psuade1.in  psuade2.in psuade3.in psuade4.in psuade5.in psuade6.in psuade7.in DESTINATION Examples/PDFTest)

