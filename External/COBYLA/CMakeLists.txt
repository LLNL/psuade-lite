add_library(cobyla2c ${LIBRARY_TYPE} cobyla.c)

install (TARGETS cobyla2c
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)