add_library(mars ${LIBRARY_TYPE} mars36_fort.f)

install (TARGETS mars 
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
