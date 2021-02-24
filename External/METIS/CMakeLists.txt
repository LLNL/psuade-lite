file(GLOB metis_SRC
  "Lib/*.c"
)
add_library(metis ${LIBRARY_TYPE} ${metis_SRC}) 

install (TARGETS metis
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
