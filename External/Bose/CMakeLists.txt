#file(GLOB bose_SRC
#  "*.c"
#)
add_library(bose ${LIBRARY_TYPE} construct.c akconst.c primes.c gfields.c galois.c memory.c boselink.c  oa.c)

install (TARGETS bose
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

