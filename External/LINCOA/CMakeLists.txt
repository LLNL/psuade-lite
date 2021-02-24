file(GLOB lincoa_SRC
  "*.f"
)

#LINCOA has a call to a function inside PSUADE, so it's undefined at library link time
#On some systems that means we can only build it statically
if(MINGW OR APPLE)
add_library(lincoa STATIC ${lincoa_SRC}) 
else(MINGW OR APPLE)
add_library(lincoa ${LIBRARY_TYPE} ${lincoa_SRC}) 
endif(MINGW OR APPLE)

install (TARGETS lincoa
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

