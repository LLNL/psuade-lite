file(GLOB lbfgs_SRC
  "src/*.c"
)

#LBFGS has a call to a function inside PSUADE, so it's undefined at library link time
#On some systems that means we can only build it statically
if(MINGW OR APPLE)
add_library(lbfgs STATIC ${lbfgs_SRC}) 
else(MINGW OR APPLE)
add_library(lbfgs ${LIBRARY_TYPE} ${lbfgs_SRC}) 
endif(MINGW OR APPLE)

install (TARGETS lbfgs
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

