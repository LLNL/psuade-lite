file(GLOB nomad_SRC
  "src/*.cpp"
)

#NOMAD has a call to a function inside PSUADE, so it's undefined at library link time
#On some systems that means we can only build it statically
if(MINGW OR APPLE)
add_library(nomad STATIC ${nomad_SRC}) 
else(MINGW OR APPLE)
add_library(nomad ${LIBRARY_TYPE} ${nomad_SRC}) 
endif(MINGW OR APPLE)

install (TARGETS nomad
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

