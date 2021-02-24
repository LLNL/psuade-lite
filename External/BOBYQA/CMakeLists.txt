file(GLOB bobyqa_SRC
  "*.f"
)

#BOBYQA has a call to a function inside PSUADE, so it's undefined at library link time
#On some systems that means we can only build it statically
if(MINGW OR APPLE)
add_library(bobyqa STATIC ${bobyqa_SRC}) 
else(MINGW OR APPLE)
add_library(bobyqa ${LIBRARY_TYPE} ${bobyqa_SRC}) 
endif(MINGW OR APPLE)

install (TARGETS bobyqa
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

