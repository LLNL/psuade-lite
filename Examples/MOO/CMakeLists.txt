
file(COPY configFile psuade.in moo.py DESTINATION ${PROJECT_BINARY_DIR}/Examples/MOO )
#${PROJECT_SOURCE_DIR}/Examples/MOO/

if (USE_BOBYQA)
  add_test (MOO ${PROJECT_SOURCE_DIR}/scripts/psuade_1test.sh ${PROJECT_BINARY_DIR}/bin/psuade ${PROJECT_SOURCE_DIR}/Examples/MOO/psScript.in ${PROJECT_SOURCE_DIR}/Examples/MOO/psScript.out)

  if(UNIX)
    install(FILES psuade.in DESTINATION Examples/MOO)
    install(PROGRAMS moo.py DESTINATION Examples/MOO)  #PROGRAMS should give moo.py execute permissions
    install(FILES configFile DESTINATION Examples/MOO)
  else(UNIX)
    install(FILES psuadeForWin.in DESTINATION Examples/MOO)
    install(PROGRAMS moo.py DESTINATION Examples/MOO)  #PROGRAMS should give moo.py execute permissions
    install(FILES configFile DESTINATION Examples/MOO)

  endif(UNIX)
endif (USE_BOBYQA)