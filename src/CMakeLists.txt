#
# Copyright (c) 2015 CNRS
# Authors: Mathieu Geisert
#
#

SET(LIBRARY_NAME ${PROJECT_NAME})

SET( SOURCES
  #quadrotor_EngineSpeed
#  quadrotorMinTime_EngineSpeed
#  rocket_c
  #quadcopter2
#  quadrotorMinTime
#  quadrotorWithSimpleArm
#  rocket
  quadWithArmPinocchio
  QuadFreeflyerPinocchio
  QuadPickUp
  MPCQuad
#  FreeflyerPinocchio_simpleControl
#  quadWithArmPinocchio_simpleControl
  )

FOREACH( source ${SOURCES} )
  ADD_EXECUTABLE( ${source} ${source}.cpp)
  PKG_CONFIG_USE_DEPENDENCY(${source} acado)
  PKG_CONFIG_USE_DEPENDENCY(${source} pinocchio)
  TARGET_LINK_LIBRARIES(${source} ${Boost_LIBRARIES} /local/mgeisert/quad/libquadspline/lib/libquadspline.so)
  INSTALL(TARGETS ${source} DESTINATION bin)
ENDFOREACH(source)



