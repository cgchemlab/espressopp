add_test(topology_manager ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/topology_manager.py)
set_tests_properties(topology_manager PROPERTIES ENVIRONMENT "${TEST_ENV}")
