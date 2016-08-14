add_test(periodic_bc ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_periodic_bc.py)
set_tests_properties(periodic_bc PROPERTIES ENVIRONMENT "${TEST_ENV}")
