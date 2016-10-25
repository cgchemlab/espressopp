add_test(change_particle_type_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/change_particle_type_test.py)
set_tests_properties(change_particle_type_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
