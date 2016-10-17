add_test(particle_group_by_type_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/particle_group_by_type_test.py)
set_tests_properties(particle_group_by_type_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
