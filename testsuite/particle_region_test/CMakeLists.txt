add_test(particle_region_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/particle_region_test.py)
set_tests_properties(particle_region_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
