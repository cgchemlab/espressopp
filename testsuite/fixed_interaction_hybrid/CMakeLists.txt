add_test(fixed_interaction_hybrid ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fixed_list_hybrid_test.py)
set_tests_properties(fixed_interaction_hybrid PROPERTIES ENVIRONMENT "${TEST_ENV}")
