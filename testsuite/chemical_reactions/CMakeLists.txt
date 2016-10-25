add_test(chemical_reactions ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/chemical_reactions.py)
add_test(atrp_activator ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/atrp_activator.py)
set_tests_properties(chemical_reactions PROPERTIES ENVIRONMENT "${TEST_ENV}")
set_tests_properties(atrp_activator PROPERTIES ENVIRONMENT "${TEST_ENV}")
