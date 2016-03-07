add_test(chemical_reactions ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/chemical_reactions.py)
set_tests_properties(chemical_reactions PROPERTIES ENVIRONMENT "${TEST_ENV}")
