# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for CohesiveInterfaceElement
add_marmot_test("TestCohesiveInterfaceElement" "${CURR_TEST_SOURCE_DIR}/TestCohesiveInterfaceElement.cpp")
