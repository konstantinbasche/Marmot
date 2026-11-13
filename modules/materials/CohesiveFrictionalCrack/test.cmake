# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for the CohesiveFrictionalCrack traction-separation law
add_marmot_test("TestCohesiveFrictionalCrack" "${CURR_TEST_SOURCE_DIR}/testCohesiveFrictionalCrack.cpp")
