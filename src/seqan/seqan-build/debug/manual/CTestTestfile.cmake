# CMake generated Testfile for 
# Source directory: /home/svilsen/seqan/seqan-src/manual
# Build directory: /home/svilsen/seqan/seqan-build/debug/manual
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(build_manual "/usr/bin/sphinx-build" "-W" "-n" "-b" "html" "-d" "/home/svilsen/seqan/seqan-build/debug/manual/doctrees" "/home/svilsen/seqan/seqan-src/manual/source" "/home/svilsen/seqan/seqan-build/debug/manual/html")
SET_TESTS_PROPERTIES(build_manual PROPERTIES  DEPENDS "build_dox")
