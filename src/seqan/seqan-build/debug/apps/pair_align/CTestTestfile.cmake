# CMake generated Testfile for 
# Source directory: /home/svilsen/seqan/seqan-src/apps/pair_align
# Build directory: /home/svilsen/seqan/seqan-build/debug/apps/pair_align
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(app_test_pair_align "/usr/bin/python" "/home/svilsen/seqan/seqan-src/apps/pair_align/tests/run_tests.py" "/home/svilsen/seqan/seqan-src" "/home/svilsen/seqan/seqan-build/debug")
SUBDIRS(lib)
