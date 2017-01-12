# CMake generated Testfile for 
# Source directory: /home/svilsen/seqan/seqan-src/dox
# Build directory: /home/svilsen/seqan/seqan-build/debug/dox
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(build_dox "/usr/bin/python" "/home/svilsen/seqan/seqan-src/dox/../util/bin/dox.py" "-b" "/home/svilsen/seqan/seqan-src/dox/.." "-i" "/home/svilsen/seqan/seqan-src/dox/../include/seqan" "-i" "/home/svilsen/seqan/seqan-src/dox/pages" "--image-dir" "/home/svilsen/seqan/seqan-src/dox/images")
