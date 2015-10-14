set -e -x

# asan makes it painful to get coredumps
MAYBE_ASAN=""
#MAYBE_ASAN=",address"

clang++ \
    --std=c++0x -W -Wall -Wno-sign-compare \
    -O2 -pipe -mmmx -msse -msse2 -msse3 \
    -ggdb \
    -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC \
    -fsanitize=integer,undefined"$MAYBE_ASAN" \
    -fno-sanitize-recover \
    main.cc -o main

time java -jar tester/tester.jar \
    -exec "./driver.sh" -seed 1 #-vis
