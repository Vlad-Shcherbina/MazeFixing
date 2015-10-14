ulimit -c unlimited

# ASAN_OPTIONS="disable_core=0:unmap_shadow_on_exit=1:abort_on_error=1" \
ASAN_SYMBOLIZER_PATH=$(which llvm-symbolizer-3.4) \
MSAN_SYMBOLIZER_PATH=$(which llvm-symbolizer-3.4) \
./main

# to extract coredump:
#   gdb main core -batch -ex bt
