#How to auto format code
## Requirements
 * clang-format shall be [installed](http://llvm.org/releases/download.html)
## Configuration
To auto-reformat changed code on each build using clang-format please add `-DWITH_FORMAT=True` argument 
to configuration command, before `..`.
## Run
Build sol