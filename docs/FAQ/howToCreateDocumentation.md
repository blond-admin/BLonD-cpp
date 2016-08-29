#How to create documentation
##Source code and .md documentation
### Requirements
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) 
* [Graphviz](http://www.graphviz.org/Download..php)
### Build
To generate html documentation with search and graphical class hierarchy's  and run:
```bash
cd BLonD-minimal-cpp/
doxygen Doxyfile
```

A `BLonD-minimal-cpp/documentation` folder will appear containing `.html` documentation.

##Unit Tests documentation (Linux only)
### Requirements
* [genhtml](http://linux.die.net/man/1/genhtml) 
* [lcov](http://ltp.sourceforge.net/coverage/lcov.php)

### Build
To generate html documentation showing which parts of your code aren’t covered by our test suite:
1) Compile project in debug mode, build and run unit tests.
2) Run:
```bash
cd BLonD-minimal-cpp/build
lcov --capture --directory .. --output-file coverage.info
genhtml coverage.info --output-directory cover
```

A `BLonD-minimal-cpp/cover` folder will appear containing `.html` documentation showing unit tests coverage.
##Benchmark documentation
### Requirements
* [FireFox](https://www.mozilla.org/en-US/firefox/new/)
### Build
1) See [how to run benchmarks](./howToRunBenchmarks.md) document to compile and run benchmarks

Static html files with Benchmark documentation are avaliable at `BLonD-minimal-cpp/benchmarks/reports` folder.
