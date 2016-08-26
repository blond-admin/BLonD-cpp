#Code Profiling
##On Windows Visual Studio
1. Set project mode to Debug
2. Debug -> Show Diagnostic Tools
3. Build Project
4. Debug -> Start Project Without Debugging
5. Performance Explorer -> Attach, Select application under investigation
##On Linux
We use [`valgrind`](http://valgrind.org/info/tools.html#callgrind) to profile and
[`kcachegrind`](http://kcachegrind.sourceforge.net/html/Home.html) to visualise. Note that profiles captured on Linux 
can be visualised on Windows using [`QCacheGrind`](https://sourceforge.net/projects/qcachegrindwin/).
1. Build Debug configuration
2. Run all tests using `ctest` or individual test application under `valgrind`:
```bash
valgrind --tool=callgrind ctest -VV
```
