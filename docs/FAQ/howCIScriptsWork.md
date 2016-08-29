#How CI scripts work
CI scripts are composed in form of YAML markup. They define build configurations that are tested each time new commit 
into repository happens.

##Select a platform
###Linux
Base image for travis ci is Ubuntu 14.04 with GCC 4.8 (one that can be found on most clusters and personal computers)
####Install Packets
We set `sources` and list `packages` in  `addons->apt`. Try to keep packages list as small as possible
###Windows
On AppVeyor we use Windows Server, MSVC 2015

###Build external libraries
We use cmake arguments to build all required external libraries (do not install them into system)
###Build blond
Use `cmake --build .` to keep it crossplatform and make system independent
###Run Unit-Tests
In CI scripts we only run unit tests - no benchmarks