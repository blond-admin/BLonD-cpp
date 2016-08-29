#How to run benchmarks
##Build project
Make sure that you have `-DWITH_BENCHMARK=True` CMake argument while creating a project. Release build is reconended but 
not obligatory. Different benchmarks may rely on specific compilers or frameworks, so please make sure you have 
inspected documentation for that specific benchmark you want to run. 
##Start benchmark
You can list benchmark arguments calling application with `--help` argument like:
```
./benchTC1_Acceleration --help
```
This would list
```
benchmark [--benchmark_list_tests={true|false}]
          [--benchmark_filter=<regex>]
          [--benchmark_min_time=<min_time>]
          [--benchmark_repetitions=<num_repetitions>]
          [--benchmark_report_aggregates_only={true|false}
          [--benchmark_format=<console|json|csv>]
          [--benchmark_out=<filename>]
          [--benchmark_out_format=<json|console|csv>]
          [--color_print={true|false}]
          [--v=<verbosity>]
```
Now you can run your benchmark calling your application with `--benchmark_out=<filename>` to save results:
```
./benchTC1_Acceleration --benchmark_out=TC1_my.json --benchmark_repetitions=3
```
You can also run with `--benchmark_repetitions=3` argument to get Standard Deviation and Mean.

##Display results
Navigate to `BLonD-minimal-cpp/benchmarks/reports`. Put your JSON file into that folder.
Generated `.json` file shall be registered in `benchmarks_list.js` arter `benchmarks = [` showing your result file name
and its relative URL before other files:
```javascript
benchmarks = [
	{name: "TC1_Acceleration", ref: "./TC1_my.json", axes: ["particles", "turns", "slices"]},
	{name: "other_data", ref: "./other.json", axes: ["other"]}
];
```
Now open `benchmark_browser.html` in [FireFox](https://www.mozilla.org/en-US/firefox/) browser if you are doing it 
locally, or share folder with it on your web
site. Note in case of mean results presence renderer filter's all other results out.

##Display multiple results
Yoe can display results created by running same benchmark and saving output into different files. This is helpful for 
testing on different hardware/platforms/compiler arguments on same machine or for different clusters. Just add multiple 
references to `.json` files into `benchmarks_list.js` `ref` in form of [JSON array](http://www.json.org/) alike:
```javascript
benchmarks = [
	{name: "TC1_Acceleration", ref: ["./TC1_my.json" ,"./TC1_clusterA.json" ,"./TC1_clusterB.json" ], axes: ["particles", "turns", "slices"]},
	{name: "other_data", ref: "./other.json", axes: ["other"]}
];
```

##Start benchmark on SLURM cluster
Create a simple scritpt called `blond.sh` like this featuring your favourite benchmark:
```bash
#!/bin/sh
#SBATCH -p phi
#SBATCH -c 48
export OMP_NUM_THREADS=48
srun ./benchTC1_Acceleration --help
srun ./benchTC1_Acceleration --benchmark_out=TC1_hybri_AfterCpu.json --benchmark_repetitions=3
srun ./benchTC1_Acceleration --benchmark_out=TC1_hybri_AfterCpu_extra.json --benchmark_repetitions=3
```
This script would create `TC1_hybri.json` file as a result.
Start it:
```bash
sbatch ./blond.sh
```
Monitor its progress calling:
```bash
squeue
```
That would show something alike:
```
JOBID   PARTITION   NAME        USER    ST  TIME        NODES   NODELIST(REASON)
20950   gpuK80      blond.sh    you     R   19:19:42    1       blade10
```
When `squeue` call will stop showing your job in it you can use result file.

To force application to stop use:
```bash
scancel $YOUR_JOBID
```
