#How to run benchmarks
##Build project
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
./benchTC1_Acceleration --benchmark_out=TC1_hybri.json
```

##Display results
Navigate to `BLonD-minimal-cpp/benchmarks/reports`. Put your JSON file into that folder.
Generated `.json` file shall be registered in `benchmarks_list.js` arter `benchmarks = [` showing your result file name
and its relative URL before other files:
```javascript
benchmarks = [
	{name: "My_TC1_Accelleration", ref: "./TC1_hybri.json"},
	{name: "other", ref: "./other.json"},
];
```
Now open `benchmark_browser.html` in [FireFox](https://www.mozilla.org/en-US/firefox/) browser if you are doing it locally, or share folder with it on your web
site.

##Start benchmark on SLURM cluster
Create a simple scritpt called `blond.sh` like this featuring your favourite benchmark:
```bash
#!/bin/sh
#SBATCH -p phi
srun ./benchTC1_Acceleration --help
srun ./benchTC1_Acceleration --benchmark_out=TC1_hybri.json
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