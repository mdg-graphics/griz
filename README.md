# Griz
A postprocessing visualization tool for Mili databases

## Build

### Configuration

From within the directory `griz/Src` exectute the command:

```
configure
```

If you need to use a version of Mili other than the default installation on LC run:

```
configure --with-mili=<path-to-mili>
```

These will generate a new build directory starting with `GRIZ4`, for example: `GRIZ4-toss_3_x86_64_ib-RZGENIE`

### Compilation

Once the build directory has been generated using `configure`, the code is compiled within this new directory using the `gmake` utility.

There are currently 4 standard builds that are commonly made:

|Build Type     |Description|
|---------------|--------------|
|debug          |Debug version of Griz with gui            |
|opt            |Optimized version of Griz with gui        |
|batchdebug     |Debug version of batch griz (no gui)      |
|batchopt       |Optimized version of batch griz (no gui)  |

These are compiled by running the command:

```
gmake <build>
```

For example:

```
gmake debug; gmake opt; gmake batchdebug; gmake batchopt
```

## Testing  (**LLNL Specific**)

A suite of tests exists within the directory `griz/Src/test`. To run the tests perform the following steps:

### 1. Checkout Scripts Directory

`Test.py` lives in `Dyna3dx/CheckOut/scripts/Tests.py`. griz_env lives there also.

To checkout locally run: `cvs co -d scripts/CheckOut/scripts` in the directory `griz/Src/test`

### 2. Running Tests

Tests are run using the script Test.py. The argument `-e griz_env` must 
be added to run the griz tests

To find all the options just run: `scripts/Test.py -h`

The tests results were baselined using the `batchopt` build of Griz.
You should be using this version when running the tests otherwise some
of the tests will fail


To run the Griz Test Suite:

All Tests:
```
scripts/Test.py -e griz_env -c <path_to_code> -p1 -s all
```

Single Suite:
```
scripts/Test.py -e griz_env -c <path_to_code> -p1 -s <suite_name>
scripts/Test.py -e griz_env -c <path_to_code> -p1 -s derived
```

Single Test:
```
scripts/Test.py -e griz_env -c <path_to_code> -p1 -t <testname>
scripts/Test.py -e griz_env -c <path_to_code> -p1 -t bar71_image_press_rmin_th
```

## Deployment on OCF (**LLNL Specific**)

To do deploy griz on the OCF perform the following steps

### 1. Update Version and Timestamp

In the director `griz/Src`, edit `configure.ac` and change the following

- If we are changing the version, change appropriately
```
    m4_define(GRIZ_VERSION_M4, V21_01)
    m4_define(GRIZ_MAJOR_M4, 20)
    m4_define(GRIZ_MINOR_M4, 1)
    m4_define(GRIZ_BUG_M4, 1)
```

- Alway change the date (you can change the time if you want, but usually we do not)
```
    m4_define(GRIZ_DATE_M4, mm/dd/yyyy)
    m4_define(GRIZ_TIME_M4, 07:00:00)
```
    
Save file and run

```
autoconf -f
```

### 2. Configure and Build

Build the `opt` and `batchopt` versions of griz

```
gmake clean; gmake opt; gmake batchopt;
```

### 3. Give executables to mdgadmin

```
give mdgadmin bin_opt/griz.linux-gnu_opt bin_batch_opt/griz.linux-gnu_opt_batch
```

### 4. Login as mdgadmin

```
xsu mdgadmin
```

### If creating a NEW version

If this is a new version then we need to create a new directory. Follow the steps below:

5. `cd /collab/usr/apps/mdg/archive/grizdir`
6. `cp -p -r [Last Version Directory name (i.e V16_01)] [New version name]`
7. `rm grizalpha`
8. `ln -sf [New version name] grizalpha`
9. `cd grizalpha/bin`
10. `cp griz.linux-gnu_opt griz.linux-gnu_opt-MM.DD.YYYY`
11. `cp griz.linux-gnu_opt_batch griz.linux-gnu_opt_batch-MM.DD.YYYY`
12. `take -f <your username>`
13. `chmod g+rx griz.linux-gnu_opt`
14. `chmod g+rx griz.linux-gnu_opt_batch`
15. **Most Important**: logout as mdgadmin and run a test using the new installed alpha version: `griz -alpha -i <filename>`


### NOT a new version

This will install griz without killing running jobs

5. `cd /collab/usr/apps/mdg/archive/grizdir/[grizalpha|grizbeta|grizversion]/bin`
6. `mv griz.linux-gnu_opt_batch griz.linux-gnu_opt_batch-MM.DD.YYYY`
7. `mv griz.linux-gnu_opt griz.linux-gnu_opt-MM.DD.YYYY` (MM.DD.YYYY should be date of version being moved)
8. `take <your username>`
9. `chmod g+rx griz.linux-gnu_opt`
10. `chmod g+rx griz.linux-gnu_opt_batch`
11. **Most Important**: logout as mdgadmin and run a test using the new installed alpha version: `griz -alpha -i <filename>`
