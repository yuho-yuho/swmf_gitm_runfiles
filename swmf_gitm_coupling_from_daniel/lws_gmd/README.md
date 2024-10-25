# LWS GMD Param Files

These files are for the NASA Living With a Star GeoMagnetic Disturbance project, led by Dr. Yue Deng. The goal is to investigate the role of the thermosphere in creating and modulating the formation of ground magnetic signatures.

## Idealized Benchmark Simulations
The basic idealized input files run the SWMF/Geospace model without GITM to create a set of drivers for GITM stand-alone.

The original specification of the runs is as follows:
1.	Resolution of GITM runs: 1 deg lat x 1 deg lon. Output every 5 second to provide time series of dB to Steve.
2.	Equinox condition and seasonal effect in future maybe. (F10.7=150 s.f.u.)
3.	5-hour simulation, 1-hour quiet, 3-hour active, 1-hour recovery; Step function for forcing changes. For the SWMF runs, we include an extra hour of quiet time to build a reliable steady state solution (6 hours simulation total).

### Solar Driving Conditions

The run numbers correspond to different solar driving conditions, given in the table below. Thermal pressures are selected to get a plasma temerature of 96831.2 $K$, which is a typical background value.

| Run # | $V$ ($km/s$) | $P$ ($nPa$) | $\rho$ ($1/cc$) | $B_Y$ ($nT$) | $B_Z$ ($nT$) |
|-------|-------------------|------------------|------------------|--------------|-------------|
| Run 00 | 400 | 0.005 | 3.74 | 0 | -1 |
| Run 01 | 400 | 0.005 | 3.74 | 0 | -1 $\rightarrow$ -20 |
| Run 02 | 400 | 0.005 $\rightarrow$ 0.025 | 3.74 $\rightarrow$ 18.68| 0 | -1 $\rightarrow$ -20 |
| Run 03 | 400 $\rightarrow$ 800 | 0.005 $\rightarrow$ 0.00625 | 3.74 $\rightarrow$ 4.67  |0 | -1 $\rightarrow$ -20 |
| Run 04 | 400 $\rightarrow$ 800 km/s | 0.005 $\rightarrow$ 0.00625 | 3.74 $\rightarrow$ 4.67 | 0 $\rightarrow$ 10 | -1 $\rightarrow$ -20|

### Code Configuration

Code config and compilation follows the SWPC configuration but with higher resolution in the IE solver:

```
./Config.pl -install= -compiler=gfortran
./Config.pl -v=GM/BATSRUS,IE/Ridley_serial,IM/RCM2
./Config.pl -o=IE:g=181,361
make SWMF PIDL rundir
```

### Running the Code

Copy all PARAMs input files into your SWMF run directory.
Replace PARAM.in with the appropriate PARAM file from the repo:

```
rm PARAM.in
ln -s PARAM.in_init_RCM PARAM.in
```

Change the `#SOLARWIND` command to use the IMF file matching the desired simulation (see above table.)

On local machines, just run the code, either parallel or serially:
```
./SWMF.exe | tee log.txt
```
or:
```
mpiexec -n 8 ./SWMF.exe | tee log.txt
```

On parallel machines, submit a jobscript. Figure it out. I believe in you.