# CreateSS

These files are for stand-alone GITM to create a steady-state solution that
other simulations can build off of.

Files in this directory:

- `UAM.in`: a GITM stand-alone input file for reference.
- `PARAM.in.gitmss`: an SWMF input file for GITM as the sole module.

GITM quirks and things to know about the UAM.in:
1a) #GRID is not smart. enough to handle a random number of processors. If lat and lon are set to 20 and 20 then you will need to run on 20x20=400 processors.

1b) Getting a desired resolution: If lat and lon are set to 20 and 20 you will have a resolution of 1 deg by 2 deg. If you need a high resolution but have fewer processors then you can adjust nlons and nlats in ModGrid.f90. So nlon (default=9) x lon (userdefined=20) = 180 longitudinal grid cells -> 2 deg resolution. If nlon (userchanged=18) x lon (userdefined=10) = 180 longitudinal grid cells -> 2 deg resolution && nlat (userchanged=18) x lat (userdefined=10) = 180 latitudinal grid cells -> 1 deg resolution == only requires 10 x 10 =100 processors. Clear as mud?

2) #SAVEPLOTS has a handful of plot output options. 3DALL and 2DMEL should have all that you would need but I've included 2DGEL as well just in case.

3) #ELECTRODYNAMICS It is "normal" to update the electrodynamics every 60s in stand alone runs. I have mixed feelings about this and usually update it more often (the cost is that it makes the run take longer; electrodynamics is a beast).

4) There are #s for constant drivers (IMF, F107, HPI) or time varying drivers. They are grouped together in the UAM.in file so you can (find them) switch between them easily. The time varying options use the same input file format that you generate for the SWMF.
