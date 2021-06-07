Install with:


```
git clone https://github.com/Astrokiwi/gtools.git gtools

cd gtools/
pip install -e . -r requirements.txt
```

Installing with `-e` means 'editable', so you can modify gtools and run without having to reinstall. You can skip `-r requirements.txt` if you think you have all the standard requirements installed (numpy, matplotlib, pynbody, scipy, etc)

The first time you run something, it will ask which directories you keep your GIZMO directories, where you want movies to be output to, and where Cloudy tables are stored. You can rechoose the directories with

```
python -m gtools.config
```

Note that the Cloudy tables are only required if you are calculating dust temperature, emissivity, heating/cooling rates, or any other value that depends on radiative transfer and isn't stored in snapshots.

To load GIZMO data and convert it into a pynbody snapshot for investigation, do:

```
import pynbody
from gtools import gizmo_tools

header,snap = gizmo_tools.load_gizmo_nbody('<GIZMO_RUN_DIRECTORY>','<DATA_OUTPUT_SUBDIRECTORY>','<SNAPSHOT_STR>'


```

where `<GIZMO_RUN_DIRECTORY>` is the directory that contains the GIZMO executable, `<DATA_OUTPUT_SUBDIRECTORY>` is the subdirectory that contains the `snapshot_???.hdf5` outputs, `<SNAPSHOT_STR>` is the number of the snapshot as a 3-digit string with leading zeros (e.g. `025`,`005`). The `header` is a `dict` containing useful information, and `snap` is the pynbody snapshot. See `https://pynbody.github.io/pynbody/` for information on pynbody. You can do almost everything just by using pynbody, once things are in pynbody format, but I have some custom routines myself.

You can animate GIZMO SPH output with:


```
python -m gtools.anim <GIZMO_RUN_DIRECTORY> <DATA_OUTPUT_SUBDIRECTORY> --nprocs [PROCS] --plot [quantity, e.g. nH, temp, dens] --rad [SIZE in pc]
```

You can create an animation of a phaseplot (e.g. nH vs T) using
```
python -m gtools.anim_phaseplots <GIZMO_RUN_DIRECTORY> <DATA_OUTPUT_SUBDIRECTORY>
```
You must edit `gtools/anim_phaseplots.py` if you want to make a different phaseplot, e.g. plotting different quantities.

You can rotate around a snapshot with:
```
python -m gtools.sph_rotanim <GIZMO_RUN_DIRECTORY> <DATA_OUTPUT_SUBDIRECTORY> <SNAPSHOT_NUMBER>
```
Again, you must edit `gtools/sph_rotanim.py` if you want to modify what quantity is being plotted.


Contact `david.john.williamson@gmail.com` if you have any problems
