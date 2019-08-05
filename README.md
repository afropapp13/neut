# NEUT

## Building and linking against

### Setup prequisitings
The configure step expects to be able to find `root-config` and expects `CERN` and `CERN_LEVEL` to be defined

### Configure and build
```
cd /path/to/neut
mkdir build; cd build;
../src/configure --prefix=$(readlink -f Linux)
make -j 8
make install
```

### Sourcing the NEUT environment
```
source /path/to/neut/build/Linux/setup.sh
```

### Using neut-config
If the environment is set up then `neut-config` should be available. To build and link an app that needs to be able to read neut vectors you could use
```
-I$(neut-config --incdir) -L$(neut-config --libdir) $(neut-config --iolibflags)
```
which might produce compiler/linker flags like:
```
-I/path/to/neut/build/Linux/include -L/path/to/neut/build/Linux/lib -lneutclassUtils -lneutclass
```

### If you need to link against the old style NEUT libraries
You can configure with the `--enable-compatnames` argument which will build the NEUT libraries in the old format. i.e. each ROOT I/O class gets its own library, the main NEUT static libraries are built per subdirectory, and the reweight library is called `libNReWeight.so`.
N.B. in this mode, highly parallel builds are not supported, if you build errors related to `XXXDict.cc` while building with `make -j X`, try running `make` with no `-j` argument and see if they are resolved.

## When you found issues or bugs in the existing code
Please do not commit to the master.
There is a branch for the next version and bug fix release.
Please use those branches to commit fixes.
(If there are not, please let us know.)

## When you try to implement a new model
Please make a branch to start including a new model implementation.
Make a new directory under src for your new model.
Try to minimize the changes in the existing part.

## Contact
Hayato-san: hayato@icrr.u-tokyo.ac.jp
Luke Pickering: picker24@msu.edu
Clarence Wret: c.wret@rochester.edu
