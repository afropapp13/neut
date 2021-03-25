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

*N.B.* If you get autotools whinging about missing LIBTOOL or AM_CONDITIONAL macros, try:

```
cd /path/to/neut
cd src
rm aclocal.m4
autoreconf -if
cd ../build;
../src/configure --prefix=$(readlink -f Linux)
make -j 8
make install
```

(Generally most autotools problems are solved by running  `autoreconf`)

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

## When you found issues or bugs in the existing code
Please do not commit to the master.<br>
There is a branch for the next version and bug fix release.<br>
Please use those branches to commit fixes.<br>
(If there are not, please let us know.)<br>

## When you try to implement a new model
Please make a branch to start including a new model implementation.<br>
Make a new directory under src for your new model.<br>
Try to minimize the changes in the existing part.<br>

## Contact
Hayato-san: hayato@icrr.u-tokyo.ac.jp<br>
Luke Pickering: picker24@msu.edu<br>
Clarence Wret: c.wret@rochester.edu<br>
