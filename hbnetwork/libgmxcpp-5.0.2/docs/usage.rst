
Usage
=====

.. highlight:: cpp

The basic idea of the library is two-fold and contains two main aspects: 1)
Reading in Gromacs files into memory using constructors and using getters to
access their information in an analysis program, and 2) a set of basic analysis
functions (see next section). Currently libgmxcpp can read in .xtc, .ndx, and
.tpr files (tpr files are limited currently to mass and charge). Below is an
example workflow which contains both of these aspects. The next two sections
contain the API details for the classes and functions.

Workflow
--------

This is a suggested workflow for using this library in constructing one's
analysis program. As an example this tutorial will walk through creating a
program that calculates the center of mass of a group of atoms from a Gromacs
simulation.

Let's say you have simulated several methanes in water. In the case of
calculating the center of mass of the methanes we'll need the .xtc file (having
the coordinates), the .ndx file (grouping the atoms), and the .tpr file (having
the masses).

The first thing to do is to construct an object associated for each file type.
First we'll read in the index file, since we'll be using it to locate the
methanes in the trajectory:::

    Index ndx("index.ndx");

Then we'll read in both the .xtc and .tpr files and associate the Index object
with it. This is optional, but we want to do it in this case since we can easily
find the methanes by our index groups:::

    Trajectory trj("traj.xtc",ndx);
    Topology top("topol.tpr",ndx);

Now all information from the simulation is available to us using object getters
from ``trj`` and ``top``. Since ``ndx`` is now associated with both of these
object we don't have to worry about calling anything from it directly. The first
thing you should do is either read in the entire trajectory, or read in some
frames. To read in the entire xtc file do::

    trj.read();

To read in only one frame do::

    trj.read_next();

To read in the next 10 frames do::

    trj.read_next(10);

``read_next`` is useful in a loop and returns the actual number of frames read
in, so you know when you are at the end of the file. It does not close the xtc
file like ``read`` does. To do so simply call::

    trj.close();

In most cases ``read()`` should be enough unless you are dealing with a large
system and run out of memory.

Now that we've called our constructors, we can get any information we want from
these objects such as atomic coordinates and masses, which is what we need for
getting the center of mass. There is a provided analysis function in the library
which gets the center of mass for a group of atoms, removing the periodic
boundary condition. For this function we need the atomic coordinates of the
atoms in the group we're interested in, the masses of those atoms, and the
simulation box for the particular frame we're interested in. Here I know that my
simulation is using a cubic box so I am using the ``cubicbox`` class instead of
the ``triclinicbox`` class. Here's how we can get that info for the methanes
from the first frame, where we have an index group with the methanes labeled as
``CH4``:::

    vector <coordinates> atom;
    vector <double> mass;
    cubicbox box;

    atom = trj.GetXYZ(0,"CH4");
    box = trj.GetCubicBox(0);
    mass = top.GetMass("CH4");

These getters are described in this documentation on the ``Trajectory`` and
``Topology`` class pages. Now to get the center of mass we just call our
analysis function:::

    coordinates com;

    com = center_of_mass(atom,mass,box);

This only works for frame 0 (the first frame), so to do this for each frame we
would put this into a loop:::

    coordinates com;
    vector <coordinates> atom;
    vector <double> mass;
    cubicbox box;

    Index ndx("index.ndx");
    Trajectory trj("traj.xtc",ndx);
    trj.read();
    Topology top("topol.tpr",ndx);

    for (int i = 0; i < trj.GetNFrames(); i++)
    {
        atom = trj.GetXYZ(i,"CH4");
        box = trj.GetCubicBox(i);
        mass = top.GetMass("CH4");
        com = center_of_mass(atom,mass,box);
    }

At this point outputting the data or averaging it, further analysis is up to
you. Note that we would have to include the appropriate header files to be able
to do this. Additionally the ``for`` loop can possibly be parallelized depending
on the analysis. A full program might be:::

    #include <vector>
    #include "gmxcpp/Index.h"
    #include "gmxcpp/Topology.h"
    #include "gmxcpp/Trajectory.h"
    #include "gmxcpp/Utils.h"
    using namespace std;

    int main()
    {

        coordinates com;
        vector <coordinates> atom;
        vector <double> mass;
        triclinicbox box;

        Index ndx("index.ndx");
        Trajectory trj("traj.xtc",ndx);
        trj.read();
        Topology top("topol.tpr",ndx);

        for (int i = 0; i < trj.GetNFrames(); i++)
        {
            atom = trj.GetXYZ(i,"CH4");
            box = trj.GetBox(i);
            mass = top.GetMass("CH4");
            com = center_of_mass(atom,mass,box);
        }

        return 0;
    }

Compiling a Program
-------------------

Say you have written the above program and saved it to ``com.cpp``. To
compile you need to link your program to libgmxcpp. Additionally if the headers
for your Gromacs installation are in a non-standard installation, which they
most probably are, you need to add that path to the ``CPLUS_INCLUDE_PATH``
environmental variable.

.. highlight:: bash
For example::

    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/gromacs/include
    g++ com.cpp -lgmxcpp

The first line needs to be changed depending on your Gromacs installation and
can be included in your bash profile so you don't have to add it every time you
compile a new program.

Other Examples
--------------

There is an example program in the ``example`` directory. Use ``make`` to compile it
and test it out on an .xtc and .ndx file from a recent simulation.

Additionally `there is an example program which calculates the radial
distribution function using this library <https://github.com/wesbarnett/rdf>`_.

An example of using ``read_next()`` in a loop along with using OpenMP for
parallelization is found `here <https://github.com/wesbarnett/tpi/blob/master/src/main.cpp>`_.
