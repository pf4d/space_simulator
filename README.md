space_simulator
===============

Space simulator!

install dependencies :
----------------------

This game requires Pylab be installed:

```bash
sudo apt-get install numpy scipy
```

On Ubuntu 13.10, install as instructed at https://code.google.com/p/pyftgl/wiki/HowToBuild :

install FTGL:

```bash
sudo apt-get install libftgl-dev
```

install boost-python:

```bash
sudo apt-get install libboost-python-dev
```

install PyFTGL:

```bash
sudo easy_install PyFTGL
```

Run the code :
--------------

In the ```src``` directory, simply type

```bash
python bounce.py
```

To move the camera, press the left mouse button and drag up and down.  There are two 'axes' in the lower-right hand side of the screen --  the right is the translational velocity (green) and acceleration (blue), while the left is the rotational velocity (yellow) and acceleration (red).  in the upper-right-hand side of the screen you see the number of particles, which in the current configuration are planets.

There is a periodic boundary on all edges of a box in the center of the screen.

The balls are normally colored orange, but when they are accelerating above a threshold, they are colored gray.
