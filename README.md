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
python spaceSimulator.py
```

Controls :
----------

* ```Page Down``` : ascend
* ```Page Up``` : descend
* ```w``` : fire rear thruster (move forward)
* ```s``` : fire forward thruster (move backward)
* ```a``` : fire left side thruster (move right)
* ```d``` : fire right side thruster (move left)
* ```q``` : fire right nose thruster (yaw left)
* ```e``` : fire left nose thruster (yaw right)
* ```arrow up``` : fire bottom nose thruster (pitch up)
* ```arrow down``` : fire top nose thruster (pitch down)
* ```arrow left``` : fire bottom right wing thruster (roll left)
* ```arrow right``` : fire bottom left wing thruster (roll right)

To move the camera, press the left mouse button and drag up and down, and zoom in and out with the mouse wheel.  

Cockpit readouts :
------------------

There are two miniture ship images in the lower-right hand side of the screen --  the right is the translational velocity (green) and acceleration (blue), while the left is the rotational velocity (yellow) and acceleration (red).  in the upper-right-hand side of the screen you see the number of particles, which in the current configuration are planets, and in the upper-left you have frames per second.

Colors of planets :
-------------------

The balls are normally colored orange, but when they are accelerating above a threshold, they are colored gray.
