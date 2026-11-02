# Four-cylinder engine simulator

A chain of MATLAB scripts that follows a four-stroke, four-cylinder spark-ignition engine from the geometry of its crank, through the gas in its cylinders, to the gear ratios of the car it drives.

![MATLAB](https://img.shields.io/badge/MATLAB-scripts-0076A8)
![GNU Octave](https://img.shields.io/badge/GNU%20Octave-runs%20unchanged-0790C0)
![Course project](https://img.shields.io/badge/RHU-Automotive%20Engineering%2C%20fall%202024-555)

[![The four cylinders running through one 720 degree cycle beside live charts of cylinder pressure, the pressure-volume loop and crank torque](docs/media/engine-simulation.gif)](docs/media/engine-simulation.mp4)

*One full cycle of the simulated engine, drawn from the scripts' own output. Click it for the 1080p60 video.*

## Watch it run

| [![Still from the engine video: cylinder 1 just after ignition, beside the pressure, P–V and torque charts](docs/media/engine-simulation-poster.jpg)](docs/media/engine-simulation.mp4) | [![Still from the run video: the console output of part2.m beside the closing results card](docs/media/simulation-run-poster.jpg)](docs/media/simulation-run.mp4) |
| --- | --- |
| **[The engine, simulated](docs/media/engine-simulation.mp4)**: 24 s, 1080p60, three cycles 133× slower than 2000 rpm. The pistons, rods and cranks move to the kinematics of `part1.m`, drawn to scale. The gas column is the cylinder volume, coloured by the cycle's temperature. The flame grows with the Wiebe burn, and each port lights up while its valve is open (on the timing in `constants.m`; the lift is drawn, not modelled). The charts follow the same crank angle: pressure in all four cylinders, the P–V loop of cylinder 1, and the crank torque. | **[A full run, every output](docs/media/simulation-run.mp4)**: 52 s, 1080p. The three scripts run in GNU Octave. The console replays every line the run printed, from the cycle converging to the gear ratios. Each of the nine figures the scripts open draws itself in as it opens, followed by a chart of the gear ratios `part2.m` prints. The video ends on the results. |

## What it is

I wrote this simulator as the main project of the Automotive Engineering course at RHU, in fall 2024. It models one engine at 2000 rpm over a full cycle of 720° of crank rotation.

An engine is one chain of dependencies. The crank geometry sets the cylinder volume, the volume sets the pressure, the pressure and the moving masses set the torque, and the torque decides the flywheel and the gearbox. The scripts are built the same way: every stage takes the previous one's output as its input, so a change to one parameter runs through every result downstream.

## Results

For the engine in the assignment:

| Quantity | Value | Notes |
| --- | --- | --- |
| Engine power | 30.6 kW | indicated, four cylinders at 2000 rpm |
| Mean torque | ~146 N·m | four cylinders, 180° apart |
| Fuel conversion efficiency | 28% | |
| Work per cylinder per cycle | ~460 J | area of the P–V loop |
| Peak cylinder pressure | ~34 atm | just after top dead centre |
| Flywheel radius | 13.6 cm | steel disc 50 mm thick, 4% speed fluctuation |
| Gear ratios (1st to 8th) | 5.38 / 3.60 / 2.49 / 1.78 / 1.31 / 1.00 / 0.79 / 0.64 | sixth is direct drive |

## How it works

```mermaid
flowchart LR
    A["constants.m<br/>every parameter"] --> B["part1.m<br/>slider-crank kinematics<br/>and cylinder volume"]
    B --> C["part2.m<br/>thermodynamic cycle"]
    C --> D["part2.m<br/>crank torque,<br/>four cylinders summed"]
    D --> E["part2.m<br/>flywheel and<br/>performance"]
    E --> F["part2.m<br/>eight-speed gearbox"]
    F --> G["animation.m<br/>four slider-cranks"]
```

**1. One file of numbers.** `constants.m` holds every input: link lengths and masses, the engine speed, the timing of intake, compression, combustion, expansion and exhaust in crank degrees, the gas properties, and the flywheel and transmission constants. Each later script loads the `.mat` file the one before it saved.

**2. The slider-crank.** `part1.m` solves the slider-crank at every crank angle: piston position, velocity and acceleration, the angle of the connecting rod, and the cylinder volume. It does this twice, once with the cylinder axis through the crankshaft (compression ratio 9.8) and once offset by −18 mm (10.1). The same pass computes the rod's angular velocity and acceleration and the accelerations of each part's centre of mass, which the torque needs later.

![Piston velocity against crank angle with and without an 18 mm crank offset](docs/media/piston-velocity.webp)

| | |
| --- | --- |
| ![Piston distance from top dead centre against crank angle, with and without the offset](docs/media/piston-displacement.webp) | ![Cylinder volume against crank angle, with and without the offset](docs/media/cylinder-volume.webp) |
| Piston travel from top dead centre. | Cylinder volume, from the clearance volume up to bottom dead centre. |

**3. A cycle that has to agree with itself.** `part2.m` runs a single-zone thermodynamic cycle over that volume. Compression and expansion are polytropic. Combustion takes 60° around top dead centre: a constant-volume heat addition from the fuel's heating value gives the end-of-combustion temperature, and a Wiebe function (a = 5, m = 2) sets how much of the charge has burned at each angle in between. Exhaust blowdown and intake are compressible flow through a valve orifice, and the leftover exhaust gas mixes with the fresh charge. That mixture sets the temperature at the start of the next compression, so the script repeats the cycle until the temperature at the end of intake matches the one it started compression with.

| | |
| --- | --- |
| ![Cylinder pressure against crank angle, with intake, compression, combustion, expansion and exhaust shaded](docs/media/pressure-crank-angle.webp) | ![Pressure-volume diagram of one cylinder over a full cycle, coloured by phase](docs/media/pv-diagram.webp) |
| Pressure over the 720° cycle. | The P–V loop of one cylinder; its area is the work per cycle. |

**4. From gas pressure to crank torque.** At every crank angle the script writes the force and moment balances for the crank, the connecting rod and the piston (their masses, the rod's inertia, gravity and the gas force) as an 8×8 linear system and solves it for the crank torque. The other three cylinders are the same curve shifted by 180°, 360° and 540°, and the sum is the torque of the engine.

| | |
| --- | --- |
| ![Crank torque of one cylinder over the cycle, with its mean](docs/media/torque-single-cylinder.webp) | ![Crank torque of four cylinders 180 degrees apart, summed, with its mean](docs/media/torque-four-cylinder.webp) |
| One cylinder: one power stroke in 720°, a mean of 36.6 N·m. | Four cylinders 180° apart: a mean of about 146 N·m. |

**5. Flywheel and performance.** The script integrates the energy of the torque above its mean between the crossings and sizes a flywheel to hold the speed fluctuation to a coefficient of 0.04, as a steel disc 50 mm thick. The same run prints the work per cycle, indicated and brake power, MEP and BMEP, efficiencies, specific fuel consumption, volumetric efficiency and mean piston speed.

| | |
| --- | --- |
| ![Torque before and after adding the flywheel, with the surplus-energy area shaded](docs/media/torque-flywheel.webp) | ![Crank speed over the cycle with the flywheel, swinging 4% around 2000 rpm](docs/media/speed-fluctuation.webp) |
| The flywheel stores the 369 J of surplus torque. | Crank speed held within 4% of 2000 rpm. |

**6. Gearbox.** The top gear comes from a top speed of 245 km/h, the first gear from what it takes to move a 2600 kg vehicle up a 7% grade while accelerating, and `fsolve` finds the two factors of a progressive series of eight ratios between them.

![Bar chart of the eight gear ratios from 5.38 in first to 0.64 in eighth](docs/media/gear-ratios.webp)

**7. Animation.** `animation.m` draws the four slider-cranks moving together and flashes each cylinder as it fires, which makes the phase shift the torque sum depends on visible. The [engine video](docs/media/engine-simulation.mp4) at the top is the same idea drawn from the saved results: the mechanism to scale, with the gas, flame, valves and charts added.

## Quick start

You need MATLAB. `part2.m` calls `fsolve`, which comes with the Optimization Toolbox.

1. Open MATLAB with the repository root as the current folder.
2. Run the scripts in order:

   ```matlab
   constants   % saves constants.mat
   part1       % kinematics and volume, saves part1.mat
   part2       % cycle, torque, flywheel, performance, gearbox; saves part2.mat
   animation   % the four-cylinder animation
   ```

Each script loads the `.mat` file the previous one saved and opens its own figures. The `.mat` files are committed, so you can also start from any step. Results print to the Command Window.

### In GNU Octave

The scripts also run unchanged in [GNU Octave](https://octave.org) (tested with 10.3), where `fsolve` is built in. From the repository root:

```bash
octave --no-gui tools/octave/run_all.m
```

`run_all.m` runs `constants`, `part1` and `part2` with two small stand-ins on the path for `xline` and `yline`, which Octave lacks. Like the scripts in MATLAB, it saves `constants.mat`, `part1.mat` and `part2.mat` to the current folder, over the committed ones. Octave can't write the function handle `equations` into a MAT file, so `part2.mat` is saved without it. [`docs/octave-run.log`](docs/octave-run.log) is the console output of the run in the video.

Octave lands within 0.1% of MATLAB: 459.3 J of work per cycle against 459.6 J, and 30.62 kW against 30.64 kW. The cycle loop stops when the temperature at the end of intake changes by less than 0.01 K, but it wobbles by a few hundredths of a kelvin from one pass to the next, so MATLAB stops after 4 passes and Octave after 6.

## Figures and videos

Everything in `docs/media` is drawn from the scripts' output by the Python scripts in `tools/`. They read `part2.mat` and never recompute the physics. You need numpy, scipy, matplotlib, Pillow and ffmpeg; the videos use the IBM Plex fonts.

```bash
python3 tools/render_figures.py                                    # the ten figures, as WebP
python3 tools/render_engine_video.py                               # engine-simulation.mp4, .gif and poster
python3 tools/render_run_video.py part2.mat docs/octave-run.log    # simulation-run.mp4 and poster
```

The figures and the engine video use the committed `part2.mat` from MATLAB. The run video uses the `part2.mat` of the Octave run it replays, so its numbers match its console.

## Repository layout

```
.
├── constants.m     every parameter -> constants.mat
├── part1.m         slider-crank kinematics and cylinder volume -> part1.mat
├── part2.m         cycle, torque, flywheel, performance, gearbox -> part2.mat
├── animation.m     four-cylinder animation
├── *.mat           saved workspaces from each step
├── submission/     the same scripts and .mat files, as handed in
├── tools/
│   ├── octave/     run_all.m, plus xline and yline stand-ins for GNU Octave
│   └── *.py        renderers for the figures and videos in docs/media
└── docs/
    ├── media/      figures and videos for this README
    └── octave-run.log   console output of the Octave run
```

## Limits

These are coursework scripts, and they make simplifications on purpose.

- The mean effective pressures the scripts print (MEP and BMEP) are about half what they should be: the displacement constant is smaller than the engine's swept volume of about 503 cc.
- Friction is a flat 10% of indicated power. There is no heat transfer to the cylinder walls, no emissions model and no knock.
- The engine turns at one fixed speed, and the speed-fluctuation curve is an idealised sine between the flywheel's limits, not a simulation of the shaft.

## Links

- Project page: [ramikronbi.com/projects/ice-engine-simulator](https://ramikronbi.com/projects/ice-engine-simulator)
- Write-up: [Simulating an engine from crank to gearbox](https://ramikronbi.com/writing/simulating-an-engine-from-crank-to-gearbox)

Built by Rami Kronbi as the main project of the Automotive Engineering course at Rafik Hariri University (RHU), fall 2024. [ramikronbi.com](https://ramikronbi.com)
