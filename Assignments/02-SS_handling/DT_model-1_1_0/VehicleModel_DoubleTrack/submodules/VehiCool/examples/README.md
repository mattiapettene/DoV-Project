# Examples

This directory contains detailed examples of how to use the `VehiCool` library.

## Test data

The directory [`data`](data) contains a set of test data files that can be used to test the library. The files are:

- [`adria.mat`](data/adria.mat) contains the left and right margins for the Adria racetrack;
- [`state_rf0.mat`](data/state_rf0.mat) contains the state of the main reference frame for the vehicle;
- [`state_chassis.mat`](data/state_chassis.mat) contains the state of the chassis;
- [`state_wheel_f.mat`](data/state_wheel_f.mat) contains the state of the front wheel (we assume that the left and right wheels are identical);
- [`state_wheel_r.mat`](data/state_wheel_r.mat) contains the state of the rear wheel (we assume that the left and right wheels are identical);
- [`trajectory.mat`](data/trajectory.mat) contains the trajectory of the vehicle;
- [`time_sim.mat`](data/time_sim.mat) contains the time vector for the simulation.

## Basic usage

The directory [`basic_usage`](basic_usage) contains a set of basic examples of how to use the library.

### Example 1: single-body object with a following camera in post-processing

The file [`example_1.m`](basic_usage/example_1.m) shows how to use the library to create a simple object and animate it in post-processing.

### Example 2: single-body object with a following camera in real-time

The file [`example_2.m`](basic_usage/example_2.m) shows how to use the library to create a simple object and animate it in “real-time” (i.e., one step at a time, as if it was attached to a simulation).

### Example 3: single-body object with a trajectory and a following camera in post-processing

The file [`example_3.m`](basic_usage/example_3.m) shows how to use the library to create a simple object and animate it in post-processing, while plotting a static trajectory.

### Example 4: single-body object with a fixed camera in post-processing

The file [`example_4.m`](basic_usage/example_4.m) shows how to use the library to create a simple object and animate it in post-processing, while using a fixed camera.

### Example 5: single-body object with a trajectory (with deferred updates) and a following camera in post-processing

The file [`example_5.m`](basic_usage/example_5.m) shows how to use the library to create a simple object and animate it in post-processing, while plotting a dynamic trajectory that is updated asynchronously.

### Example 6: single-body object with a trajectory (with deferred updates) and a birdeye camera in post-processing

The file [`example_6.m`](basic_usage/example_6.m) shows how to use the library to create a simple object and animate it in post-processing, while plotting a dynamic trajectory that is updated asynchronously and using a birdeye camera.

## Advanced usage

The directory [`advanced_usage`](advanced_usage) contains a set of advanced examples of how to use the library.

### Example 1: multi-body system with a following camera in post-processing

The file [`example_1.m`](advanced_usage/example_1.m) shows how to use the library to create a simple multi-body system and animate it in post-processing.

### Example 2: multi-body system with a following camera in real-time

The file [`example_2.m`](advanced_usage/example_2.m) shows how to use the library to create a simple multi-body system and animate it in “real-time” (i.e., one step at a time, as if it was attached to a simulation).
