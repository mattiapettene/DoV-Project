# Features <!-- omit in toc -->

This file describes the implemented features of the library in detail.

## Table of contents <!-- omit in toc -->

- [Visualize single/multi-body objects in 3D](#visualize-singlemulti-body-objects-in-3d)
  - [Custom objects](#custom-objects)
  - [Custom camera objects](#custom-camera-objects)
- [Visualize 2D/3D tracks](#visualize-2d3d-tracks)
  - [Default tracks](#default-tracks)
- [Animate and/or save a video in post-processing/step-by-step](#animate-andor-save-a-video-in-post-processingstep-by-step)
  - [Post-processing](#post-processing)
  - [Step-by-step](#step-by-step)
- [Deferred update of objects](#deferred-update-of-objects)


## Visualize single/multi-body objects in 3D

The library supports the visualization of both single and multi-body objects. The former are quite self-explanatory, if you add to the scene a single object (e.g., a single STL object) you will get a single-body object. Whilst these objects could be anything (since the library allows you to add STL objects, or even your own), they are limited to the 9 Degrees of Freedom (DoF) of MATLAB's [`hgtransform`](https://it.mathworks.com/help/matlab/creating_plots/transform-objects-within-a-group.html), which are rotation, translation, and scaling. If you desire greater freedom with your objects, you can use the multi-body objects. They are a collection of single-body objects connected together by means of kinship. It is customary in multi-body modelling to define one or more forefathers from which a tree of kinship stems. For example, when modelling a car one could define a moving reference frame as the forefather and then define the chassis and each individual wheel as its children. Defined a car in this way, you will have 9 DoF per body, which amounts to 9×6 = 54 DoF in total. Keep in mind though that when two bodies are related, the child will move relative to its father, this is intuitive but may lead to unexpected results (believe us).

### Custom objects

Whilst the library provides a set of default objects, you can easily add your own. To do so, you need to create a class that inherits from the [`BaseObject`](src/abstract/@BaseObject/BaseObject.m) abstract class. The class must implement the method `out = create_descriptor( obj )`. The library calls this method when first rendering the scene and is responsible for creating the object's descriptor and passing its handle. The descriptor can be any of MATLAB's built-in graphical objects (e.g., scatter3, patch) or something else you come up with (as long as it is attachable to the axes of the animation).

Even though the only method you are required to implement is `create_descriptor( obj )`, you can also override the [`BaseObject`](src/abstract/@BaseObject/BaseObject.m) built-in methods and change how an object behaves. We have done that for the [`ScatterObject`](src/objects/ScatterObject.m), check it out for inspiration.

### Custom camera objects

Although a camera object may seem rather different from an ordinary object, it shares a lot of the same ideas. Even though it is implemented separately, users can still create custom camera objects to suit their needs. To do so, you need to create a class that inherits from the [`BaseCamera`](src/abstract/@BaseCamera/BaseCamera.m) abstract class. The class must implement the method `out = create_state( obj )`. The library calls this method every time it updates the camera position and is responsible for creating the camera's state. Depending on the use case, you may want the camera to stay still (as implemented in the [`FixedCamera`](src/cameras/FixedCamera.m) class), follow a moving object (as implemented in the [`FollowerCamera`](src/cameras/FollowerCamera.m) class), or even have its own independent motion. You can easily implement your own camera class by following the examples of the built-in ones.

## Visualize 2D/3D tracks

The library was born with 3D in mind, so when we talk about a 2D track what we mean is a flat track visualized in 3D. We provide a set of default tracks, but you can easily add your own. To do so, you need to create a class that inherits from the [`BaseTrack`](src/abstract/BaseTrack.m) abstract class. The class must implement the method `plot( obj, ax )` and assign the property `ax` which sets the axes of the plot. The library calls this method when first rendering the scene and is responsible plotting the track in the scene. We do not enforce any restrictions on the track's plot, since it is up to the user to decide how to visualize it.

### Default tracks

The library provides a set of default tracks, which are:

- [`RaceTrack`](src/tracks/RaceTrack.m) represents a 2D racetrack which is defined by its left and right margins. We plot the track as a 3D surface (it is actually a patch for MATLAB's enthusiasts). (IT MIGHT BE CAPABLE OF 3D TRACKS, BUT WE HAVE NOT TESTED IT YET)
- [`RDFTrack`](src/tracks/RDFTrack.m) represents a 3D track which is defined by an RDF file (i.e., a file containing a set of faces and vertices in 3D). We plot the track as a 3D surface (it is actually a patch for MATLAB's enthusiasts).
- [`RingTrack`](src/tracks/RingTrack.m) represents a 3D testing track which is defined by its radius and number of rings. We construct the track as a number of concentric rings that can have a certain banking (specified as height of the centreline). We plot the track as a 3D surface (it is actually a patch for MATLAB's enthusiasts). (UNDER ACTIVE DEVELOPMENT)

## Animate and/or save a video in post-processing/step-by-step

The library offers two options to animate your scenario: in post-processing or step-by-step. By default the library works in the post-processing mode, expecting the user to provide all the data at once and then rendering the animation. However, it might also be useful to generate the animation during a simulation (i.e., in “real-time”). For this purpose, the library has some tools that enable you to attach the live animation to your simulation (although it will require some more effort from the user side). In both modes the library is capable of saving the animation as a video file, but the quality of the video (e.g., resolution, frame rate) is guaranteed only in the post-processing mode.

The reason the post-processing is the default mode is mainly due to the poor performance that MATLAB has with animations (especially when saving a video at the same time). Imagine that even in this mode, where your machine is probably doing nothing else, the library is still not able to guarantee a smooth animation (although we try our best). This is why we do not recommend using the library in the step-by-step mode, unless you have a very powerful machine or you do not care if the animation is not smooth.

### Post-processing

This is the easiest way to animate your scenario and the default one. To animate your scenario, all you need to do is simply create the objects in the scene with the entire dataset from the start (e.g., the entire state evolution) and animate it. You can find examples of this in the [`examples`](../examples) folder.

### Step-by-step

This is the other visualization mode we provide for the sole purpose of animating your scenario one step at a time, in “real-time”. This proves useful if you need to visualize your scenario during a simulation for example. To animate your scenario, you need to create the objects in the scene with the initial state and then update them at each time step. You can find examples of this in the [`examples`](../examples) folder.

## Deferred update of objects

This is a fairly advanced feature that allows you to defer the update of objects only when it is needed. By default, the library already renders the scene at a specific frame-rate (be it the default one or a user-specified one), which already improves performance considerably. In fact, we discourage to use this feature to improve performance. Before we explain what it is, let us first explain how the library updates the scene.

- First, the library renders the very first frame of the animation.
- Then, the library starts a time-based loop from the start time $t = 0$ to the end time specified by the user $t = t_f$ with a step size equal to the sampling time $t_s$.
- At each time step, the library updates the objects (e.g., the objects' state) in the scene but does **not** render the scene.
- When the time step is a multiple of the frame-rate, the library renders the scene.

A deferred update will only changes when the library updates an object, not when the library renders the scene. Given that updating objects usually takes little to no time, this feature may not improve performance at all (hence our discouragement). However, if you need to update the position of an object only at specific multiples of the time step, then this feature enables you to do just that. For example, assume you have the following scene:

- $t_s = 0.001\;s$;
- a car object that needs to be updated at every time step;
- a planned trajectory object that needs to be updated only every $10$ time steps (i.e., every $0.01\;s$).

You can achieve this scenario very easily by using the deferred update feature. You can find some examples of this in the [`examples`](../examples) folder.
