# Cameras

This directory contains the various types of camera objects that can be used in a scene.

All the cameras inherit from the [BaseCamera](../abstract/@BaseCamera/BaseCamera.m) class, which can be used by the user as a template for creating new camera objects and extend the existing functionality.

## Camera classes

- [FixedCamera](FixedCamera.m) - Camera object with a fixed position.
- [FollowingCamera](FollowingCamera.m) - Camera object that follows a moving object.