# Objects

This directory contains the various types of objects that can be represented in a scene.

All the objects inherit from the [BaseObject](../abstract/BaseObject/BaseObject.m) class, which can be used by the user as a template for creating new objects and extend the existing functionality.

## Object classes

- [STLObject](STLObject.m) - Object for STL objects. It renders a patch object with the STL data.
- [ScatterObject](ScatterObject.m) - Object for scatter objects. It renders a scatter object with the provided data.