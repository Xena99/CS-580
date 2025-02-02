# CS580 Homework - Rendering Utah Teapot Implementation

<table>
  <tr>
    <td><img src="./images/Phong%20Shading.png" alt="Phong Shading" width="350"/></td>
    <td><img src="./images/Texture%20Mapping.png" alt="Texture Mapping" width="350"/></td>
  </tr>
</table>

## Overview
This project is part of CS580, where a simple 3D rendering pipeline is implemented to render triangles using various transformations, Phong shading, texture mapping, and depth buffering (Z-buffer). The main features of this project include:

- **Matrix-based transformations:** for world, camera, and projection transformations.
- **Phong shading:** for realistic lighting effects.
- **Texture mapping:** to map textures onto 3D surfaces.
- **Camera transformations:** to handle view changes using transformation stacks.

---

## Key Components

### 1. **Matrix Transformations**
The rendering pipeline relies on matrix operations to transform vertices through different spaces:

- **Model Space**: Local object coordinates.
- **World Space**: Global scene coordinates.
- **Camera/View Space**: Coordinates relative to the camera.
- **Screen Space**: 2D coordinates mapped to the display.

The following transformation matrices are used:
- **Rotation Matrices (X, Y, Z)**: For rotating objects around the respective axes.
- **Translation Matrix**: For moving objects within the scene.
- **Scaling Matrix**: For resizing objects.
- **Perspective Projection Matrix**: To simulate the perspective effect.

These matrices are combined using **matrix multiplication** and stored in a transformation stack.

---

### 2. **Phong Shading Model**
Phong shading is used to compute pixel colors using the following components:

- **Ambient Reflection:** Simulates indirect light.
- **Diffuse Reflection:** Based on the angle between the surface normal and light direction (Lambertian reflection).
- **Specular Reflection:** Simulates shiny surfaces, using the reflection direction and view direction.

#### Formula:
The final color \( C \) for a pixel is computed as:

\[
C = k_a \cdot I_a + \sum_{i=1}^{n} \left[ k_d (N \cdot L_i) I_i + k_s (R_i \cdot V)^n I_i \right]
\]

Where:
- \( k_a, k_d, k_s \) are the material coefficients for ambient, diffuse, and specular reflections.
- \( I_a \) is the ambient light intensity.
- \( I_i \) is the intensity of the directional light.
- \( N \) is the surface normal.
- \( L_i \) is the light direction.
- \( R_i \) is the reflection vector.
- \( V \) is the view direction.
- \( n \) is the shininess coefficient (specular power).

---

### 3. **Texture Mapping**
Textures are applied to objects using a mapping from the object’s surface to a 2D texture image.

**Steps:**
1. **Texture Coordinates (u, v)**: Provided per vertex and interpolated across the triangle using barycentric interpolation.
2. **Perspective Correction**: The coordinates are adjusted using the perspective projection to avoid distortion.
3. **Bilinear Interpolation**: To sample the texture color smoothly between adjacent pixels.

#### Texture Function Implementation (`tex_fun`)
- Reads texture data from a binary PPM file.
- Performs bilinear interpolation to get the texture color.

---

### 4. **Camera Transformation**
The camera transformation matrix \( X_{iw} \) is constructed to transform world coordinates into camera coordinates.

**Steps to Compute \( X_{iw} \):**
1. Calculate the **Z-axis** of the camera space: \( Z = \text{normalize}(\text{lookat} - \text{position}) \).
2. Calculate the **up' vector**: \( \text{up'} = \text{up} - (\text{up} \cdot Z) Z \).
3. Calculate the **X-axis** using the cross product: \( X = \text{normalize}(\text{up'} \times Z) \).
4. Calculate the **Y-axis**: \( Y = Z \times X \).
5. Fill the rotation and translation components of \( X_{iw} \).

---

### 5. **Rendering Pipeline Workflow**
1. **Vertex Processing**: Transform object vertices through the model, view, and projection matrices.
2. **Rasterization**: Convert triangles into pixel fragments.
3. **Depth Testing (Z-buffer)**: Ensure that only the closest fragments are rendered.
4. **Shading**: Compute the pixel color using the Phong shading model.
5. **Texture Mapping**: Apply texture colors to the fragments.

---

## Major Functions

### **Matrix Operations**
- `GzRotXMat(float degree, GzMatrix mat)`: Creates a rotation matrix around the X-axis.
- `GzRotYMat(float degree, GzMatrix mat)`: Creates a rotation matrix around the Y-axis.
- `GzRotZMat(float degree, GzMatrix mat)`: Creates a rotation matrix around the Z-axis.
- `GzTrxMat(GzCoord translate, GzMatrix mat)`: Creates a translation matrix.
- `GzScaleMat(GzCoord scale, GzMatrix mat)`: Creates a scaling matrix.

### **Camera Setup**
- `GzPutCamera(GzRender *render, GzCamera *camera)`: Initializes the camera and its transformation matrix.
- `ComputeXiw(GzMatrix& Xiw, const GzCamera& camera)`: Computes the world-to-camera transformation.

### **Shading and Lighting**
- `GzShadingEquation(GzRender *render, GzColor& color, const GzCoord normal)`: Computes the color for a fragment using Phong shading.

### **Rasterization**
- `GzPutTriangle(GzRender *render, int numParts, GzToken *nameList, GzPointer *valueList)`: Rasterizes a triangle and computes the final pixel colors using texture mapping and shading.
- `sortVerticesAndNormals`: Sorts the vertices of a triangle in counter-clockwise order.

### **Texture Handling**
- `tex_fun(float u, float v, GzColor color)`: Performs texture lookup and bilinear interpolation.

---

## Data Structures
- **`GzRender`**: Holds rendering state, including the transformation matrices, lighting parameters, and Z-buffer.
- **`GzCamera`**: Defines the camera's position, orientation, and field of view.
- **`GzDisplay`**: Manages the pixel data for the final rendered image.

---

## How to Run the Program
1. **Compile** the project using a compatible C++ compiler.
2. **Run** the program with a scene configuration to load models, set lighting parameters, and start rendering.
3. **Output**: The rendered image will be displayed or saved as a PPM file.
