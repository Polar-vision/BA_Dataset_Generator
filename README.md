# 🚀 Building from Source
**BA_Dataset_Generator:** This tool is designed to convert camera and reconstruction results exported from **COLMAP**
(`cameras.txt`, `images.txt`, `points3D.txt`) into a unified dataset format suitable for **Bundle Adjustment (BA)** and related optimization pipelines.

## 📂 Data Format  
Each dataset is organized using the following structure: 
### 📐 Intrinsics (`cal.txt`)

Contains the camera intrinsic parameters for each camera:

- **fx, fy** — focal lengths  
- **cx, cy** — principal point coordinates  

Each entry corresponds to one camera.

---

### 🎥 Extrinsics (`Cam.txt`)

Contains the camera pose parameters, including:

- **Euler angles**
  - **ez** — rotation around the z-axis  
  - **ey** — rotation around the y-axis  
  - **ex** — rotation around the x-axis  
- **Perspective center**
  - **Xc, Yc, Zc**
- **Camera ID**

---

### 🌍 3D Points (`XYZ.txt`)

Lists the reconstructed 3D object points:

- **X, Y, Z** coordinates for each point.

---

### 🔗 Feature Tracks (`Feature.txt`)

Each line represents a feature track observed across multiple images, including:

- Number of observing views  
- Corresponding image indices  
- Image coordinates **(u, v)** in each image  

This file encodes the 2D–3D observation relationships required for Bundle Adjustment.

## 🖥 Tested Platforms
This project has been tested on:

- **Operating System:** Windows 11  
- **Compiler / IDE:** Visual Studio 2022 (Win32)  
- **CMake version:** ≥ 3.10  

> Note: Other platforms may work, but they have not been verified.


## 📦 Dependencies
- **OpenCV**: 4.11.0   
- **Eigen**: 3.3.5
---

## ⚙️ Installation

### 1️⃣ Clone the repository
```bash
git clone https://github.com/Polar-vision/BA_Dataset_Generator.git
cd BA_Dataset_Generator
```

### 2️⃣ Build the library
```bash
mkdir build
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DOpenCV_DIR=<path_to_opencv>/build
cmake --build build --config Release
```

### 3️⃣ Run the demo
```bash
cd build/Release
colmap2pba.exe
```

💡 Note:
Before building and running the demo, make sure that OpenCV has been successfully built and installed.  
In addition, please add the following directory to your system PATH environment variable:  

```bash
<path_to_opencv>/build/x64/vc64/bin
```

This ensures that the required OpenCV runtime DLLs can be found at execution time.
