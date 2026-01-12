# 🚀 Building from Source
**BA_Dataset_Generator:** A flexible and extensible tool for generating synthetic and real-world  datasets for Bundle Adjustment (BA) research and evaluation.

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
example.exe
```

💡 Note:
Before building and running the demo, make sure that OpenCV has been successfully built and installed. 
In addition, please add the following directory to your system PATH environment variable: 

```bash
<path_to_opencv>/build/x64/vc64/bin
```bash

This ensures that the required OpenCV runtime DLLs can be found at execution time.
