## 🛠️ How to Build

This project uses **CMake** for cross-platform builds.

### Prerequisites

Make sure you have the following installed:

- A modern C and C++ compiler
- CMake 3.1+

### Build Instructions

Repeat the following steps for each tool (`cache-simulator`, `branch-predictor`, or `raw-dependency-analyzer`):

```bash
cd <tool_name>
mkdir build && cd build
cmake ..
make
```
### On windows:
if you do have make on windows you should use:
```bash
mkdir build && cd build
cmake .. -G "MinGW Makefiles"
make
```
otherwise:
```bash
mkdir build && cd build
cmake ..
cmake --build . --config Release
```

exe should now be inside /build/Release folder 
# Computer Architecture Simulation Suite

A modular suite of tools designed to simulate and analyze fundamental components of computer architecture.  
The project includes:

- ✅ Multi-Level Cache Simulator (L1/L2)
- ✅ Branch Prediction Simulator
- ✅ RAW Dependency Analyzer for Out-of-Order Execution

Each tool is developed in C++ and designed for trace-based benchmarking and analysis.

---

## 📦 Project Structure

Simulation-Suite/ \
│ \
├── cache-simulator/ # L1/L2 cache hierarchy simulator\
├── branch-predictor/ # 2-level branch prediction simulator\
├── raw-dependency-analyzer/ # RAW dependency graph analyzer for OoO pipelines

---

## 🧰 Tool Descriptions

### 1. 🧠 Cache Simulator

Simulates a multi-level (L1/L2) inclusive cache system with:

- LRU eviction policy
- Write-back and write-allocate
- Non-blocking write-back handling
- Miss rate and access time statistics

### 2. 🧮 Branch Predictor

Implements a configurable two-level predictor with:

- Direct-mapped BTB
- Bimodal FSM states
- Global/Local (G-Share / L-Share) history schemes
- Trace-based speculative execution and flush-rate analysis

### 3. 🕸️ Dependency Analyzer

Analyzes Read-After-Write (RAW) hazards in Out-of-Order (OoO) execution:

- Instruction-level dependency graph construction
- Critical path detection
- Program depth estimation using latency-aware topological sort
- Resolves false dependencies via register renaming
