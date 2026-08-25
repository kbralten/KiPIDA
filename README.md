# Ki-PIDA (KiCad Power Integrity & Delivery Analyzer)

Ki-PIDA is a native KiCad plugin for DC, AC, and thermal Power Integrity (PI) analysis. It allows PCB designers to simulate voltage drops (IR drop), current densities, rail-to-ground impedance, and steady-state 3D board temperatures directly within the KiCad Pcbnew environment, eliminating complex external workflows.

## 🚀 Why Ki-PIDA?

Modern electronics operate with tight voltage margins. An IR drop of just 30mV can lead to system instability in sub-1.0V Socs. High current densities also pose thermal risks and reliability hazards like electromigration. 

Ki-PIDA democratizes high-end PI analysis by:
- **Ensuring Stability:** Detect voltage violations at the layout stage.
- **Reducing Iterations:** Identify "neck-down" regions and hotspots before prototyping.
- **Seamless Workflow:** Interactive layout-driven analysis without leaving KiCad.

## ✨ Key Features

- **Native Integration:** Built for KiCad 9.0+ using the Python Scripting API.
- **Power Tree Management:** Auto-discover power rails and manage complex hierarchies including VRM efficiency modeling.
- **Hybrid 2.5D Solver:** Fast and accurate simulation using an optimized resistive mesh approach.
- **Multi-Physics Support:** Coupled electro-thermal simulation to account for temperature-dependent copper resistivity.
- **Multi-Rail Analysis:** Simulate complex power trees with nested regulators (Buck, LDO) and enforce correct dependency solving.
- **Project Persistence:** Automatically saves your power tree configuration (sources, loads, regulators) in the project directory, so you don't lose your setup.
- **Visual Feedback:** Interactive heatmaps for voltage and current density, with dedicated tabs for each power rail in the system.
- **AC Impedance:** Sweep rail-to-ground impedance magnitude and phase over a logarithmic frequency range.
- **Decoupling Optimization:** Rank values for existing unpopulated/DNP capacitor footprints against a target-impedance envelope.
- **3D Thermal Model:** Solve through-stack and lateral heat conduction using the extracted PCB stackup, spatial copper coverage, and thermal vias.
- **Airflow Convection:** Configure natural, forced, or custom heat-transfer coefficients with exposed top, bottom, and edge surfaces.
- **Electro-Thermal Coupling:** Iterate DC copper loss and temperature-dependent copper resistance until the configured convergence threshold is reached.

## 📦 Installation

Ki-PIDA is designed to run within the KiCad 9.0+ environment. Follow these steps to install and enable the plugin:

### 1. Enable the KiCad API
Ki-PIDA communicates with KiCad via the new IPC-based API.
1. Open KiCad.
2. Go to **Preferences** > **Common**.
3. Under the **API** section, check the box for **Enable API**.
4. Restart KiCad if prompted.
![alt text](image-3.png)

### 2. Install the Plugin
1. Locate your KiCad plugins directory:
   - **Windows:** `%APPDATA%\kicad\9.0\plugins`
   - **Linux:** `~/.local/share/kicad/9.0/plugins`
   - **macOS:** `~/Library/Application Support/kicad/9.0/plugins`
2. Download or clone this repository.
3. Copy the `KiPIDA` folder into the `plugins` directory.

> [!NOTE]
> Ki-PIDA includes a self-contained dependency manager that will automatically install required libraries (NumPy, SciPy, Shapely, Matplotlib) upon first launch if they are missing from your KiCad Python environment.

---

## 📖 Tutorial: Your First IR Drop Analysis

Follow these steps to perform a DC Power Integrity analysis on your board.

### 1. Launch the Plugin
Open your PCB layout in KiCad Pcbnew and click the **Ki-PIDA** icon in the top toolbar to open the analyzer.

### 2. Review Discovered Power Rails
- **Auto-Discovery:** Upon launch, Ki-PIDA scans your board for power rails and attempts to load any existing configuration from `kipida_config.json` in your project folder.
- **Add Roots:** Identify your main input rails (e.g., `+12V_IN`, `VBUS`).
- **Define Regulators:** Use the **+ Regulator** button to create relationships between rails (e.g., `12V -> 5V`). Ki-PIDA supports:
    - **Linear Regulators (LDOs):** Pass current 1:1 from input to output.
    - **Switching Regulators (Buck/Boost):** Conserve power based on efficiency (e.g., 90%).
    - **Multi-Output Support:** Handle PMICs where one component drives multiple output rails.

![alt text](image-4.png)

### 3. Add Sources (VRMs / Power Inputs)
Identify where power enters this net:
1. Click **+ Source**.
2. Select the source component (e.g., a regulator `U1` or connector `J1`).
3. In the dialog, check the **Pads** that are connected to the power net.
4. Click **OK**.

### 4. Add Loads (Integrated Circuits / Sinks)
Identify the components consuming power:
1. Click **+ Load**.
2. Select the sink component (e.g., MCU `U2` or FPGA `U3`).
3. Enter the **Total Current (A)** consumed by this component (e.g., `0.5` for 500mA).
4. Check the **Pads** through which the current is drawn.
5. Click **OK**.

![alt text](image-5.png)

### 5. Run the Simulation
Before running, you can adjust the **Mesh Resolution (mm)**. A value of `0.1mm` is usually sufficient for accurate results.
- Click **Run Simulation**.
- The solver analyzes the power tree topology to determine the correct solution order (Leaf-to-Root) for current propagation.
- Example: For a `12V -> 5V -> 3.3V` chain, it solves 3.3V first, applies that load to the 5V rail, solves 5V, and finally solves the 12V input.

### 6. Analyze Results
Once "Simulation Success" appears, the UI will jump to the **Results** tab.

- **Rail Selection:** A tab will be created for each power rail in your system.
- **Per-Rail Visualization:** Inside each rail's tab, you can view:
    - **3D View:** A 3D voltage plot of the entire net.
    - **Layer Views:** Individual 2D heatmaps for every layer containing copper for that net.

> [!TIP]
> Use the **Enable Debug Log** checkbox if you encounter issues during meshing or solving to see more detail in the Log tab.

## Tutorial: AC Impedance and Decoupling Optimization

1. Define the rail's sources and loads in **Power Tree & Config**, then open **AC Impedance**.
2. Select the power rail, return net, source component, and measurement component. Regulator outputs that feed the selected rail are offered as source components.
3. Set the logarithmic sweep range, point count, target impedance, and the source's small-signal resistance/inductance.
4. Review the detected rail-to-ground capacitors. Ki-PIDA reads common KiCad values such as `100n`, `4u7`, and `10uF`; package-derived ESR/ESL values remain engineering estimates.
5. Choose **Run AC Analysis** to plot impedance magnitude and phase, or **Optimize Decoupling** to evaluate the detected disabled/DNP capacitor footprints.
6. Save the project configuration to persist the AC profile in `<project>.kipida.json` beside the `.kicad_pro` file.

The optimizer is intentionally non-destructive: it reports footprint/value recommendations but does not modify the PCB or schematic. It only uses candidate capacitor locations already present on the board.

## Tutorial: 3D Thermal and Airflow Analysis

1. Define rail voltages, loads, and regulator efficiencies in **Power Tree & Config**, then open **3D Thermal**.
2. Select the ambient temperature, thermal grid size, exposed surfaces, and an airflow mode:
   - **Natural:** uses a conservative natural-convection coefficient.
   - **Forced:** derives the coefficient from air speed and applies the configured flow direction across the board surface.
   - **Custom:** uses a user-supplied heat-transfer coefficient.
3. Choose **Refresh Power Estimates**. Load dissipation is estimated as `V × I`; regulator dissipation uses LDO voltage drop or switching efficiency. Double-click any component to enter a reviewed power and compact package thermal model.
4. Choose **Run Thermal** for a single steady-state solve. Enable **Include DC copper losses** to reuse losses from the DC branch solution.
5. Choose **Run Coupled** to iterate copper resistance, DC branch loss, and board temperature. Review the hotspot, energy balance, component junction estimates, and 3D/top/bottom plots in **Results**.
6. Save the project configuration to persist the thermal profile in `<project>.kipida.json`.

The airflow model applies convective boundary conditions to the 3D solid board mesh. It is intended for board-level design comparison and hotspot screening; it is not a volumetric CFD enclosure or fan model. Component junction temperatures use the configured compact `theta-JB` estimate and therefore require engineering review before sign-off.

## 🛠️ Technical Overview (For Developers)

Ki-PIDA is built on a modular architecture designed for performance and maintainability.

### Architecture
- **Extractor (`extractor.py`):** Interfaces with the KiCad API to pull filled zone geometry, track layouts, and physical stackup data.
- **Mesher (`mesh.py`):** Discretizes continuous copper geometry into a 2D/3D resistive grid (Rasterization).
- **Solver (`solver.py`):** Uses an Admittance Matrix (Stamps method) and optimized SciPy sparse solvers (SuperLU/CG) to solve the electrical system.
- **AC Model (`ac_model.py`):** Builds coupled rail/return meshes and maps sources, measurement ports, and rail-to-ground capacitors.
- **AC Solver (`ac_solver.py`):** Stamps frequency-dependent sparse complex admittance matrices for copper/via RL branches and capacitor RLC models.
- **Decoupling Optimizer (`decoupling_optimizer.py`):** Deterministically searches existing candidate footprints against the target-impedance score.
- **Thermal Model (`thermal_model.py`):** Combines stackup, board copper, vias, component placement, power-tree dissipation, and optional DC branch losses.
- **Thermal Mesh (`thermal_mesh.py`):** Builds a finite-volume 3D solid mesh with anisotropic FR-4, spatial copper conductivity, thermal-via branches, radiation, and convective surface boundaries.
- **Thermal Solver (`thermal_solver.py`):** Solves the sparse steady-state heat equation and reports hotspot, component junction estimates, and energy balance.
- **Electro-Thermal Solver (`electrothermal.py`):** Iterates DC branch resistance and loss with the solved copper temperature field.
- **Visualizer (`visualizer.py`):** Generates heatmaps via Matplotlib and renders them as overlays in KiCad.

### Methodology
Electrical analysis utilizes a **Hybrid 2.5D Finite Difference Method (FDM)**. It represents PCB layers as 2D grids connected vertically by via/PTH elements. DC analysis uses resistive branches; AC analysis retains the same topology and adds stackup-sensitive branch inductance plus lumped source/capacitor RLC models. Thermal analysis uses a separate **3D finite-volume solid-conduction model** through the physical stackup. This is a board-level engineering model, not a full-wave electromagnetic solver or volumetric CFD solver.

### Stack
- **Languages:** Python 3.9+
- **UI:** wxPython
- **Math:** NumPy & SciPy
- **Geometry:** Shapely

## � Current State (Alpha)

As of the current version, Ki-PIDA implements end-to-end DC IR drop, AC target-impedance, and steady-state 3D thermal analysis.

### Capabilities:
- **Comprehensive Extraction:** Extracts tracks, pads, and filled zones (respecting thermal reliefs and voids) from KiCad 9.0+ boards.
- **3D Meshing Engine:** Converts geometry into a resistive mesh across multiple layers, correctly modeling via and PTH conductances.
- **Robust Linear Solver:** Solves the circuit using SciPy's sparse matrix backend. Includes island detection to warn about floating sections of copper that could cause numerical issues.
- **Automated Diagnostics:** Detects isolated copper nodes and disjoint electrical islands during the solve phase.
- **Target-Impedance Sweep:** Reports worst-case impedance/frequency and PASS/FAIL against the configured envelope.
- **Deterministic Decoupling Search:** Recommends values for existing DNP/candidate footprints without editing the board.
- **Thermal and Airflow Solve:** Models 3D board conduction, natural/forced/custom convection, radiation, vias, and component heat sources.
- **Coupled Iteration:** Feeds branch-level DC `I²R` losses into the thermal solve and updates copper resistance with temperature.

### User Experience:
- **Automated Rail Discovery:** Instantly find power nets based on zone connectivity.
- **Granular Control:** Assign sources and loads down to the individual pad level.
- **In-Memory Visualization:** Instant generation of color-coded heatmaps to inspect voltage distribution without exporting files.

## �🗺️ Roadmap

- **Phase 1:** DC IR Drop, basic thermal checks, and power tree UI.
- **Phase 2:** AC Impedance Analysis ($Z$ vs Frequency) and decoupling capacitor optimization.
- **Phase 3 (Current):** Full 3D board thermal modeling with airflow convection and iterative DC coupling.
