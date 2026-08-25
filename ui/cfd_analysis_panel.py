"""wxPython configuration panel for Phase 4 enclosure CFD."""

import math

import wx

try:
    from models import CFDBoundaryPatch, EnclosureCFDSettings
except (ImportError, ValueError):
    from models import CFDBoundaryPatch, EnclosureCFDSettings


class CFDBoundaryDialog(wx.Dialog):
    def __init__(self, parent, patch=None):
        super().__init__(parent, title="Enclosure Boundary Patch")
        patch = patch or CFDBoundaryPatch(name="Patch")
        grid = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        grid.AddGrowableCol(1, 1)
        self.name = wx.TextCtrl(self, value=patch.name)
        self.kind = wx.Choice(self, choices=["INLET", "OUTLET", "FAN", "VENT", "WALL"])
        self.kind.SetStringSelection(patch.kind)
        self.face = wx.Choice(self, choices=["XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"])
        self.face.SetStringSelection(patch.face)
        self.center_u = wx.TextCtrl(self, value=f"{patch.center_u:g}")
        self.center_v = wx.TextCtrl(self, value=f"{patch.center_v:g}")
        self.size_u = wx.TextCtrl(self, value=f"{patch.size_u:g}")
        self.size_v = wx.TextCtrl(self, value=f"{patch.size_v:g}")
        self.velocity = wx.TextCtrl(self, value=f"{patch.velocity_m_s:g}")
        self.temperature = wx.TextCtrl(self, value=f"{patch.temperature_c:g}")
        self.pressure = wx.TextCtrl(self, value=f"{patch.pressure_pa:g}")
        for label, control in (
            ("Name:", self.name), ("Type:", self.kind), ("Face:", self.face),
            ("Center U (0..1):", self.center_u), ("Center V (0..1):", self.center_v),
            ("Size U (0..1):", self.size_u), ("Size V (0..1):", self.size_v),
            ("Velocity (m/s):", self.velocity), ("Temperature (C):", self.temperature),
            ("Gauge pressure (Pa):", self.pressure),
        ):
            grid.Add(wx.StaticText(self, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(control, 1, wx.EXPAND)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(grid, 1, wx.EXPAND | wx.ALL, 12)
        sizer.Add(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL), 0, wx.EXPAND | wx.ALL, 12)
        self.SetSizerAndFit(sizer)

    def get_patch(self):
        patch = CFDBoundaryPatch(
            name=self.name.GetValue().strip() or "Patch",
            kind=self.kind.GetStringSelection() or "VENT",
            face=self.face.GetStringSelection() or "XMIN",
            center_u=float(self.center_u.GetValue()), center_v=float(self.center_v.GetValue()),
            size_u=float(self.size_u.GetValue()), size_v=float(self.size_v.GetValue()),
            velocity_m_s=float(self.velocity.GetValue()),
            temperature_c=float(self.temperature.GetValue()),
            pressure_pa=float(self.pressure.GetValue()),
        )
        if any(value < 0 or value > 1 for value in (
            patch.center_u, patch.center_v, patch.size_u, patch.size_v
        )):
            raise ValueError("Patch center and size values must be between 0 and 1.")
        if patch.kind in {"INLET", "FAN"} and patch.velocity_m_s <= 0:
            raise ValueError("An inlet or fan patch requires a positive velocity.")
        return patch


class CFDAnalysisPanel(wx.Panel):
    def __init__(self, parent, log_callback=None):
        super().__init__(parent)
        self.log_callback = log_callback
        self.settings = EnclosureCFDSettings()
        self._init_ui()

    def _init_ui(self):
        main = wx.BoxSizer(wx.VERTICAL)
        geometry_box = wx.StaticBoxSizer(wx.VERTICAL, self, "Enclosure and PCB")
        parent = geometry_box.GetStaticBox()
        grid = wx.FlexGridSizer(cols=6, hgap=6, vgap=6)
        for column in (1, 3, 5):
            grid.AddGrowableCol(column, 1)
        self.width = wx.TextCtrl(parent, value="120")
        self.depth = wx.TextCtrl(parent, value="100")
        self.height = wx.TextCtrl(parent, value="50")
        self.orientation = wx.Choice(parent, choices=["XY", "XZ", "YZ"])
        self.orientation.SetSelection(0)
        self.offset_x = wx.TextCtrl(parent, value="0")
        self.offset_y = wx.TextCtrl(parent, value="0")
        self.offset_z = wx.TextCtrl(parent, value="15")
        self.wall_h = wx.TextCtrl(parent, value="5")
        rows = [
            ("Width (mm):", self.width, "Depth (mm):", self.depth, "Height (mm):", self.height),
            ("Orientation:", self.orientation, "Offset X (mm):", self.offset_x, "Offset Y (mm):", self.offset_y),
            ("Board Z (mm):", self.offset_z, "Wall h (W/m2K):", self.wall_h, "", wx.StaticText(parent)),
        ]
        for row in rows:
            for index in range(0, 6, 2):
                grid.Add(wx.StaticText(parent, label=row[index]), 0, wx.ALIGN_CENTER_VERTICAL)
                grid.Add(row[index + 1], 1, wx.EXPAND)
        geometry_box.Add(grid, 0, wx.EXPAND | wx.ALL, 8)
        main.Add(geometry_box, 0, wx.EXPAND | wx.ALL, 5)

        solver_box = wx.StaticBoxSizer(wx.VERTICAL, self, "Volumetric CFD Solver")
        solver_parent = solver_box.GetStaticBox()
        solver_grid = wx.FlexGridSizer(cols=6, hgap=6, vgap=6)
        for column in (1, 3, 5):
            solver_grid.AddGrowableCol(column, 1)
        self.ambient = wx.TextCtrl(solver_parent, value="25")
        self.cell = wx.TextCtrl(solver_parent, value="5")
        self.iterations = wx.TextCtrl(solver_parent, value="250")
        self.tolerance = wx.TextCtrl(solver_parent, value="1e-4")
        self.relaxation = wx.TextCtrl(solver_parent, value="0.45")
        self.time_step = wx.TextCtrl(solver_parent, value="0.02")
        self.pressure_iterations = wx.TextCtrl(solver_parent, value="60")
        self.buoyancy = wx.CheckBox(solver_parent, label="Boussinesq buoyancy")
        self.buoyancy.SetValue(True)
        self.phase3_heat = wx.CheckBox(solver_parent, label="Use Phase 3 heat sources")
        self.phase3_heat.SetValue(True)
        self.dc_losses = wx.CheckBox(solver_parent, label="Include DC copper losses")
        self.dc_losses.SetValue(True)
        for row in (
            ("Ambient (C):", self.ambient, "Cell (mm):", self.cell, "Iterations:", self.iterations),
            ("Tolerance:", self.tolerance, "Relaxation:", self.relaxation, "Pseudo dt (s):", self.time_step),
            ("Pressure iterations:", self.pressure_iterations, "", self.buoyancy, "", self.phase3_heat),
        ):
            for index in range(0, 6, 2):
                solver_grid.Add(wx.StaticText(solver_parent, label=row[index]), 0, wx.ALIGN_CENTER_VERTICAL)
                solver_grid.Add(row[index + 1], 1, wx.EXPAND)
        solver_box.Add(solver_grid, 0, wx.EXPAND | wx.ALL, 8)
        solver_box.Add(self.dc_losses, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
        self.estimate = wx.StaticText(solver_parent, label="Estimated cells: --")
        solver_box.Add(self.estimate, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
        main.Add(solver_box, 0, wx.EXPAND | wx.ALL, 5)

        patch_box = wx.StaticBoxSizer(wx.VERTICAL, self, "Inlets, Outlets, Fans, and Vents")
        patch_parent = patch_box.GetStaticBox()
        self.patch_list = wx.ListCtrl(patch_parent, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        for index, (label, width) in enumerate((
            ("Name", 110), ("Type", 75), ("Face", 70), ("Center", 100),
            ("Size", 100), ("Velocity", 80), ("Temp", 70),
        )):
            self.patch_list.InsertColumn(index, label, width=width)
        patch_box.Add(self.patch_list, 1, wx.EXPAND | wx.ALL, 5)
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        self.add_patch = wx.Button(patch_parent, label="Add Patch")
        self.edit_patch = wx.Button(patch_parent, label="Edit Patch")
        self.remove_patch = wx.Button(patch_parent, label="Remove Patch")
        self.default_pair = wx.Button(patch_parent, label="Add Fan Pair")
        for button in (self.add_patch, self.edit_patch, self.remove_patch, self.default_pair):
            buttons.Add(button, 0, wx.RIGHT, 5)
        patch_box.Add(buttons, 0, wx.ALL, 5)
        main.Add(patch_box, 1, wx.EXPAND | wx.ALL, 5)
        main.Add(wx.StaticText(self, label=(
            "Steady incompressible laminar CFD with Boussinesq buoyancy and conjugate "
            "solid-air energy. Fans are represented as boundary-flow patches."
        )), 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        self.SetSizer(main)

        self.add_patch.Bind(wx.EVT_BUTTON, self._on_add)
        self.edit_patch.Bind(wx.EVT_BUTTON, self._on_edit)
        self.remove_patch.Bind(wx.EVT_BUTTON, self._on_remove)
        self.default_pair.Bind(wx.EVT_BUTTON, self._on_default_pair)
        self.patch_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self._on_edit)
        for control in (self.width, self.depth, self.height, self.cell):
            control.Bind(wx.EVT_TEXT, self._update_estimate)
        self._update_estimate()

    def _update_estimate(self, event=None):
        try:
            width, depth, height = float(self.width.GetValue()), float(self.depth.GetValue()), float(self.height.GetValue())
            cell = max(0.001, float(self.cell.GetValue()))
            count = (max(3, math.ceil(width / cell)) *
                     max(3, math.ceil(depth / cell)) *
                     max(3, math.ceil(height / cell)))
            self.estimate.SetLabel(f"Estimated cells: {count:,}")
        except ValueError:
            self.estimate.SetLabel("Estimated cells: invalid settings")
        if event:
            event.Skip()

    def _selected_index(self):
        return self.patch_list.GetFirstSelected()

    def _on_add(self, event):
        dialog = CFDBoundaryDialog(self)
        try:
            if dialog.ShowModal() == wx.ID_OK:
                self.settings.patches.append(dialog.get_patch())
                self._refresh_patches()
        except ValueError as exc:
            wx.MessageBox(str(exc), "Invalid CFD Patch", wx.OK | wx.ICON_ERROR)
        finally:
            dialog.Destroy()

    def _on_edit(self, event):
        index = self._selected_index()
        if index < 0:
            return
        dialog = CFDBoundaryDialog(self, self.settings.patches[index])
        try:
            if dialog.ShowModal() == wx.ID_OK:
                self.settings.patches[index] = dialog.get_patch()
                self._refresh_patches()
        except ValueError as exc:
            wx.MessageBox(str(exc), "Invalid CFD Patch", wx.OK | wx.ICON_ERROR)
        finally:
            dialog.Destroy()

    def _on_remove(self, event):
        index = self._selected_index()
        if index >= 0:
            self.settings.patches.pop(index)
            self._refresh_patches()

    def _on_default_pair(self, event):
        ambient = float(self.ambient.GetValue())
        self.settings.patches.extend([
            CFDBoundaryPatch("Fan Inlet", "FAN", "XMIN", 0.5, 0.5, 0.35, 0.35, 1.0, ambient),
            CFDBoundaryPatch("Outlet", "OUTLET", "XMAX", 0.5, 0.5, 0.35, 0.35, 0.0, ambient),
        ])
        self._refresh_patches()

    def _refresh_patches(self):
        self.patch_list.DeleteAllItems()
        for patch in self.settings.patches:
            row = self.patch_list.InsertItem(self.patch_list.GetItemCount(), patch.name)
            self.patch_list.SetItem(row, 1, patch.kind)
            self.patch_list.SetItem(row, 2, patch.face)
            self.patch_list.SetItem(row, 3, f"{patch.center_u:g}, {patch.center_v:g}")
            self.patch_list.SetItem(row, 4, f"{patch.size_u:g}, {patch.size_v:g}")
            self.patch_list.SetItem(row, 5, f"{patch.velocity_m_s:g}")
            self.patch_list.SetItem(row, 6, f"{patch.temperature_c:g}")

    def get_settings(self):
        result = self.settings
        result.ambient_c = float(self.ambient.GetValue())
        geometry = result.geometry
        geometry.width_mm, geometry.depth_mm, geometry.height_mm = (
            float(self.width.GetValue()), float(self.depth.GetValue()), float(self.height.GetValue())
        )
        geometry.board_orientation = self.orientation.GetStringSelection() or "XY"
        geometry.board_offset_x_mm = float(self.offset_x.GetValue())
        geometry.board_offset_y_mm = float(self.offset_y.GetValue())
        geometry.board_offset_z_mm = float(self.offset_z.GetValue())
        geometry.wall_heat_transfer_w_m2k = float(self.wall_h.GetValue())
        solver = result.solver
        solver.cell_size_mm = float(self.cell.GetValue())
        solver.max_iterations = int(self.iterations.GetValue())
        solver.tolerance = float(self.tolerance.GetValue())
        solver.relaxation = float(self.relaxation.GetValue())
        solver.pseudo_time_step_s = float(self.time_step.GetValue())
        solver.pressure_iterations = int(self.pressure_iterations.GetValue())
        solver.include_buoyancy = self.buoyancy.GetValue()
        result.use_phase3_heat_sources = self.phase3_heat.GetValue()
        result.include_dc_copper_losses = self.dc_losses.GetValue()
        if min(geometry.width_mm, geometry.depth_mm, geometry.height_mm, solver.cell_size_mm) <= 0:
            raise ValueError("Enclosure dimensions and CFD cell size must be greater than zero.")
        if solver.max_iterations < 1 or solver.pressure_iterations < 1:
            raise ValueError("CFD iteration limits must be greater than zero.")
        return result

    def set_settings(self, settings):
        self.settings = settings or EnclosureCFDSettings()
        geometry, solver = self.settings.geometry, self.settings.solver
        for control, value in (
            (self.width, geometry.width_mm), (self.depth, geometry.depth_mm),
            (self.height, geometry.height_mm), (self.offset_x, geometry.board_offset_x_mm),
            (self.offset_y, geometry.board_offset_y_mm), (self.offset_z, geometry.board_offset_z_mm),
            (self.wall_h, geometry.wall_heat_transfer_w_m2k), (self.ambient, self.settings.ambient_c),
            (self.cell, solver.cell_size_mm), (self.iterations, solver.max_iterations),
            (self.tolerance, solver.tolerance), (self.relaxation, solver.relaxation),
            (self.time_step, solver.pseudo_time_step_s),
            (self.pressure_iterations, solver.pressure_iterations),
        ):
            control.SetValue(f"{value:g}")
        index = self.orientation.FindString(geometry.board_orientation)
        self.orientation.SetSelection(index if index != wx.NOT_FOUND else 0)
        self.buoyancy.SetValue(solver.include_buoyancy)
        self.phase3_heat.SetValue(self.settings.use_phase3_heat_sources)
        self.dc_losses.SetValue(self.settings.include_dc_copper_losses)
        self._refresh_patches()
        self._update_estimate()
