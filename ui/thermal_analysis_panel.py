import wx

try:
    from models import AirflowSettings, ThermalAnalysisSettings
    from thermal_model import PowerLossEstimator
except (ImportError, ValueError):
    from models import AirflowSettings, ThermalAnalysisSettings
    from thermal_model import PowerLossEstimator


class ThermalComponentDialog(wx.Dialog):
    def __init__(self, parent, component):
        super().__init__(parent, title=f"Thermal Model: {component.ref_des}")
        sizer = wx.BoxSizer(wx.VERTICAL)
        grid = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        grid.AddGrowableCol(1, 1)
        self.txt_power = wx.TextCtrl(self, value=f"{component.power_w:g}")
        self.txt_width = wx.TextCtrl(self, value=f"{component.width_mm:g}")
        self.txt_depth = wx.TextCtrl(self, value=f"{component.depth_mm:g}")
        self.txt_height = wx.TextCtrl(self, value=f"{component.height_mm:g}")
        self.txt_theta = wx.TextCtrl(self, value=f"{component.theta_jb_c_per_w:g}")
        self.txt_tjmax = wx.TextCtrl(self, value=f"{component.max_junction_c:g}")
        self.chk_enabled = wx.CheckBox(self, label="Include as heat source")
        self.chk_enabled.SetValue(component.enabled)
        for label, control in (
            ("Power (W):", self.txt_power),
            ("Width (mm):", self.txt_width),
            ("Depth (mm):", self.txt_depth),
            ("Height (mm):", self.txt_height),
            ("Theta JB (C/W):", self.txt_theta),
            ("Max junction (C):", self.txt_tjmax),
        ):
            grid.Add(wx.StaticText(self, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(control, 1, wx.EXPAND)
        grid.Add(wx.StaticText(self, label="Enabled:"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.chk_enabled, 0)
        sizer.Add(grid, 1, wx.EXPAND | wx.ALL, 12)
        sizer.Add(self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL), 0, wx.EXPAND | wx.ALL, 12)
        self.SetSizerAndFit(sizer)

    def apply_to(self, component):
        values = [
            float(self.txt_power.GetValue()),
            float(self.txt_width.GetValue()),
            float(self.txt_depth.GetValue()),
            float(self.txt_height.GetValue()),
            float(self.txt_theta.GetValue()),
            float(self.txt_tjmax.GetValue()),
        ]
        if any(value < 0 for value in values[:-1]) or values[-1] <= -273.15:
            raise ValueError("Thermal model values must be non-negative and physically valid.")
        component.power_w, component.width_mm, component.depth_mm = values[:3]
        component.height_mm, component.theta_jb_c_per_w, component.max_junction_c = values[3:]
        component.enabled = self.chk_enabled.GetValue()
        component.model_source = "user"


class ThermalAnalysisPanel(wx.Panel):
    def __init__(self, parent, rails_provider, log_callback=None):
        super().__init__(parent)
        self.rails_provider = rails_provider
        self.log_callback = log_callback
        self.settings = ThermalAnalysisSettings()
        self._init_ui()

    def _init_ui(self):
        main = wx.BoxSizer(wx.VERTICAL)
        settings_box = wx.StaticBoxSizer(wx.VERTICAL, self, "3D Thermal and Airflow Setup")
        settings_parent = settings_box.GetStaticBox()
        grid = wx.FlexGridSizer(cols=4, hgap=6, vgap=6)
        grid.AddGrowableCol(1, 1)
        grid.AddGrowableCol(3, 1)

        self.txt_ambient = wx.TextCtrl(settings_parent, value="25")
        self.txt_grid = wx.TextCtrl(settings_parent, value="1.0")
        self.choice_airflow = wx.Choice(settings_parent, choices=["NATURAL", "FORCED", "CUSTOM"])
        self.choice_airflow.SetSelection(0)
        self.txt_velocity = wx.TextCtrl(settings_parent, value="0")
        self.txt_direction = wx.TextCtrl(settings_parent, value="0")
        self.txt_custom_h = wx.TextCtrl(settings_parent, value="10")
        self.txt_iterations = wx.TextCtrl(settings_parent, value="6")
        self.txt_convergence = wx.TextCtrl(settings_parent, value="0.1")

        rows = [
            ("Ambient (C):", self.txt_ambient, "Grid size (mm):", self.txt_grid),
            ("Airflow mode:", self.choice_airflow, "Air speed (m/s):", self.txt_velocity),
            ("Air direction (deg):", self.txt_direction, "Custom h (W/m2K):", self.txt_custom_h),
            ("Coupled iterations:", self.txt_iterations, "Convergence (C):", self.txt_convergence),
        ]
        for left_label, left_control, right_label, right_control in rows:
            grid.Add(wx.StaticText(settings_parent, label=left_label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(left_control, 1, wx.EXPAND)
            grid.Add(wx.StaticText(settings_parent, label=right_label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(right_control, 1, wx.EXPAND)
        settings_box.Add(grid, 0, wx.EXPAND | wx.ALL, 8)

        checks = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_top = wx.CheckBox(settings_parent, label="Top exposed")
        self.chk_bottom = wx.CheckBox(settings_parent, label="Bottom exposed")
        self.chk_edges = wx.CheckBox(settings_parent, label="Edges exposed")
        self.chk_radiation = wx.CheckBox(settings_parent, label="Include radiation")
        self.chk_dc_loss = wx.CheckBox(settings_parent, label="Include DC copper losses")
        for check in (self.chk_top, self.chk_bottom, self.chk_edges, self.chk_radiation, self.chk_dc_loss):
            check.SetValue(True)
            checks.Add(check, 0, wx.RIGHT, 14)
        settings_box.Add(checks, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 8)
        main.Add(settings_box, 0, wx.EXPAND | wx.ALL, 5)

        component_box = wx.StaticBoxSizer(wx.VERTICAL, self, "Component Heat Sources")
        component_parent = component_box.GetStaticBox()
        self.component_list = wx.ListCtrl(component_parent, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        for index, (label, width) in enumerate((
            ("Use", 50), ("Ref", 70), ("Power (W)", 90), ("Size (mm)", 105),
            ("Theta JB", 80), ("Tj max", 75), ("Source", 150),
        )):
            self.component_list.InsertColumn(index, label, width=width)
        component_box.Add(self.component_list, 1, wx.EXPAND | wx.ALL, 5)
        buttons = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_estimate = wx.Button(component_parent, label="Refresh Power Estimates")
        self.btn_toggle = wx.Button(component_parent, label="Enable / Disable Selected")
        self.btn_edit = wx.Button(component_parent, label="Edit Selected Model")
        buttons.Add(self.btn_estimate, 0, wx.RIGHT, 5)
        buttons.Add(self.btn_toggle, 0, wx.RIGHT, 5)
        buttons.Add(self.btn_edit, 0)
        component_box.Add(buttons, 0, wx.ALL, 5)
        main.Add(component_box, 1, wx.EXPAND | wx.ALL, 5)
        main.Add(wx.StaticText(
            self,
            label=("Steady-state 3D solid conduction with convective boundaries. "
                   "Estimated powers and package models must be reviewed before sign-off."),
        ), 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        self.SetSizer(main)

        self.btn_estimate.Bind(wx.EVT_BUTTON, self._on_estimate)
        self.btn_toggle.Bind(wx.EVT_BUTTON, self._on_toggle)
        self.btn_edit.Bind(wx.EVT_BUTTON, self._on_edit)
        self.component_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self._on_edit)

    def refresh_components(self, preserve_user=True):
        saved = {component.ref_des: component for component in self.settings.components}
        estimated = PowerLossEstimator.estimate(list(self.rails_provider() or []))
        if preserve_user:
            for index, component in enumerate(estimated):
                prior = saved.get(component.ref_des)
                if prior is not None and prior.model_source == "user":
                    estimated[index] = prior
        self.settings.components = estimated
        self._update_component_list()

    def _update_component_list(self):
        self.component_list.DeleteAllItems()
        for component in self.settings.components:
            row = self.component_list.InsertItem(
                self.component_list.GetItemCount(), "Yes" if component.enabled else "No"
            )
            self.component_list.SetItem(row, 1, component.ref_des)
            self.component_list.SetItem(row, 2, f"{component.power_w:g}")
            self.component_list.SetItem(row, 3, f"{component.width_mm:g} x {component.depth_mm:g}")
            self.component_list.SetItem(row, 4, f"{component.theta_jb_c_per_w:g}")
            self.component_list.SetItem(row, 5, f"{component.max_junction_c:g}")
            self.component_list.SetItem(row, 6, component.model_source)

    def _selected_component(self):
        selected = self.component_list.GetFirstSelected()
        if selected == -1 or selected >= len(self.settings.components):
            return None
        return self.settings.components[selected]

    def _on_estimate(self, event):
        self.refresh_components(preserve_user=True)

    def _on_toggle(self, event):
        component = self._selected_component()
        if component is not None:
            component.enabled = not component.enabled
            self._update_component_list()

    def _on_edit(self, event):
        component = self._selected_component()
        if component is None:
            return
        dialog = ThermalComponentDialog(self, component)
        try:
            if dialog.ShowModal() == wx.ID_OK:
                dialog.apply_to(component)
                self._update_component_list()
        except ValueError as exc:
            wx.MessageBox(str(exc), "Invalid Thermal Model", wx.OK | wx.ICON_ERROR)
        finally:
            dialog.Destroy()

    def get_settings(self):
        mode = self.choice_airflow.GetStringSelection() or "NATURAL"
        self.settings.ambient_c = float(self.txt_ambient.GetValue())
        self.settings.grid_size_mm = float(self.txt_grid.GetValue())
        self.settings.airflow = AirflowSettings(
            mode=mode,
            velocity_m_s=float(self.txt_velocity.GetValue()),
            direction_deg=float(self.txt_direction.GetValue()),
            custom_h_w_m2k=float(self.txt_custom_h.GetValue()),
            expose_top=self.chk_top.GetValue(),
            expose_bottom=self.chk_bottom.GetValue(),
            expose_edges=self.chk_edges.GetValue(),
        )
        self.settings.include_radiation = self.chk_radiation.GetValue()
        self.settings.include_dc_copper_losses = self.chk_dc_loss.GetValue()
        self.settings.coupled_iterations = int(self.txt_iterations.GetValue())
        self.settings.convergence_c = float(self.txt_convergence.GetValue())
        if self.settings.grid_size_mm <= 0 or self.settings.coupled_iterations < 1:
            raise ValueError("Grid size and coupled iterations must be greater than zero.")
        return self.settings

    def set_settings(self, settings):
        self.settings = settings or ThermalAnalysisSettings()
        airflow = self.settings.airflow
        self.txt_ambient.SetValue(f"{self.settings.ambient_c:g}")
        self.txt_grid.SetValue(f"{self.settings.grid_size_mm:g}")
        index = self.choice_airflow.FindString(airflow.mode)
        self.choice_airflow.SetSelection(index if index != wx.NOT_FOUND else 0)
        self.txt_velocity.SetValue(f"{airflow.velocity_m_s:g}")
        self.txt_direction.SetValue(f"{airflow.direction_deg:g}")
        self.txt_custom_h.SetValue(f"{airflow.custom_h_w_m2k:g}")
        self.txt_iterations.SetValue(str(self.settings.coupled_iterations))
        self.txt_convergence.SetValue(f"{self.settings.convergence_c:g}")
        self.chk_top.SetValue(airflow.expose_top)
        self.chk_bottom.SetValue(airflow.expose_bottom)
        self.chk_edges.SetValue(airflow.expose_edges)
        self.chk_radiation.SetValue(self.settings.include_radiation)
        self.chk_dc_loss.SetValue(self.settings.include_dc_copper_losses)
        # Recompute all automatic entries when loading a project so configs
        # saved by older versions cannot retain V*I connector heat or output-
        # inductor regulator placement. Explicit user models remain untouched.
        self.refresh_components(preserve_user=True)
