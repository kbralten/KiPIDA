import wx

try:
    from ac_model import ACModelBuilder, format_capacitance, parse_capacitance
    from models import ACAnalysisSettings, ACMeasurementPort, ACSourceModel
except (ImportError, ValueError):
    from ac_model import ACModelBuilder, format_capacitance, parse_capacitance
    from models import ACAnalysisSettings, ACMeasurementPort, ACSourceModel


class CapacitorModelDialog(wx.Dialog):
    def __init__(self, parent, capacitor):
        super().__init__(parent, title=f"Capacitor Model: {capacitor.ref_des}")
        sizer = wx.BoxSizer(wx.VERTICAL)
        grid = wx.FlexGridSizer(cols=2, hgap=8, vgap=8)
        grid.AddGrowableCol(1, 1)
        self.txt_capacitance = wx.TextCtrl(self, value=format_capacitance(capacitor.capacitance_f))
        self.txt_esr = wx.TextCtrl(self, value=f"{capacitor.esr_ohm * 1e3:g}")
        self.txt_esl = wx.TextCtrl(self, value=f"{capacitor.esl_h * 1e9:g}")
        for label, control in (
            ("Capacitance:", self.txt_capacitance),
            ("ESR (mOhm):", self.txt_esr),
            ("ESL (nH):", self.txt_esl),
        ):
            grid.Add(wx.StaticText(self, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(control, 1, wx.EXPAND)
        sizer.Add(grid, 1, wx.EXPAND | wx.ALL, 12)

        buttons = self.CreateStdDialogButtonSizer(wx.OK | wx.CANCEL)
        sizer.Add(buttons, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 12)
        self.SetSizerAndFit(sizer)

    def apply_to(self, capacitor):
        capacitance = parse_capacitance(self.txt_capacitance.GetValue())
        if capacitance is None:
            raise ValueError("Enter a valid capacitance such as 100n, 4u7, or 10uF.")
        esr = float(self.txt_esr.GetValue()) / 1e3
        esl = float(self.txt_esl.GetValue()) * 1e-9
        if esr < 0 or esl < 0:
            raise ValueError("ESR and ESL must be zero or greater.")
        capacitor.capacitance_f = capacitance
        capacitor.esr_ohm = esr
        capacitor.esl_h = esl
        capacitor.model_source = "user"


class ACAnalysisPanel(wx.Panel):
    """Configuration surface for impedance sweeps and decoupling optimization."""

    def __init__(self, parent, board, rails_provider, log_callback=None):
        super().__init__(parent)
        self.board = board
        self.rails_provider = rails_provider
        self.log_callback = log_callback
        self.profiles = {}
        self.capacitors = []
        self._updating = False
        self._profile_rail_name = ""
        self.builder = ACModelBuilder(board, log_callback=log_callback)
        self._init_ui()

    def _init_ui(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)

        settings_box = wx.StaticBoxSizer(wx.VERTICAL, self, "AC Impedance Setup")
        settings_parent = settings_box.GetStaticBox()
        grid = wx.FlexGridSizer(cols=4, hgap=6, vgap=6)
        grid.AddGrowableCol(1, 1)
        grid.AddGrowableCol(3, 1)

        self.choice_rail = wx.Choice(settings_parent)
        self.choice_ground = wx.Choice(settings_parent)
        self.choice_source = wx.Choice(settings_parent)
        self.choice_port = wx.Choice(settings_parent)
        self.txt_f_start = wx.TextCtrl(settings_parent, value="1e3")
        self.txt_f_stop = wx.TextCtrl(settings_parent, value="1e8")
        self.txt_points = wx.TextCtrl(settings_parent, value="121")
        self.txt_target_mohm = wx.TextCtrl(settings_parent, value="50")
        self.txt_source_r_mohm = wx.TextCtrl(settings_parent, value="10")
        self.txt_source_l_nh = wx.TextCtrl(settings_parent, value="1")
        self.txt_max_additions = wx.TextCtrl(settings_parent, value="8")

        rows = [
            ("Power rail:", self.choice_rail, "Return net:", self.choice_ground),
            ("Source component:", self.choice_source, "Measurement component:", self.choice_port),
            ("Start frequency (Hz):", self.txt_f_start, "Stop frequency (Hz):", self.txt_f_stop),
            ("Frequency points:", self.txt_points, "Target |Z| (mOhm):", self.txt_target_mohm),
            ("Source R (mOhm):", self.txt_source_r_mohm, "Source L (nH):", self.txt_source_l_nh),
            ("Max capacitor additions:", self.txt_max_additions, "", wx.StaticText(settings_parent, label="")),
        ]
        for left_label, left_ctrl, right_label, right_ctrl in rows:
            grid.Add(wx.StaticText(settings_parent, label=left_label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(left_ctrl, 1, wx.EXPAND)
            grid.Add(wx.StaticText(settings_parent, label=right_label), 0, wx.ALIGN_CENTER_VERTICAL)
            grid.Add(right_ctrl, 1, wx.EXPAND)
        settings_box.Add(grid, 1, wx.EXPAND | wx.ALL, 8)
        main_sizer.Add(settings_box, 0, wx.EXPAND | wx.ALL, 5)

        cap_box = wx.StaticBoxSizer(wx.VERTICAL, self, "Detected Rail-to-Ground Capacitors")
        cap_parent = cap_box.GetStaticBox()
        self.cap_list = wx.ListCtrl(cap_parent, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        for index, (title, width) in enumerate((
            ("Use", 55), ("Ref", 70), ("Capacitance", 100), ("ESR (mOhm)", 95),
            ("ESL (nH)", 85), ("Status", 100), ("Model", 190),
        )):
            self.cap_list.InsertColumn(index, title, width=width)
        cap_box.Add(self.cap_list, 1, wx.EXPAND | wx.ALL, 5)

        cap_buttons = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_refresh = wx.Button(cap_parent, label="Refresh Detection")
        self.btn_toggle = wx.Button(cap_parent, label="Enable / Disable Selected")
        self.btn_edit = wx.Button(cap_parent, label="Edit Selected Model")
        cap_buttons.Add(self.btn_refresh, 0, wx.RIGHT, 5)
        cap_buttons.Add(self.btn_toggle, 0, wx.RIGHT, 5)
        cap_buttons.Add(self.btn_edit, 0)
        cap_box.Add(cap_buttons, 0, wx.ALL, 5)
        main_sizer.Add(cap_box, 1, wx.EXPAND | wx.ALL, 5)

        note = wx.StaticText(
            self,
            label="ESR/ESL values inferred from package are estimates. Review them before sign-off.",
        )
        main_sizer.Add(note, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, 10)
        self.SetSizer(main_sizer)

        self.choice_rail.Bind(wx.EVT_CHOICE, self._on_context_changed)
        self.choice_ground.Bind(wx.EVT_CHOICE, self._on_context_changed)
        self.btn_refresh.Bind(wx.EVT_BUTTON, self._on_refresh)
        self.btn_toggle.Bind(wx.EVT_BUTTON, self._on_toggle)
        self.btn_edit.Bind(wx.EVT_BUTTON, self._on_edit)
        self.cap_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self._on_edit)

    def _rails(self):
        return list(self.rails_provider() or [])

    def _selected_text(self, choice):
        selection = choice.GetSelection()
        return choice.GetString(selection) if selection != wx.NOT_FOUND else ""

    def _select_text(self, choice, value):
        if not value:
            if choice.GetCount():
                choice.SetSelection(0)
            return
        index = choice.FindString(value)
        if index != wx.NOT_FOUND:
            choice.SetSelection(index)
        elif choice.GetCount():
            choice.SetSelection(0)

    def _active_rail(self, name=None):
        name = name or self._selected_text(self.choice_rail)
        return next((rail for rail in self._rails() if rail.net_name == name), None)

    def _source_options(self, rail):
        """Return direct sources plus regulator outputs feeding this rail."""
        options = {
            source.component_ref.ref_des: list(source.pad_names)
            for source in rail.sources
        }
        for parent_rail in self._rails():
            for regulator in parent_rail.child_regulators:
                if regulator.output_rail_name == rail.net_name and regulator.output_ref_des:
                    options.setdefault(regulator.output_ref_des, list(regulator.output_pad_names))
        return options

    def refresh(self, force_discovery=False):
        if self._updating:
            return
        self._updating = True
        try:
            previous_rail = self._selected_text(self.choice_rail)
            rail_names = [rail.net_name for rail in self._rails()]
            self.choice_rail.Set(rail_names)
            self._select_text(self.choice_rail, previous_rail)

            ground_names = self.builder.discover_ground_nets()
            previous_ground = self._selected_text(self.choice_ground)
            self.choice_ground.Set(ground_names)
            self._select_text(self.choice_ground, previous_ground or "GND")

            rail = self._active_rail()
            if rail is None:
                self.capacitors = []
                self._update_cap_list()
                return

            profile = self.profiles.get(rail.net_name)
            if profile is None:
                profile = ACAnalysisSettings(rail_name=rail.net_name)
                self.profiles[rail.net_name] = profile

            self._select_text(self.choice_ground, profile.ground_net_name or "GND")
            source_refs = sorted(self._source_options(rail))
            load_refs = sorted({load.component_ref.ref_des for load in rail.loads})
            self.choice_source.Set(source_refs)
            self.choice_port.Set(load_refs)
            self._select_text(self.choice_source, profile.source.ref_des)
            self._select_text(self.choice_port, profile.measurement_port.ref_des)

            self.txt_f_start.SetValue(f"{profile.frequency_start_hz:g}")
            self.txt_f_stop.SetValue(f"{profile.frequency_stop_hz:g}")
            self.txt_points.SetValue(str(profile.frequency_points))
            self.txt_target_mohm.SetValue(f"{profile.target_impedance_ohm * 1e3:g}")
            self.txt_source_r_mohm.SetValue(f"{profile.source.resistance_ohm * 1e3:g}")
            self.txt_source_l_nh.SetValue(f"{profile.source.inductance_h * 1e9:g}")
            self.txt_max_additions.SetValue(str(profile.optimizer_max_additions))

            ground_name = self._selected_text(self.choice_ground)
            if force_discovery or not profile.capacitors:
                discovered = self.builder.discover_capacitors(rail.net_name, ground_name)
                existing = {cap.ref_des: cap for cap in profile.capacitors}
                for capacitor in discovered:
                    if capacitor.ref_des in existing:
                        saved = existing[capacitor.ref_des]
                        capacitor.enabled = saved.enabled
                        capacitor.capacitance_f = saved.capacitance_f
                        capacitor.esr_ohm = saved.esr_ohm
                        capacitor.esl_h = saved.esl_h
                        capacitor.model_source = saved.model_source
                profile.capacitors = discovered
            self.capacitors = profile.capacitors
            self._profile_rail_name = rail.net_name
            self._update_cap_list()
        finally:
            self._updating = False

    def _update_cap_list(self):
        self.cap_list.DeleteAllItems()
        for capacitor in self.capacitors:
            row = self.cap_list.InsertItem(self.cap_list.GetItemCount(), "Yes" if capacitor.enabled else "No")
            self.cap_list.SetItem(row, 1, capacitor.ref_des)
            self.cap_list.SetItem(row, 2, format_capacitance(capacitor.capacitance_f))
            self.cap_list.SetItem(row, 3, f"{capacitor.esr_ohm * 1e3:g}")
            self.cap_list.SetItem(row, 4, f"{capacitor.esl_h * 1e9:g}")
            self.cap_list.SetItem(row, 5, "Candidate/DNP" if capacitor.candidate else "Populated")
            self.cap_list.SetItem(row, 6, capacitor.model_source)

    def _on_context_changed(self, event):
        if self._updating:
            return
        self._save_current_profile()
        self.refresh(force_discovery=True)

    def _on_refresh(self, event):
        self._save_current_profile()
        self.refresh(force_discovery=True)

    def _on_toggle(self, event):
        selected = self.cap_list.GetFirstSelected()
        if selected == -1 or selected >= len(self.capacitors):
            return
        self.capacitors[selected].enabled = not self.capacitors[selected].enabled
        self._update_cap_list()

    def _on_edit(self, event):
        selected = self.cap_list.GetFirstSelected()
        if selected == -1 or selected >= len(self.capacitors):
            return
        dialog = CapacitorModelDialog(self, self.capacitors[selected])
        try:
            if dialog.ShowModal() != wx.ID_OK:
                return
            dialog.apply_to(self.capacitors[selected])
            self._update_cap_list()
        except ValueError as exc:
            wx.MessageBox(str(exc), "Invalid Capacitor Model", wx.OK | wx.ICON_ERROR)
        finally:
            dialog.Destroy()

    def _save_current_profile(self):
        rail = self._active_rail(self._profile_rail_name)
        if rail is None:
            rail = self._active_rail()
        if rail is None:
            return None
        profile = self.profiles.get(rail.net_name, ACAnalysisSettings(rail_name=rail.net_name))
        profile.rail_name = rail.net_name
        profile.ground_net_name = self._selected_text(self.choice_ground) or "GND"
        profile.frequency_start_hz = float(self.txt_f_start.GetValue())
        profile.frequency_stop_hz = float(self.txt_f_stop.GetValue())
        profile.frequency_points = int(self.txt_points.GetValue())
        profile.target_impedance_ohm = float(self.txt_target_mohm.GetValue()) / 1e3
        profile.optimizer_max_additions = int(self.txt_max_additions.GetValue())

        source_ref = self._selected_text(self.choice_source)
        source_pad_names = self._source_options(rail).get(source_ref, [])
        profile.source = ACSourceModel(
            ref_des=source_ref,
            rail_pad_names=source_pad_names,
            ground_pad_names=self.builder.pad_names_for_net(source_ref, profile.ground_net_name),
            resistance_ohm=float(self.txt_source_r_mohm.GetValue()) / 1e3,
            inductance_h=float(self.txt_source_l_nh.GetValue()) * 1e-9,
        )

        port_ref = self._selected_text(self.choice_port)
        load = next((item for item in rail.loads if item.component_ref.ref_des == port_ref), None)
        profile.measurement_port = ACMeasurementPort(
            ref_des=port_ref,
            rail_pad_names=list(load.pad_names) if load else [],
            ground_pad_names=self.builder.pad_names_for_net(port_ref, profile.ground_net_name),
        )
        profile.capacitors = self.capacitors
        self.profiles[rail.net_name] = profile
        return profile

    def get_settings(self):
        return self._save_current_profile()

    def get_profiles(self):
        self._save_current_profile()
        return self.profiles

    def set_profiles(self, profiles):
        self.profiles = dict(profiles or {})
        self.refresh(force_discovery=False)
