import os
import sys
from pathlib import Path
from venv import logger

if sys.platform == "darwin":
    # Prefer the managed virtual environment's installed packages
    # over paths inherited from KiCad's bundled Python.
    venv_site = (
        Path(sys.prefix)
        / "lib"
        / f"python{sys.version_info.major}.{sys.version_info.minor}"
        / "site-packages"
    )

    sys.path.insert(0, str(venv_site))


import wx
import sys
import os
import logging
import tempfile
import kipy
import kipy.board

# Add current directory to path so we can import local modules
plugin_dir = os.path.dirname(os.path.abspath(__file__))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

from ui.main_dialog import KiPIDA_MainDialog


def main():
    # Setup logging to file for debugging
    log_file = os.path.join(tempfile.gettempdir(), "kipida_entry.log")
    logging.basicConfig(filename=log_file, level=logging.DEBUG, filemode="w")
    logger = logging.getLogger("KiPIDA")
    logger.info("Ki-PIDA Entry Point Started")

    # 1. Connect to KiCad
    # Use our robust connection logic (same as in extractor.py)
    socket_path = os.environ.get("KICAD_API_SOCKET")
    if not socket_path:
        temp_dir = tempfile.gettempdir()
        potential_path = os.path.join(temp_dir, "kicad", "api.sock")

        # On Windows, skip existence check for named pipes
        if sys.platform == "win32":
            socket_path = potential_path
        elif os.path.exists(potential_path):
            socket_path = potential_path

    # Ensure protocol prefix
    if (
        socket_path
        and not socket_path.startswith("ipc://")
        and not socket_path.startswith("\\\\.\\pipe\\")
    ):
        socket_path = f"ipc://{socket_path}"

    initial_logs = []

    def early_log(msg):
        logger.info(msg)
        initial_logs.append(msg)

    # early_log(f"Connecting to KiCad at {socket_path}")

    client = None
    board = None
    project = None

    try:
        # early_log("Initializing kipy client...")
        client = kipy.KiCad(socket_path=socket_path, timeout_ms=3000)
        # early_log("Requesting board object...")
        board = client.get_board()

        # Try to get project info from board.document without calling get_project
        # to avoid breaking board API access
        project = None
        if board and hasattr(board, "document"):
            try:
                doc = board.document
                # Document is a protobuf DocumentSpecifier with a nested 'project' field

                # Create a simple object with path and name attributes
                class ProjectInfo:
                    def __init__(self, doc):
                        self.path = None
                        self.name = None

                        # Extract the 'project' field from DocumentSpecifier
                        if hasattr(doc, "ListFields"):
                            fields = doc.ListFields()

                            for field_desc, value in fields:
                                field_name = field_desc.name

                                # The 'project' field contains the project info
                                if field_name == "project":
                                    # value is a protobuf with 'name' and 'path' attributes
                                    if hasattr(value, "path"):
                                        self.path = value.path
                                    if hasattr(value, "name"):
                                        self.name = value.name
                                    break

                        # If we have path but no name, extract from path
                        if self.path and not self.name:
                            import os

                            self.name = os.path.basename(self.path)

                project = ProjectInfo(doc)
                if project.path:
                    early_log(f"Project: {project.name} at {project.path}")
                else:
                    early_log("Could not find project path in document")
                    project = None
            except Exception as e:
                early_log(f"Failed to get project info from board.document: {e}")
                import traceback

                early_log(traceback.format_exc())
                project = None
        else:
            early_log("Board has no document attribute")

        if board:
            early_log("Ki-PIDA connected successfully.")
        else:
            early_log("ERROR: client.get_board() returned None.")
    except Exception as e:
        early_log(f"Connection failed with error: {e}")
        import traceback

        early_log(traceback.format_exc())

    # 2. Initialize wx App
    app = wx.App(False)

    # 3. Create and Show Dialog
    try:
        dlg = KiPIDA_MainDialog(None, board_adapter=board, project=project)
    except Exception as e:
        early_log(f"CRITICAL: Failed to create KiPIDA_MainDialog: {e}")
        import traceback

        early_log(traceback.format_exc())
        # Show a simple message box if possible
        dlg_err = wx.MessageDialog(
            None,
            f"Failed to start Ki-PIDA: {e}\n\nCheck logs for details.",
            "Startup Error",
            wx.OK | wx.ICON_ERROR,
        )
        dlg_err.ShowModal()
        dlg_err.Destroy()
        return

    # Pass initial logs to UI
    for msg in initial_logs:
        dlg.log(f"[INIT] {msg}")

    try:
        if sys.platform == "darwin":

            def activate_app():
                try:
                    from AppKit import NSApplication

                    ns_app = NSApplication.sharedApplication()
                    ns_app.activateIgnoringOtherApps_(True)
                    dlg.Raise()

                    logger.info(
                        "Delayed activation: active %s",
                        ns_app.isActive(),
                    )

                except Exception:
                    logger.exception("macOS application activation failed")

            wx.CallLater(100, activate_app)

        dlg.ShowModal()

    except Exception:
        logger.exception("KiPIDA dialog runtime error")

    finally:
        dlg.Destroy()


if __name__ == "__main__":
    main()
