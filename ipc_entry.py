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
    log_file = os.path.join(tempfile.gettempdir(), 'kipida_entry.log')
    logging.basicConfig(filename=log_file, level=logging.DEBUG, filemode='w')
    logger = logging.getLogger("KiPIDA")
    logger.info("Ki-PIDA Entry Point Started")
    
    # 1. Connect to KiCad
    # Use our robust connection logic (same as in extractor.py)
    socket_path = os.environ.get("KICAD_API_SOCKET")
    if not socket_path:
        temp_dir = tempfile.gettempdir()
        potential_path = os.path.join(temp_dir, 'kicad', 'api.sock')
        
        # On Windows, skip existence check for named pipes
        if sys.platform == 'win32':
             socket_path = potential_path
        elif os.path.exists(potential_path):
             socket_path = potential_path
             
    # Ensure protocol prefix
    if socket_path and not socket_path.startswith('ipc://') and not socket_path.startswith('\\\\.\\pipe\\'):
         socket_path = f"ipc://{socket_path}"

    initial_logs = []
    def early_log(msg):
        logger.info(msg)
        initial_logs.append(msg)

    # early_log(f"Connecting to KiCad at {socket_path}")
    
    client = None
    board = None
    
    try:
        # early_log("Initializing kipy client...")
        client = kipy.KiCad(socket_path=socket_path, timeout_ms=3000)
        # early_log("Requesting board object...")
        board = client.get_board()
        if board:
            early_log("Ki-PIDA connected successfully.")
            # Quick sanity check
            if hasattr(board, 'footprints'):
                 pass
            else:
                early_log("WARNING: Board object retrieved but missing 'footprints' attribute.")
        else:
            early_log("ERROR: client.get_board() returned None.")
    except Exception as e:
        early_log(f"Connection failed with error: {e}")
        import traceback
        early_log(traceback.format_exc())
        
    # 2. Initialize wx App
    app = wx.App()
    
    # 3. Create and Show Dialog
    try:
        dlg = KiPIDA_MainDialog(None, board_adapter=board) 
    except Exception as e:
        early_log(f"CRITICAL: Failed to create KiPIDA_MainDialog: {e}")
        import traceback
        early_log(traceback.format_exc())
        # Show a simple message box if possible
        dlg_err = wx.MessageDialog(None, f"Failed to start Ki-PIDA: {e}\n\nCheck logs for details.", "Startup Error", wx.OK | wx.ICON_ERROR)
        dlg_err.ShowModal()
        dlg_err.Destroy()
        return
    
    # Pass initial logs to UI
    for msg in initial_logs:
        dlg.log(f"[INIT] {msg}")
    
    try:
        dlg.ShowModal()
    except Exception as e:
        logger.error(f"Runtime error: {e}")
        import traceback
        traceback.print_exc()
    finally:
        dlg.Destroy()
        
    app.MainLoop()

if __name__ == "__main__":
    main()
