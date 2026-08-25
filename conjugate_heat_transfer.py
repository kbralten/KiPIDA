"""Phase 4 orchestration for PCB solid conduction and enclosure airflow."""

try:
    from .cfd_mesh import CFDMeshGenerator
    from .cfd_model import EnclosureModelBuilder
    from .cfd_solver import EnclosureCFDSolver
except (ImportError, ValueError):
    from cfd_mesh import CFDMeshGenerator
    from cfd_model import EnclosureModelBuilder
    from cfd_solver import EnclosureCFDSolver


class ConjugateHeatTransferSolver:
    """Build and solve the enclosure using Phase 3 board heat-source data."""

    def __init__(self, debug=False, log_callback=None):
        self.debug = debug
        self.log_callback = log_callback

    def solve(
        self, board_model, settings, progress_callback=None, cancel_callback=None
    ):
        enclosure = EnclosureModelBuilder(
            debug=self.debug, log_callback=self.log_callback
        ).build(board_model, settings)
        mesh = CFDMeshGenerator(
            debug=self.debug, log_callback=self.log_callback
        ).generate_mesh(enclosure, settings)
        result = EnclosureCFDSolver(
            debug=self.debug, log_callback=self.log_callback
        ).solve(
            mesh,
            settings,
            progress_callback=progress_callback,
            cancel_callback=cancel_callback,
        )
        return mesh, result
