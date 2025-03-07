import logging
import rich.console
import rich.logging


def setup_logging(log_level=logging.INFO):
    # Set up logging
    console = rich.console.Console(
        force_jupyter=False,
        # stderr=True,
        theme=rich.theme.Theme(
            {"logging.level.debug": "cyan", "logging.level.info": "green"}
        ),
    )

    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        handlers=[rich.logging.RichHandler(
            console=console,
            enable_link_path=False)
        ],
    )
    
    log = logging.getLogger(__name__)
    log.info("Logging initialized")
    return log
