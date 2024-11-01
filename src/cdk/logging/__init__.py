import logging
import rich.console
import rich.logging


def setup_logging():
    # Set up logging
    console = rich.console.Console(
        force_jupyter=False,
        # stderr=True,
        theme=rich.theme.Theme(
            {"logging.level.debug": "cyan", "logging.level.info": "green"}
        ),
    )

    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    log.handlers.clear()
    log.addHandler(
        rich.logging.RichHandler(console=console, enable_link_path=False)
    )
    log.info("Logging initialized")
    return log
