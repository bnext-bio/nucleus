import logging
import rich.console
import rich.logging

from logging import ERROR, WARNING, INFO, DEBUG


def set_module_logging_level(level):
    # Get the name of the current package
    package_name = __name__.split(".")[0]

    # Iterate over all loggers
    for name, logger in logging.root.manager.loggerDict.items():
        # Check if the logger belongs to our package
        if name == package_name or name.startswith(f"{package_name}."):
            if isinstance(
                logger, logging.Logger
            ):  # Ensure it's actually a logger
                logger.setLevel(level)


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
        level=INFO,
        format="%(message)s",
        handlers=[
            rich.logging.RichHandler(console=console, enable_link_path=False)
        ],
    )

    log = logging.getLogger(__name__)

    set_module_logging_level(log_level)

    log.info("Logging initialized")
    return log
