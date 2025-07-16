"""
Opentrons Utilities

This module provides utility functions for use in Opentrons protocols.

Usage:
*tbd*
"""

import logging

from opentrons import protocol_api

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def get_logger(protocol: protocol_api.ProtocolContext):
    def log_comment(msg, **kwargs):
        structured_log = " ".join([f"{k}={v}" for k, v in kwargs.items()])

        log.info(f"msg={msg} {structured_log}")
        protocol.comment(
            f"> {msg} {f'({structured_log})' if structured_log != '' else ''}"
        )

    return log_comment
