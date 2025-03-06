import os

from ebi_eva_common_pyutils.config import cfg


def load_config(*args):
    """Load a config file from any path provided.
    If none are provided then read from a file path provided in the environment variable RELEASE_CONFIG.
    If not provided then default to .release_config.yml place in the current users' home"""
    cfg.load_config_file(
        *args,
        os.getenv('RELEASE_CONFIG'),
        os.path.expanduser('~/.release_config.yml')
    )
