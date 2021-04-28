# from https://stackoverflow.com/questions/50877387/
# how-to-define-a-logger-in-python-once-for-the-whole-program

import logging
import socket
import datetime


def get_timestamp():
    """Simple function to return current time"""
    now = datetime.datetime.now()
    return now.strftime("%Y%m%d_%H%M%S")


def setup_log(output_dir="./", prefix="", suffix="", quiet=False):
    """Helper function to set up logs

    This method sets up a log based on the level of verbosity set
    by the user. This log should be accessible across all modules.

    Parameters
    ----------
    output_dir: string
        path where the files should be written
    prefix: string
        prefix to prepend to output info files
    suffix: string
        prefix to append to output info files
    quiet: bool
        whether to suppress logs by setting the logging level to ERROR

    Returns
    -------
    None

    """

    # log file name will start with a user prefix if provided, otherwise date
    if prefix == "":
        log_prefix = get_timestamp()
    else:
        log_prefix = prefix

    log_file = output_dir + log_prefix + suffix + ".log"

    # setup logs
    if quiet:
        # only alert on errors and do NOT write out to file
        setup_logging("", "ERROR")
    else:
        # set logging to info and write to file
        # currently no option for warnings only
        setup_logging(log_file, "INFO")


def setup_logging(log_file="", level="DEBUG"):
    """Adds a configured stream handler to the root logger

    Parameters
    ----------
    output_dir: string
        path where the files should be written
    level: string
        logging level, must be "ERROR","WARNING","INFO", or "DEBUG"

    Returns
    -------
    None

    """
    syslog = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s : %(levelname)s : %(name)s : %(message)s"
    )
    syslog.setFormatter(formatter)

    logger = logging.getLogger()
    logger.addHandler(syslog)

    if log_file != "":
        fh = logging.FileHandler(log_file)
        logger.addHandler(fh)

    if level == "DEBUG":
        logger.setLevel(logging.DEBUG)
    elif level == "INFO":
        logger.setLevel(logging.INFO)
    elif level == "WARNING":
        logger.setLevel(logging.WARNING)
    elif level == "ERROR":
        logger.setLevel(logging.ERROR)
    else:
        raise Exception(
            "Logging level:" + str(level) + " not recognized.\n"
            'Must be one of "ERROR","WARNING","INFO", or "DEBUG"'
        )


def host_log_adapter(logger):
    """Helper function to retrieve logger"""
    hostname = {"hostname": socket.gethostname()}
    return logging.LoggerAdapter(logger, hostname)


# creates an object to enable other modules to retrieve from other modules
# may need to evaluate if this should be refactored to make sure the desired
# log is setup before retrieving
logger = host_log_adapter(logging.getLogger())
